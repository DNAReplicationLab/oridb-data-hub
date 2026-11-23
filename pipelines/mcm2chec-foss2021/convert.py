#!/usr/bin/env python3
"""Convert the Foss 2021 Mcm2-ChEC pombe read depths for the OriDB hub.

    python3 convert.py --hub-root /path/to/oridb-data-hub [--raw-dir DIR] [--bin 25]

GSE150800 covers five species; its *S. pombe* arm is deposited as per-base read-depth
tables - three files, 12.57 million rows each - rather than as bedGraph, so there is
nothing to serve until they are binned and converted.

Two decisions are made here and recorded in dataset.json:

  - **Binning.** Depth is averaged into 25 bp bins. At one base per interval a single
    track is 38 MB of bigWig; at 25 bp it is 3.6 MB, and 25 bp is still an order of
    magnitude finer than the 300 bp - 1 kb bins of the Pu-seq tracks it sits beside.
    ChEC's own resolution is finer than this: for footprint-level work, go to the
    archive files, which this pipeline will reproduce at --bin 1.
  - **Scale.** Each track is divided by its own genome-wide mean depth, so a value of
    1 is average coverage for that library. The libraries differ several-fold in depth,
    and a container that overlays HU against G1 is only readable if the two are on the
    same scale. The raw mean depth of every column is recorded, so the published counts
    can be recovered by multiplying back.

The free-MNase control is converted the same way and served beside the samples: ChEC
signal is read against it, and a hub that served only the samples would invite the
comparison to be skipped.

Requires bedGraphToBigWig (UCSC) on PATH, plus numpy and pandas.
"""
import argparse
import json
import os
import subprocess
import sys
import time
import urllib.error
import urllib.request

import numpy as np
import pandas as pd

GEO = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150800/suppl/'
UA = 'oridb-data-hub/1.0 (+https://github.com/DNAReplicationLab/oridb-data-hub)'
RENAME = {'chrI': 'I', 'chrII': 'II', 'chrIII': 'III'}

# One entry per output track: which archive file it comes from, which column, and what
# it is. The column names are the submitter's own, from the file's header.
TRACKS = [
    {'out': 'mcm2chec.hu_rep1.all.bw',
     'file': 'GSE150800_pombe_MCM2ChEC_HU_G1_read_depths_3_size_ranges_081121_1.txt.gz',
     'column': 'pombe_MCM2ChEC_HU_rep1_library_fragment_size_range_all',
     'label': 'Mcm2-ChEC HU rep1, all fragment sizes'},
    {'out': 'mcm2chec.g1_rep1.all.bw',
     'file': 'GSE150800_pombe_MCM2ChEC_HU_G1_read_depths_3_size_ranges_081121_1.txt.gz',
     'column': 'pombe_MCM2ChEC_G1_rep1_library_fragment_size_range_all',
     'label': 'Mcm2-ChEC G1 rep1, all fragment sizes'},
    {'out': 'mcm2chec.hu_rep1.51_100bp.bw',
     'file': 'GSE150800_pombe_MCM2ChEC_HU_G1_read_depths_3_size_ranges_081121_1.txt.gz',
     'column': 'pombe_MCM2ChEC_HU_rep1_library_fragment_size_range_51bpto100bp',
     'label': 'Mcm2-ChEC HU rep1, 51-100 bp fragments'},
    {'out': 'mcm2chec.g1_rep1.51_100bp.bw',
     'file': 'GSE150800_pombe_MCM2ChEC_HU_G1_read_depths_3_size_ranges_081121_1.txt.gz',
     'column': 'pombe_MCM2ChEC_G1_rep1_library_fragment_size_range_51bpto100bp',
     'label': 'Mcm2-ChEC G1 rep1, 51-100 bp fragments'},
    {'out': 'mcm2chec.hu_rep1.151_200bp.bw',
     'file': 'GSE150800_pombe_MCM2ChEC_HU_G1_read_depths_3_size_ranges_081121_1.txt.gz',
     'column': 'pombe_MCM2ChEC_HU_rep1_library_fragment_size_range_151bpto200bp',
     'label': 'Mcm2-ChEC HU rep1, 151-200 bp fragments'},
    {'out': 'mcm2chec.g1_rep1.151_200bp.bw',
     'file': 'GSE150800_pombe_MCM2ChEC_HU_G1_read_depths_3_size_ranges_081121_1.txt.gz',
     'column': 'pombe_MCM2ChEC_G1_rep1_library_fragment_size_range_151bpto200bp',
     'label': 'Mcm2-ChEC G1 rep1, 151-200 bp fragments'},
    {'out': 'mcm2chec.hu_rep2.all.bw',
     'file': 'GSE150800_processed_pombe_MCM2ChEC_HU_rep2_080621_1.txt.gz',
     'column': None,
     'label': 'Mcm2-ChEC HU rep2'},
    {'out': 'free_mnase.rep1.bw',
     'file': 'GSE150800_processed_pombe_free_MNase_rep1_081221_1.txt.gz',
     'column': None,
     'label': 'Free MNase control rep1'},
]


def fetch(name, raw_dir, attempts=4):
    """NCBI answers a burst of requests with 403, so back off rather than fail."""
    path = os.path.join(raw_dir, name)
    if os.path.exists(path):
        return path
    os.makedirs(raw_dir, exist_ok=True)
    request = urllib.request.Request(GEO + name, headers={'User-Agent': UA})
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(request, timeout=900) as r, open(path, 'wb') as out:
                out.write(r.read())
            return path
        except urllib.error.HTTPError as err:
            if err.code not in (403, 429, 500, 503) or attempt == attempts - 1:
                raise
            time.sleep(5 * (attempt + 1))
    return path


def read_sizes(path):
    return {c: int(n) for c, n in (l.split()[:2] for l in open(path) if l.strip())}


def load(path, columns):
    """Per-base depth table -> {column: {chrom: array indexed by 1-based coordinate}}."""
    has_header = 'read_depths_3_size_ranges' in os.path.basename(path)
    if has_header:
        frame = pd.read_csv(path, sep='\t')
        frame = frame.rename(columns={frame.columns[0]: 'chr', frame.columns[1]: 'coord'})
    else:
        frame = pd.read_csv(path, sep='\t', header=None, names=['chr', 'coord', columns[0]])
    frame['chr'] = frame['chr'].map(RENAME)
    frame = frame[frame['chr'].notna()]
    return frame


def to_bigwig(frame, column, sizes, bin_size, bedgraph, bigwig, chrom_sizes_path):
    """Average depth into bins, scale to the genome-wide mean, and convert."""
    total = float(frame[column].sum())
    bases = sum(sizes.values())
    mean_depth = total / bases
    written = 0
    with open(bedgraph, 'w') as out:
        for chrom in sorted(sizes, key=lambda c: (len(c), c)):
            size = sizes[chrom]
            depth = np.zeros(size, dtype=np.float64)
            part = frame[frame['chr'] == chrom]
            coords = part['coord'].to_numpy() - 1          # the tables are 1-based
            values = part[column].to_numpy(dtype=np.float64)
            keep = (coords >= 0) & (coords < size)
            depth[coords[keep]] = values[keep]
            n = size // bin_size * bin_size
            binned = depth[:n].reshape(-1, bin_size).mean(axis=1) / mean_depth
            for i in np.flatnonzero(binned):
                out.write(f'{chrom}\t{i * bin_size}\t{i * bin_size + bin_size}\t{binned[i]:.4f}\n')
                written += 1
            tail = depth[n:]
            if tail.size and tail.mean():
                out.write(f'{chrom}\t{n}\t{size}\t{tail.mean() / mean_depth:.4f}\n')
                written += 1
    subprocess.run(['bedGraphToBigWig', bedgraph, chrom_sizes_path, bigwig], check=True)
    return {'bins': written, 'mean_depth_per_base': round(mean_depth, 4),
            'total_depth': total, 'bytes': os.path.getsize(bigwig)}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--hub-root', required=True)
    ap.add_argument('--raw-dir', default='raw')
    ap.add_argument('--genome', default='ASM294v2')
    ap.add_argument('--bin', type=int, default=25, help='bin size in bp')
    a = ap.parse_args()

    chrom_sizes_path = os.path.join(a.hub_root, f'{a.genome}.chrom.sizes')
    if not os.path.exists(chrom_sizes_path):
        sys.exit(f'no {chrom_sizes_path}; the hub must serve {a.genome} first')
    sizes = read_sizes(chrom_sizes_path)

    dest = os.path.join(a.hub_root, a.genome, 'foss2021')
    os.makedirs(dest, exist_ok=True)
    work = os.path.join(a.raw_dir, 'work')
    os.makedirs(work, exist_ok=True)

    record = {}
    frames = {}
    for spec in TRACKS:
        path = fetch(spec['file'], a.raw_dir)
        if spec['file'] not in frames:
            column = [spec['column'] or spec['out'].replace('.bw', '')]
            frames[spec['file']] = load(path, column)
        frame = frames[spec['file']]
        column = spec['column'] or frame.columns[2]
        stats = to_bigwig(frame, column, sizes, a.bin,
                          os.path.join(work, spec['out'].replace('.bw', '.bedgraph')),
                          os.path.join(dest, spec['out']), chrom_sizes_path)
        stats.update({'source_file': spec['file'], 'source_column': column,
                      'label': spec['label'], 'bin_size_bp': a.bin})
        record[spec['out']] = stats
        print(f"{spec['out']}: {stats['bins']} bins, mean depth "
              f"{stats['mean_depth_per_base']:.3f} per base, {stats['bytes'] // 1024} kB")

    dataset = {
        'dataset_id': 'foss2021-asm294v2',
        'title': 'Mcm2-ChEC in S. pombe (Foss et al. 2021), per-base read depths binned '
                 f'to {a.bin} bp on ASM294v2',
        'organism': 'Schizosaccharomyces pombe',
        'assembly': a.genome,
        'bin_size_bp': a.bin,
        'measurement': 'MCM2-7 occupancy by chromatin endogenous cleavage (ChEC)',
        'scale': 'read depth divided by the genome-wide mean depth of that library; '
                 '1 = average coverage',
        'source': {
            'doi': '10.1371/journal.pgen.1009714',
            'pubmed': '34473702',
            'archive_series': ['GSE150800'],
            'files_used': sorted({s['file'] for s in TRACKS}),
            'note': 'The pombe arm of a five-species series. The archive deposits '
                    'per-base read-depth tables for pombe, not bedGraph.',
        },
        'provenance': {
            'pipeline': 'pipelines/mcm2chec-foss2021/convert.py',
            'pipeline_revision': None,
            'transformations': [
                'sequences renamed chrI/chrII/chrIII -> I/II/III',
                f'per-base depth averaged into {a.bin} bp bins',
                'each track divided by its own genome-wide mean depth',
            ],
        },
        'tracks': record,
    }
    with open(os.path.join(dest, 'dataset.json'), 'w') as fh:
        json.dump(dataset, fh, indent=2)
        fh.write('\n')
    print(f'-> {dest}')


if __name__ == '__main__':
    main()
