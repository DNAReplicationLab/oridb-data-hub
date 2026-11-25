#!/usr/bin/env python3
"""Add the Reuter 2024 ORC/MCM ChIP-exo tracks to the OriDB hub.

    python3 convert.py --hub-root /path/to/oridb-data-hub [--raw-dir DIR]

GSE240779 is a SuperSeries; the G1- and G2-arrest ChIP-exo samples are in its
subseries GSE173527 (Reuter et al. 2024, Nat Commun 15:7306). Six samples are taken:
Mcm4-Flag in G1, Orc2-Flag in G1 and Orc2-Flag in G2, two biological replicates each.
The G1-release time course in the other subseries (GSE240472) is not included.

Unlike the other pipelines here, this one does not re-derive coverage from reads. The
submitters' bigWigs are already on sacCer3 - all 17 sequence names and lengths match
the hub's chrom.sizes exactly - so they are served as deposited, and the pipeline's job
is to verify that, name the files for the hub, and record provenance.

Two file flavours are deposited per sample and they are NOT on the same scale:

  - the combined `*wt_G1_1.bw` track is CPM-normalised with reads extended to about
    209 bp (deepTools --normalizeUsing CPM); every one of the six has a genome-wide
    mean of 17.1, so they are directly comparable. These are copied byte-for-byte.
  - the strand-split `*_pos.bw` / `*_neg.bw` pair is RAW integer read coverage at the
    50 bp read length, signed by strand, with no depth normalisation at all. Library
    depths across these six samples differ about 1.5-fold, so raw strand tracks would
    put a depth difference on screen as an occupancy difference.

So the strand tracks are the one thing this pipeline changes: each is multiplied by
1e6 / (mapped reads in that SAMPLE), a single factor per sample applied to both
strands. Using a per-file factor would have rescaled the two strands independently and
destroyed the strand asymmetry, which is the whole reason for serving them. Mapped
reads are recovered as (sum of |coverage|) / 50 over both strands; the factor and the
read count for every sample are recorded in dataset.json, so the deposited integers can
be recovered by dividing back. After scaling, the strand tracks are per million reads
like the combined track, but about 4.2x lower in absolute height because they carry a
50 bp read footprint rather than a 209 bp extended one; they therefore get their own
viewLimits in the container spec rather than sharing the combined track's.

Between-replicate Pearson correlation is computed for each condition from the combined
tracks and recorded in dataset.json as the quality statement for the dataset. No
mean-of-replicates track is written: both replicates are served and the reader compares
them.

Requires pyBigWig and numpy. Writes <hub-root>/sacCer3/reuter2024/ and refuses to run
unless sacCer3.chrom.sizes is already there.
"""
import argparse
import datetime
import hashlib
import json
import os
import shutil
import sys
import time
import urllib.error
import urllib.request

import numpy as np
import pyBigWig

GEO = 'https://ftp.ncbi.nlm.nih.gov/geo/samples/'
UA = 'oridb-data-hub/1.0 (+https://github.com/DNAReplicationLab/oridb-data-hub)'
READ_LEN = 50          # single-end read length, from the series' data processing fields
GENOME = 'sacCer3'
DATASET = 'reuter2024'

# One entry per sample: the GEO accession, the condition this pipeline calls it, and the
# three deposited bigWig basenames. The strand files carry the submitters' own flowcell
# lane suffixes, which is why they cannot be derived from the accession.
SAMPLES = [
    {'gsm': 'GSM5269283', 'out': 'mcm4_g1_rep1',
     'title': 'Mcm4-Flag ChIP-Exo from cells arrested in G1, rep1',
     'combined': 'GSM5269283_Mcm4wt_G1_1.bw',
     'pos': 'GSM5269283_MCM4_G1_1_S39_L005_R1_001_pos.bw',
     'neg': 'GSM5269283_MCM4_G1_1_S39_L005_R1_001_neg.bw'},
    {'gsm': 'GSM5269284', 'out': 'mcm4_g1_rep2',
     'title': 'Mcm4-Flag ChIP-Exo from cells arrested in G1, rep2',
     'combined': 'GSM5269284_Mcm4wt_G1_2.bw',
     'pos': 'GSM5269284_MCM4_G1_2_S40_L005_R1_001_pos.bw',
     'neg': 'GSM5269284_MCM4_G1_2_S40_L005_R1_001_neg.bw'},
    {'gsm': 'GSM5269285', 'out': 'orc2_g1_rep1',
     'title': 'Orc2-Flag ChIP-Exo from cells arrested in G1, rep1',
     'combined': 'GSM5269285_Orc2wt_G1_1.bw',
     'pos': 'GSM5269285_ORC2_G1_1_S43_L005_R1_001_pos.bw',
     'neg': 'GSM5269285_ORC2_G1_1_S43_L005_R1_001_neg.bw'},
    {'gsm': 'GSM5269286', 'out': 'orc2_g1_rep2',
     'title': 'Orc2-Flag ChIP-Exo from cells arrested in G1, rep2',
     'combined': 'GSM5269286_Orc2wt_G1_2.bw',
     'pos': 'GSM5269286_ORC2_G1_2_S44_L005_R1_001_pos.bw',
     'neg': 'GSM5269286_ORC2_G1_2_S44_L005_R1_001_neg.bw'},
    {'gsm': 'GSM5269287', 'out': 'orc2_g2_rep1',
     'title': 'Orc2-Flag ChIP-Exo from cells arrested in G2, rep1',
     'combined': 'GSM5269287_Orc2wt_G2_1.bw',
     'pos': 'GSM5269287_ORC2_G2_1_S47_L005_R1_001_pos.bw',
     'neg': 'GSM5269287_ORC2_G2_1_S47_L005_R1_001_neg.bw'},
    {'gsm': 'GSM5269288', 'out': 'orc2_g2_rep2',
     'title': 'Orc2-Flag ChIP-Exo from cells arrested in G2, rep2',
     'combined': 'GSM5269288_Orc2wt_G2_2.bw',
     'pos': 'GSM5269288_ORC2_G2_2_S48_L005_R1_001_pos.bw',
     'neg': 'GSM5269288_ORC2_G2_2_S48_L005_R1_001_neg.bw'},
]

CONDITIONS = [('mcm4_g1', 'Mcm4-Flag, G1 arrest'),
              ('orc2_g1', 'Orc2-Flag, G1 arrest'),
              ('orc2_g2', 'Orc2-Flag, G2 arrest')]


def fetch(gsm, basename, raw_dir):
    """Download one supplementary file, with the back-off GEO's FTP host needs.

    GEO answers a burst of requests with 403 rather than 429, so a plain retry loop that
    only handles 429 gives up on the first rate limit.
    """
    dest = os.path.join(raw_dir, basename)
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    url = f'{GEO}{gsm[:-3]}nnn/{gsm}/suppl/{basename}'
    for attempt in range(6):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': UA})
            with urllib.request.urlopen(req, timeout=180) as r, open(dest + '.part', 'wb') as f:
                shutil.copyfileobj(r, f)
            os.replace(dest + '.part', dest)
            return dest
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as e:
            code = getattr(e, 'code', None)
            if code is not None and code not in (403, 429, 500, 502, 503):
                raise
            wait = 5 * 2 ** attempt
            print(f'  {basename}: {e}; retrying in {wait}s', file=sys.stderr)
            time.sleep(wait)
    sys.exit(f'{basename}: still failing after 6 attempts')


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def check_assembly(path, sizes):
    """Refuse any file whose sequences are not exactly the hub's assembly."""
    bw = pyBigWig.open(path)
    chroms = bw.chroms()
    bw.close()
    if chroms != sizes:
        extra = set(chroms) - set(sizes)
        missing = set(sizes) - set(chroms)
        bad = {k: (chroms[k], sizes[k]) for k in set(chroms) & set(sizes)
               if chroms[k] != sizes[k]}
        sys.exit(f'{os.path.basename(path)} is not on {GENOME}: '
                 f'extra={sorted(extra)} missing={sorted(missing)} lengths={bad}')


def total_abs_coverage(path):
    bw = pyBigWig.open(path)
    total = abs(bw.header()['sumData'])
    bw.close()
    return total


def write_scaled(src, dest, factor, sizes, order):
    """Copy a bigWig multiplying every value by `factor`, interval structure kept."""
    src_bw = pyBigWig.open(src)
    out = pyBigWig.open(dest, 'w')
    out.addHeader([(c, sizes[c]) for c in order])
    for chrom in order:
        iv = src_bw.intervals(chrom, 0, sizes[chrom])
        if not iv:
            continue
        starts = np.fromiter((i[0] for i in iv), dtype=np.int64, count=len(iv))
        ends = np.fromiter((i[1] for i in iv), dtype=np.int64, count=len(iv))
        values = np.fromiter((i[2] for i in iv), dtype=np.float64, count=len(iv)) * factor
        out.addEntries([chrom] * len(iv), starts.tolist(), ends=ends.tolist(),
                       values=values.astype(np.float32).tolist())
    out.close()
    src_bw.close()


def whole_genome(path, order, sizes):
    bw = pyBigWig.open(path)
    v = np.concatenate([np.nan_to_num(bw.values(c, 0, sizes[c], numpy=True)) for c in order])
    bw.close()
    return v


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--hub-root', required=True, help='repository root holding hub.txt')
    ap.add_argument('--raw-dir', default=None,
                    help='cache for the downloaded GEO files (default <hub-root>/../raw/reuter2024)')
    args = ap.parse_args()

    chrom_sizes = os.path.join(args.hub_root, f'{GENOME}.chrom.sizes')
    if not os.path.exists(chrom_sizes):
        sys.exit(f'{chrom_sizes} not found: add the {GENOME} assembly to the hub first')
    order, sizes = [], {}
    for line in open(chrom_sizes):
        if line.strip():
            name, length = line.split()
            order.append(name)
            sizes[name] = int(length)

    raw_dir = args.raw_dir or os.path.join(args.hub_root, os.pardir, 'raw', DATASET)
    os.makedirs(raw_dir, exist_ok=True)
    out_dir = os.path.join(args.hub_root, GENOME, DATASET)
    os.makedirs(out_dir, exist_ok=True)

    records = []
    for s in SAMPLES:
        print(f'{s["gsm"]}  {s["out"]}')
        paths = {k: fetch(s['gsm'], s[k], raw_dir) for k in ('combined', 'pos', 'neg')}
        for p in paths.values():
            check_assembly(p, sizes)

        # One scale factor per SAMPLE, from both strands together, so the pos:neg ratio
        # survives. Coverage of a 50 bp read contributes 50 to the sum of |coverage|.
        reads = (total_abs_coverage(paths['pos']) + total_abs_coverage(paths['neg'])) / READ_LEN
        factor = 1e6 / reads

        combined_out = os.path.join(out_dir, f'chipexo.{s["out"]}.bw')
        shutil.copyfile(paths['combined'], combined_out)
        for strand in ('pos', 'neg'):
            write_scaled(paths[strand], os.path.join(out_dir, f'chipexo.{s["out"]}.{strand}.bw'),
                         factor, sizes, order)

        records.append({
            'sample': s['gsm'],
            'title': s['title'],
            'track': s['out'],
            'combined_file': f'chipexo.{s["out"]}.bw',
            'combined_source': s['combined'],
            'combined_sha256': sha256(paths['combined']),
            'combined_transformation': 'none; copied as deposited',
            'strand_files': [f'chipexo.{s["out"]}.pos.bw', f'chipexo.{s["out"]}.neg.bw'],
            'strand_sources': [s['pos'], s['neg']],
            'strand_sha256': [sha256(paths['pos']), sha256(paths['neg'])],
            'mapped_reads': round(reads),
            'strand_scale_factor': round(factor, 6),
            'strand_transformation': (
                f'raw read coverage multiplied by {factor:.6f} = 1e6 / {round(reads)} mapped '
                f'reads (both strands, from sum|coverage|/{READ_LEN}); divide by this factor '
                f'to recover the deposited integers'),
        })

    # Between-replicate agreement, from the combined tracks at base resolution.
    correlations = {}
    for cond, label in CONDITIONS:
        a = whole_genome(os.path.join(out_dir, f'chipexo.{cond}_rep1.bw'), order, sizes)
        b = whole_genome(os.path.join(out_dir, f'chipexo.{cond}_rep2.bw'), order, sizes)
        r = float(np.corrcoef(a, b)[0, 1])
        correlations[cond] = round(r, 4)
        print(f'  {label}: rep1 vs rep2 Pearson r = {r:.4f}')

    dataset = {
        'dataset_id': f'{DATASET}-saccer3',
        'title': ('ORC and MCM2-7 ChIP-exo in G1 and G2 arrest '
                  '(Reuter et al. 2024), as deposited on sacCer3'),
        'organism': 'Saccharomyces cerevisiae',
        'assembly': GENOME,
        'bin_size_bp': 1,
        'measurement': 'Orc2 and Mcm4 occupancy by ChIP-exo',
        'scale': ('combined tracks: CPM with reads extended to ~209 bp, genome-wide mean '
                  '17.1 in every sample; strand tracks: coverage per million mapped reads '
                  'at the 50 bp read length, about 4.2x lower in height for that reason'),
        'source': {
            'doi': '10.1038/s41467-024-51538-9',
            'pubmed': '39181881',
            'citation': ('Reuter LM, Khadayate SP, Mossler A, Liebl K, Faull SV, Karimi MM, '
                         'Speck C. MCM2-7 loading-dependent ORC release ensures genome-wide '
                         'origin licensing. Nat Commun 15:7306 (2024)'),
            'archive_series': ['GSE240779', 'GSE173527'],
            'bioproject': 'PRJNA1005236',
            'samples_used': [f'{s["gsm"]} ({s["title"]})' for s in SAMPLES],
            'samples_excluded': ['GSE240472 subseries (Orc2 ChIP-exo 24 min after G1 release)',
                                 'ChExMix peak calls deposited beside each bigWig'],
            'exclusion_reason': ('the hub carries the G1 and G2 arrest conditions; the release '
                                 'time course and the per-replicate peak calls were not '
                                 'requested and the peak calls are unfiltered'),
        },
        'provenance': {
            'pipeline': f'pipelines/chipexo-{DATASET}',
            'generated': datetime.datetime.now(datetime.timezone.utc)
                                 .strftime('%Y-%m-%dT%H:%M:%SZ'),
            'submitter_processing': [
                'Bowtie2 v2.3.5 to sacCer3, 50 bp single-end reads, Illumina HiSeq 2500',
                'deepTools v3.2.1 bigWig generation, --normalizeUsing CPM (combined tracks only)',
            ],
            'validation': (
                f'all 18 deposited bigWigs carry exactly the {len(order)} sequences of '
                f'{GENOME}.chrom.sizes at identical lengths, so no renaming or lift was '
                'needed; combined tracks were copied byte-for-byte and their sha256 of the '
                'deposited file is recorded per sample'),
            'replicate_correlation': {
                'statistic': 'Pearson r between replicates, combined tracks, base resolution, '
                             'whole genome',
                'values': correlations,
            },
            'caveat': (
                'the deposited combined track is not the sum of the deposited strand pair '
                '(r = 0.83-0.86 between them): the combined track extends reads to ~209 bp '
                'while the strand tracks keep the 50 bp read footprint, so the combined '
                'track is the smoother of the two views of the same reads'),
        },
        'samples': records,
    }
    with open(os.path.join(out_dir, 'dataset.json'), 'w') as f:
        json.dump(dataset, f, indent=2)
        f.write('\n')
    print(f'wrote {out_dir}')


if __name__ == '__main__':
    main()
