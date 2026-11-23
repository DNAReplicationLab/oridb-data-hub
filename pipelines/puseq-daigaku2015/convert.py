#!/usr/bin/env python3
"""Convert the Daigaku 2015 Pu-seq profiles for the OriDB hub.

    python3 convert.py --hub-root /path/to/oridb-data-hub [--raw-dir DIR]

GSE62108 publishes four derived genome-wide profiles at series level, as wiggle and
bedGraph on ASM294v2. This renames their sequences, fixes their edges and converts them
to bigWig, and turns the origin-efficiency bedGraph into a bigBed call set as well as a
profile, because a tier rule needs discrete calls and not a signal.

Three things need doing to the files as published, and each is recorded in dataset.json:

  - **Sequence names.** The wigs use `chr1`, `chr2`, `chr3`; this hub uses PomBase's
    `I`, `II`, `III`, which is also what OriDB's pombe records use.
  - **Edges.** The last bin of a chromosome is written at the full span even where the
    chromosome ends inside it, so a few intervals run past the end of the sequence.
    Those are clipped. One interval in the leftward-fork file starts past the end of
    chromosome I altogether and is dropped: there is nothing to clip it to.
  - **Nothing else.** Values are published as they are - no rescaling, no smoothing.

Requires bedGraphToBigWig and bedToBigBed (UCSC) on PATH.
"""
import argparse
import gzip
import json
import os
import subprocess
import sys
import time
import urllib.error
import urllib.request

GEO = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE62nnn/GSE62108/suppl/'
RENAME = {'chr1': 'I', 'chr2': 'II', 'chr3': 'III'}

# Each profile, with how it should be drawn. viewLimits are set from the published
# range rather than from a guess, and the comment says what the numbers are.
PROFILES = {
    'trep': {
        'file': 'GSE62108_PU-seq_Trep.wig.gz',
        'measurement': 'replication timing, Trep in minutes after release',
        'out': 'puseq.trep.bw',
    },
    'fork_direction': {
        'file': 'GSE62108_PU-seq_leftward_moving_fork.wig.gz',
        'measurement': 'fraction of forks moving leftward, 0 to 1',
        'out': 'puseq.leftward_fork.bw',
    },
    'termination': {
        'file': 'GSE62108_PU-seq_terminaton-events.wig.gz',
        'measurement': 'percentage of termination events per bin',
        'out': 'puseq.termination.bw',
        'recentre': True,
    },
    # Origin efficiency is served ONLY as a call set, not as a profile. The published
    # bedGraph is 1,145 discrete 300 bp origin calls rather than a continuous signal,
    # and six positions are listed twice with different efficiencies, which a bigWig
    # cannot hold (it forbids overlaps) and which must not be silently merged: both
    # values are what the study published. A bigBed keeps every record as it stands.
    'origin_efficiency': {
        'file': 'GSE62108_origin-effici.bedgraph.gz',
        'measurement': 'origin efficiency, percentage of cells firing the origin',
        'calls': 'puseq.origins.bb',
    },
}

ORIGINS_AS = """table puseqOrigins
"Daigaku 2015 Pu-seq origin calls on ASM294v2"
    (
    string  chrom;       "Reference sequence"
    uint    chromStart;  "Start, 0-based"
    uint    chromEnd;    "End, exclusive"
    string  name;        "Call identifier, chromosome and position"
    uint    score;       "Efficiency scaled to 0-1000 over the published range"
    char[1] strand;      "Not applicable"
    string  efficiency;  "Published efficiency, percentage of cells"
    )
"""


# NCBI returns 403 to a burst of requests rather than 429, so a run that fetches
# several files in a row needs to back off and retry rather than fail.
UA = 'oridb-data-hub/1.0 (+https://github.com/DNAReplicationLab/oridb-data-hub)'


def fetch(name, raw_dir, attempts=4):
    path = os.path.join(raw_dir, name)
    if os.path.exists(path):
        return path
    os.makedirs(raw_dir, exist_ok=True)
    request = urllib.request.Request(GEO + name, headers={'User-Agent': UA})
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(request, timeout=300) as r, open(path, 'wb') as out:
                out.write(r.read())
            return path
        except urllib.error.HTTPError as err:
            if err.code not in (403, 429, 500, 503) or attempt == attempts - 1:
                raise
            time.sleep(5 * (attempt + 1))
    return path


def read_sizes(path):
    sizes = {}
    for line in open(path):
        c, n = line.split()[:2]
        sizes[c] = int(n)
    return sizes


def wig_intervals(path):
    """variableStep wiggle -> (chrom, start0, end, value). Positions are 1-based."""
    chrom, span = None, 1
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('track'):
                continue
            if line.startswith('variableStep'):
                fields = dict(f.split('=') for f in line.split()[1:])
                chrom = RENAME[fields['chrom']]
                span = int(fields.get('span', 1))
                continue
            pos, value = line.split()
            start = int(pos) - 1
            yield chrom, start, start + span, float(value)


def bedgraph_intervals(path):
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            chrom, start, end, value = line.split()
            yield RENAME[chrom], int(start), int(end), float(value)


def recentre(intervals):
    """A sliding window cannot be a bigWig, so give each value its own bin.

    The termination profile is published as a 1 kb window stepped every 300 bp, so
    consecutive intervals overlap by 700 bp and bigWig refuses them. Each value is
    therefore written to a 300 bp bin centred on its window: no value is changed and
    the sampling grid is preserved, but the bins tile instead of overlapping. The
    width of the window the value was computed over is recorded in dataset.json.
    """
    rows = list(intervals)
    steps = sorted(b[1] - a[1] for a, b in zip(rows, rows[1:])
                   if a[0] == b[0] and b[1] > a[1])
    step = steps[len(steps) // 2]
    out = []
    for chrom, start, end, value in rows:
        centre = (start + end) // 2
        out.append((chrom, max(0, centre - step // 2), centre - step // 2 + step, value))
    return out, step, rows[0][2] - rows[0][1]


def write_bedgraph(intervals, sizes, out_path):
    """Clip to the chromosome, drop what starts past it, and report both."""
    clipped = dropped = kept = 0
    rows = []
    for chrom, start, end, value in intervals:
        limit = sizes[chrom]
        if start >= limit:
            dropped += 1
            continue
        if end > limit:
            end = limit
            clipped += 1
        rows.append((chrom, start, end, value))
        kept += 1
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(out_path, 'w') as fh:
        for chrom, start, end, value in rows:
            fh.write(f'{chrom}\t{start}\t{end}\t{value:g}\n')
    return {'intervals': kept, 'clipped_to_chromosome_end': clipped,
            'dropped_past_chromosome_end': dropped}, rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--hub-root', required=True)
    ap.add_argument('--raw-dir', default='raw', help='where the GEO files are cached')
    ap.add_argument('--genome', default='ASM294v2')
    a = ap.parse_args()

    sizes_path = os.path.join(a.hub_root, f'{a.genome}.chrom.sizes')
    if not os.path.exists(sizes_path):
        sys.exit(f'no {sizes_path}; the hub must serve {a.genome} before a track can be added')
    sizes = read_sizes(sizes_path)

    dest = os.path.join(a.hub_root, a.genome, 'daigaku2015')
    os.makedirs(dest, exist_ok=True)
    work = os.path.join(a.raw_dir, 'work')
    os.makedirs(work, exist_ok=True)

    record = {}
    for key, spec in PROFILES.items():
        src = fetch(spec['file'], a.raw_dir)
        reader = bedgraph_intervals if spec['file'].endswith('.bedgraph.gz') else wig_intervals
        source = reader(src)
        window = None
        if spec.get('recentre'):
            source, step, window = recentre(source)
        bg = os.path.join(work, key + '.bedgraph')
        stats, rows = write_bedgraph(source, sizes, bg)
        if window is not None:
            stats['published_window_bp'] = window
            stats['window_step_bp'] = step
        values = [r[3] for r in rows]
        stats.update({
            'source_file': spec['file'],
            'measurement': spec['measurement'],
            'bin_size_bp': rows[0][2] - rows[0][1],
            'value_min': round(min(values), 4),
            'value_max': round(max(values), 4),
        })
        if 'out' in spec:
            bw = os.path.join(dest, spec['out'])
            subprocess.run(['bedGraphToBigWig', bg, sizes_path, bw], check=True)
            stats['track'] = spec['out']
            stats['bytes'] = os.path.getsize(bw)

        # The efficiency profile is also the study's origin call set: each interval is
        # one called origin, not a bin of a continuous signal.
        if 'calls' in spec:
            lo, hi = min(values), max(values)
            bed = os.path.join(work, 'origins.bed')
            with open(bed, 'w') as fh:
                for chrom, start, end, value in rows:
                    score = int(round(1000 * (value - lo) / (hi - lo))) if hi > lo else 0
                    fh.write(f'{chrom}\t{start}\t{end}\t{chrom}-{start + 1}\t{score}\t.\t{value:g}\n')
            asf = os.path.join(dest, 'puseq.origins.as')
            with open(asf, 'w') as fh:
                fh.write(ORIGINS_AS)
            bb = os.path.join(dest, spec['calls'])
            subprocess.run(['bedToBigBed', f'-as={asf}', '-type=bed6+1', '-tab',
                            bed, sizes_path, bb], check=True)
            repeated = len(rows) - len({(r[0], r[1]) for r in rows})
            stats['call_set'] = {'track': spec['calls'], 'calls': len(rows),
                                 'positions_listed_twice': repeated,
                                 'bytes': os.path.getsize(bb)}
        record[key] = stats
        print(f"{key}: {stats['intervals']} intervals "
              f"({stats['clipped_to_chromosome_end']} clipped, "
              f"{stats['dropped_past_chromosome_end']} dropped) -> "
              f"{spec.get('out') or spec['calls']}")

    dataset = {
        'dataset_id': 'daigaku2015-asm294v2',
        'title': 'Polymerase usage sequencing (Pu-seq) in S. pombe (Daigaku et al. 2015), '
                 'the published derived profiles on ASM294v2',
        'organism': 'Schizosaccharomyces pombe',
        'assembly': a.genome,
        'bin_size_bp': record['trep']['bin_size_bp'],
        'measurement': 'replication timing, fork direction, termination and origin efficiency',
        'scale': 'see each track; the units differ between them',
        'source': {
            'doi': '10.1038/nsmb.2962',
            'pubmed': '25664722',
            'archive_series': ['GSE62108'],
            'files_used': [p['file'] for p in PROFILES.values()],
            'note': 'Series-level derived profiles. The 25 per-sample strand-resolved '
                    'wigs in GSE62108_RAW.tar are not served.',
        },
        'provenance': {
            'pipeline': 'pipelines/puseq-daigaku2015/convert.py',
            'pipeline_revision': None,
            'transformations': [
                'sequences renamed chr1/chr2/chr3 -> I/II/III',
                'intervals clipped to the chromosome end; intervals starting past it dropped',
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
