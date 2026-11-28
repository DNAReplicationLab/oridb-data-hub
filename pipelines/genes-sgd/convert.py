#!/usr/bin/env python3
"""Build the SGD gene annotation tracks for the OriDB hub.

    python3 convert.py --hub-root /path/to/oridb-data-hub [--raw-dir DIR]

Two tracks on sacCer3, from SGD's own chromosomal-feature GFF3:

  - sgd.protein_coding.bb - all 6,613 ORFs, coloured by SGD's orf_classification
    (Verified, Uncharacterized, Dubious) in the same green/blue/grey the site uses for
    origin evidence tiers, so a Dubious ORF beside an origin is visibly dubious rather
    than quietly absent.
  - sgd.rna_genes.bb - all 438 non-coding RNA genes: 299 tRNA, 77 snoRNA, 31 ncRNA,
    24 rRNA, 6 snRNA and 1 telomerase RNA, coloured by class. One track rather than
    six because five of the classes hold 31 features or fewer, and colour tells them
    apart without costing five more rows in the panel. Transposable-element genes and
    pseudogenes are deliberately left out: neither is an RNA gene.

Why this is here at all. Every other dataset in the hub is a measurement; this one is
reference annotation, and UCSC's sacCer3 already carries SGD Genes natively. The
embedded JBrowse panel does not: it draws only what the hub declares, and an origin
without its neighbouring genes is hard to read. So the tracks live in an `annotation`
section of their own and are marked as panel defaults rather than joining the UCSC
overview, where they would be drawn on top of UCSC's own copy of the same annotation.

The assembly question, and why this is safe. SGD's current release is R64-4-1 while
sacCer3 is built on R64-1-1. Those differ in annotation, not in sequence: all 17
chromosome lengths are identical. Rather than trust that, the pipeline reads the
chromosome features out of the GFF and refuses to run unless every one matches
sacCer3.chrom.sizes exactly - if SGD ever does change the sequence, this stops instead
of silently writing features at wrong coordinates. Only `chrmt` is renamed, to `chrM`.

Gene models are real models, not spans: BED12 with blocks from the CDS features for an
ORF and from the noncoding_exon features for a tRNA, so the 378 annotated introns are
drawn as introns. Where SGD gives a gene two transcript records - 4,512 of the 6,613
do - the CDS features name both as parents and the block set is therefore the same, so
the union is taken rather than one isoform being picked arbitrarily.

Requires bedToBigBed (UCSC) on PATH. The conda build of it links libpng by its Homebrew
path; if it fails to start, run with DYLD_FALLBACK_LIBRARY_PATH=$CONDA_PREFIX/lib.
"""
import argparse
import datetime
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

GFF_URL = ('https://downloads.yeastgenome.org/curation/chromosomal_feature/'
           'saccharomyces_cerevisiae.gff.gz')
UA = 'oridb-data-hub/1.0 (+https://github.com/DNAReplicationLab/oridb-data-hub)'
GENOME = 'sacCer3'
DATASET = 'sgd'
RENAME = {'chrmt': 'chrM'}

# The site's evidence-tier palette, reused so that one colour means one thing across
# the hub: green for a call with support, blue for one that is probably real, grey for
# one that is not believed.
ORF_COLOUR = {'Verified': '26,127,75',
              'Uncharacterized': '47,125,155',
              'Dubious': '138,143,152'}
# One colour per RNA class, from a qualitative palette that stays distinguishable
# under deuteranopia. tRNA keeps the purple it carried when it was a track of its own.
# The second element of each pair is the gene feature's transcript child, which is
# where the exons hang: SGD writes gene -> transcript -> noncoding_exon throughout.
RNA_CLASSES = {
    'tRNA_gene': ('tRNA', '117,112,179'),
    'snoRNA_gene': ('snoRNA', '217,95,2'),
    'ncRNA_gene': ('ncRNA', '27,158,119'),
    'rRNA_gene': ('rRNA', '231,41,138'),
    'snRNA_gene': ('snRNA', '230,171,2'),
    'telomerase_RNA_gene': ('telomerase RNA', '102,102,102'),
}
RNA_TRANSCRIPTS = {'tRNA', 'snoRNA', 'ncRNA', 'rRNA', 'snRNA', 'telomerase_RNA'}

GENE_AS = """table sgdGene
"Protein-coding genes from the SGD chromosomal-feature GFF3"
    (
    string  chrom;              "Reference sequence chromosome"
    uint    chromStart;         "Start position in chromosome"
    uint    chromEnd;           "End position in chromosome"
    string  name;               "Systematic name, e.g. YFL039C"
    uint    score;              "Unused, always 0"
    char[1] strand;             "+ or -"
    uint    thickStart;         "Start of coding region"
    uint    thickEnd;           "End of coding region"
    uint    reserved;           "Colour, by ORF classification"
    int     blockCount;         "Number of coding blocks"
    int[blockCount] blockSizes; "Block sizes"
    int[blockCount] chromStarts; "Block starts relative to chromStart"
    string  geneName;           "Common gene name, e.g. ACT1, or empty"
    string  orfClass;           "SGD orf_classification"
    string  sgdid;              "SGD accession"
    string  description;        "SGD display description"
    )
"""

TRNA_AS = """table sgdRnaGene
"Non-coding RNA genes from the SGD chromosomal-feature GFF3"
    (
    string  chrom;              "Reference sequence chromosome"
    uint    chromStart;         "Start position in chromosome"
    uint    chromEnd;           "End position in chromosome"
    string  name;               "Systematic name, e.g. YNCA0002W"
    uint    score;              "Unused, always 0"
    char[1] strand;             "+ or -"
    uint    thickStart;         "Equal to chromStart: an RNA gene is not translated"
    uint    thickEnd;           "Equal to chromStart"
    uint    reserved;           "Colour, by RNA class"
    int     blockCount;         "Number of exons"
    int[blockCount] blockSizes; "Block sizes"
    int[blockCount] chromStarts; "Block starts relative to chromStart"
    string  geneName;           "Common gene name, e.g. TRN1, or empty"
    string  rnaClass;           "tRNA, snoRNA, ncRNA, rRNA, snRNA or telomerase RNA"
    string  sgdid;              "SGD accession"
    string  description;        "SGD display description"
    )
"""


def fetch(url, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': UA})
            with urllib.request.urlopen(req, timeout=300) as r, open(dest + '.part', 'wb') as f:
                shutil.copyfileobj(r, f)
            os.replace(dest + '.part', dest)
            return dest
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as e:
            wait = 5 * 2 ** attempt
            print(f'  {url}: {e}; retrying in {wait}s', file=sys.stderr)
            time.sleep(wait)
    sys.exit(f'{url}: still failing after 5 attempts')


def attrs(field):
    """GFF3 attribute column to a dict, percent-decoded."""
    out = {}
    for kv in field.split(';'):
        if '=' in kv:
            k, v = kv.split('=', 1)
            out[k] = urllib.parse.unquote(v)
    return out


def read_gff(path):
    rows = []
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue            # the GFF carries a FASTA section after the features
            rows.append((f[0], f[2], int(f[3]), int(f[4]), f[6], attrs(f[8])))
    return rows


def check_assembly(rows, sizes):
    seen = {}
    for chrom, kind, start, end, strand, a in rows:
        if kind == 'chromosome':
            seen[RENAME.get(chrom, chrom)] = end
    if not seen:
        sys.exit('no chromosome features in the GFF; cannot verify the assembly')
    bad = {c: (n, sizes.get(c)) for c, n in seen.items() if sizes.get(c) != n}
    missing = set(sizes) - set(seen)
    if bad or missing:
        sys.exit(f'SGD release does not match {GENOME}: lengths differing (SGD, hub) '
                 f'{bad}; absent from the GFF {sorted(missing)}. The sequence has '
                 f'changed, so these features cannot be placed on this assembly.')
    return len(seen)


def bed12(chrom, strand, blocks, name, colour, thick, extra):
    """One BED12+ line from a set of 1-based inclusive blocks."""
    blocks = sorted(blocks)
    start = blocks[0][0] - 1
    end = blocks[-1][1]
    sizes = [b[1] - b[0] + 1 for b in blocks]
    starts = [(b[0] - 1) - start for b in blocks]
    t0, t1 = (thick if thick else (start, start))
    return '\t'.join([chrom, str(start), str(end), name, '0', strand,
                      str(t0), str(t1), colour, str(len(blocks)),
                      ','.join(map(str, sizes)), ','.join(map(str, starts))] + extra)


def build(rows):
    """Return (protein-coding BED lines, tRNA BED lines, counts by class)."""
    genes, trna_genes = {}, {}
    transcript_of = {}          # transcript id -> gene id
    cds, exons = {}, {}

    for chrom, kind, start, end, strand, a in rows:
        if kind == 'gene':
            genes[a['ID']] = (chrom, strand, a)
        elif kind in RNA_CLASSES:
            trna_genes[a['ID']] = (chrom, strand, a, kind, (start, end))
        elif kind == 'mRNA' or kind in RNA_TRANSCRIPTS:
            for parent in a.get('Parent', '').split(','):
                transcript_of[a['ID']] = parent
        elif kind == 'CDS':
            for parent in a.get('Parent', '').split(','):
                cds.setdefault(parent, set()).add((start, end))
        elif kind == 'noncoding_exon':
            for parent in a.get('Parent', '').split(','):
                exons.setdefault(parent, set()).add((start, end))

    def blocks_for(gene_id, table):
        """Union of a gene's blocks across its transcripts."""
        found = set()
        for tid, gid in transcript_of.items():
            if gid == gene_id:
                found |= table.get(tid, set())
        return found

    coding, classes = [], {}
    for gid, (chrom, strand, a) in genes.items():
        blocks = blocks_for(gid, cds)
        if not blocks:
            continue
        cls = a.get('orf_classification', '')
        classes[cls] = classes.get(cls, 0) + 1
        span = (min(b[0] for b in blocks) - 1, max(b[1] for b in blocks))
        coding.append(bed12(RENAME.get(chrom, chrom), strand, blocks, gid,
                            ORF_COLOUR.get(cls, '0,0,0'), span,
                            [a.get('gene', ''), cls,
                             a.get('dbxref', '').replace('SGD:', ''),
                             a.get('display', '')]))

    trna, rna_classes = [], {}
    for gid, (chrom, strand, a, kind, span) in trna_genes.items():
        label, colour = RNA_CLASSES[kind]
        # Most RNA genes hang their exons off a transcript child; where SGD gives none,
        # the gene feature's own span is the only model there is, and falling back to
        # it beats silently dropping a gene the hub says it serves.
        blocks = blocks_for(gid, exons) or {span}
        rna_classes[label] = rna_classes.get(label, 0) + 1
        trna.append(bed12(RENAME.get(chrom, chrom), strand, blocks, gid,
                          colour, None,
                          [a.get('gene', ''), label,
                           a.get('dbxref', '').replace('SGD:', ''),
                           a.get('display', '')]))
    return coding, trna, classes, rna_classes


def write_bigbed(lines, as_text, n_extra, out_path, chrom_sizes, work):
    bed = os.path.join(work, os.path.basename(out_path) + '.bed')
    asf = os.path.join(work, os.path.basename(out_path) + '.as')
    with open(bed, 'w') as f:
        # bedToBigBed needs the file sorted by sequence then start, and will not sort
        # it for us; the sort is the usual reason a conversion fails late.
        for line in sorted(lines, key=lambda l: (l.split('\t')[0], int(l.split('\t')[1]))):
            f.write(line + '\n')
    with open(asf, 'w') as f:
        f.write(as_text)
    cmd = ['bedToBigBed', f'-type=bed12+{n_extra}', '-tab', f'-as={asf}',
           bed, chrom_sizes, out_path]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f'bedToBigBed failed: {r.stderr.strip()}\n  {" ".join(cmd)}')
    return bed


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--hub-root', required=True)
    ap.add_argument('--raw-dir', default=None)
    args = ap.parse_args()

    chrom_sizes = os.path.join(args.hub_root, f'{GENOME}.chrom.sizes')
    if not os.path.exists(chrom_sizes):
        sys.exit(f'{chrom_sizes} not found: add the {GENOME} assembly to the hub first')
    sizes = {}
    for line in open(chrom_sizes):
        if line.strip():
            n, s = line.split()
            sizes[n] = int(s)

    raw_dir = args.raw_dir or os.path.join(args.hub_root, os.pardir, 'raw', DATASET)
    os.makedirs(raw_dir, exist_ok=True)
    out_dir = os.path.join(args.hub_root, GENOME, DATASET)
    os.makedirs(out_dir, exist_ok=True)

    gff = fetch(GFF_URL, os.path.join(raw_dir, 'saccharomyces_cerevisiae.gff.gz'))
    release, produced = '', ''
    with gzip.open(gff, 'rt') as fh:
        for line in fh:
            if not line.startswith('#'):
                break
            if line.startswith('#!assembly'):
                release = line.split(maxsplit=1)[1].strip()
            if line.startswith('#!date-produced'):
                produced = line.split(maxsplit=1)[1].strip()

    rows = read_gff(gff)
    n_seq = check_assembly(rows, sizes)
    coding, trna, classes, rna_classes = build(rows)
    print(f'SGD {release} ({produced}): {n_seq} sequences verified against {GENOME}')
    print(f'  protein-coding {len(coding)} ' +
          ', '.join(f'{k} {v}' for k, v in sorted(classes.items())))
    print(f'  RNA genes {len(trna)} ' +
          ', '.join(f'{k} {v}' for k, v in sorted(rna_classes.items())))

    write_bigbed(coding, GENE_AS, 4, os.path.join(out_dir, 'sgd.protein_coding.bb'),
                 chrom_sizes, raw_dir)
    write_bigbed(trna, TRNA_AS, 4, os.path.join(out_dir, 'sgd.rna_genes.bb'),
                 chrom_sizes, raw_dir)
    for name, text in (('sgdGene.as', GENE_AS), ('sgdRnaGene.as', TRNA_AS)):
        with open(os.path.join(out_dir, name), 'w') as f:
            f.write(text)
    # The plain BED is shipped beside the bigBed: it is what a reader downloads.
    for src, name in ((coding, 'sgd.protein_coding.bed'), (trna, 'sgd.rna_genes.bed')):
        with open(os.path.join(out_dir, name), 'w') as f:
            for line in sorted(src, key=lambda l: (l.split('\t')[0], int(l.split('\t')[1]))):
                f.write(line + '\n')

    with open(gff, 'rb') as f:
        digest = hashlib.sha256(f.read()).hexdigest()
    dataset = {
        'dataset_id': f'{DATASET}-saccer3',
        'title': f'Gene annotation from the Saccharomyces Genome Database ({release})',
        'organism': 'Saccharomyces cerevisiae',
        'assembly': GENOME,
        'bin_size_bp': None,
        'measurement': 'reference gene annotation, not a measurement',
        'scale': 'features, not a signal',
        'source': {
            'archive_series': [],
            'url': GFF_URL,
            'release': release,
            'date_produced': produced,
            'sha256': digest,
            'citation': ('Saccharomyces Genome Database (SGD), '
                         'https://www.yeastgenome.org/'),
        },
        'provenance': {
            'pipeline': f'pipelines/genes-{DATASET}',
            'generated': datetime.datetime.now(datetime.timezone.utc)
                                 .strftime('%Y-%m-%dT%H:%M:%SZ'),
            'validation': (
                f'the GFF is {release} while {GENOME} is built on R64-1-1; those differ '
                f'in annotation and not in sequence, and the pipeline verifies it rather '
                f'than assuming it - all {n_seq} chromosome features match '
                f'{GENOME}.chrom.sizes exactly, and the run aborts otherwise. Only '
                f'chrmt is renamed, to chrM.'),
            'caveat': ('UCSC serves its own SGD Genes track on sacCer3, built from a '
                       'different release by a different pipeline; these tracks exist '
                       'for the embedded JBrowse panel, which has no native annotation, '
                       'and are not added to the hub overview for that reason'),
        },
        'calls': None,
        'counts': {'protein_coding': len(coding), 'rna_genes': len(trna),
                   'orf_classification': classes, 'rna_class': rna_classes},
    }
    with open(os.path.join(out_dir, 'dataset.json'), 'w') as f:
        json.dump(dataset, f, indent=2)
        f.write('\n')
    print(f'wrote {out_dir}')


if __name__ == '__main__':
    main()
