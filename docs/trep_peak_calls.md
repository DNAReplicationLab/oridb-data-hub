# The Trep peak calls

`sacCer3/muller2014/syncseq.trep.peaks.bed` and its indexed form `.bb` are **derived**
data, not a profile deposited by the pipeline: they are origin calls made from
`syncseq.trep.bw` for OriDB. Everything needed to reproduce them is here.

## What a call is

In the Trep track an origin is a local **minimum** — the earliest-replicating point of a
valley — not a maximum. That was checked against the OriDB Confirmed origins rather than
assumed: the 350 confirmed cerevisiae origins that fall in a fitted bin sit at a median
Trep of 38.2 min against a genome-wide median of 43.3 min, the 25th percentile of the
genome-wide distribution (Mann-Whitney p = 1.1e-36). The track therefore stores minutes
with early low, which is why the `trep` container sets `negateValues on` for display.

## How they were called

1. The profile is read in its native 1 kb bins. Unfitted bins (2,096 of 12,079) are
   linearly interpolated across gaps of up to 3 bins; a longer gap splits the chromosome,
   so no peak is invented across missing data. Calls with an interpolated bin within 1 kb
   of the summit carry `interpolated = 1` in the last column (46 of 212).
2. The profile is smoothed with a Gaussian of sigma = 2 bins (2 kb) and peaks are called
   on its negation with `scipy.signal.find_peaks`, prominence >= 2 minutes, minimum
   separation 5 kb.
3. The feature extent is the band around the summit in which smoothed Trep stays within
   1 minute of the trough minimum — roughly the distance a fork covers in that minute.
   `thickStart`/`thickEnd` mark the summit bin.

Prominence is the parameter that matters. The bin-to-bin Trep step has a robust sigma of
1.36 minutes, which smoothing at 2 kb reduces to about 0.51 minutes, so a 2 minute
prominence is a 4 sigma trough. Minimum separation is nearly inert at this smoothing:
distances from 1 kb to 8 kb all give 212 calls.

## Rebuilding the indexed file

```bash
bedToBigBed -as=sacCer3/muller2014/trepPeaks.as -type=bed8+3 \
    sacCer3/muller2014/syncseq.trep.peaks.bed sacCer3.chrom.sizes \
    sacCer3/muller2014/syncseq.trep.peaks.bb
```

`trepPeaks.as` is reproduced here in full, so the repository is self-contained if that
file is ever lost:

```
table trepPeaks
"Replication origin calls from the Muller et al. (2014) sync-seq Trep profile (local minima of replication time)"
    (
    string  chrom;        "Reference sequence"
    uint    chromStart;   "Start of the call, 0-based"
    uint    chromEnd;     "End of the call: the band within 1 minute of the trough minimum"
    string  name;         "Call identifier, numbered in genome order"
    uint    score;        "Prominence scaled to 0-1000 for shading"
    char[1] strand;       "Not applicable"
    uint    thickStart;   "Summit: midpoint of the 1 kb bin holding the earliest Trep"
    uint    thickEnd;     "Summit end (thickStart + 1)"
    float   trep;         "Smoothed Trep at the summit, minutes from alpha-factor release"
    float   prominence;   "Depth of the trough below its enclosing saddle, minutes"
    uint    interpolated; "1 if an unfitted bin was interpolated within 1 kb of the summit"
    )
```

The BED is written in the sort order `bedToBigBed` requires (chromosome name
lexicographic, then start), while the `name` column is numbered in genome order.
`trepPeaks.as` declares the three fields beyond BED8 — `trep`, `prominence` and
`interpolated` — and is a build input, not a served file.
