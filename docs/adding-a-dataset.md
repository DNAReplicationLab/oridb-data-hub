# Adding a dataset

A dataset arrives as a set of bigWig files on an assembly the hub already serves, plus a
`dataset.json` describing where it came from.

## 1. Name it

Use a short, stable identifier: first author and year is usually right (`muller2014`), or
the source accession if there is no paper. It becomes the directory name and the prefix of
every container track, so it should not change later.

## 2. Generate the stanza

Containers are **comparisons**, not file groupings. Which profiles belong side by side is
a judgement — a haploid against a diploid, one platform against another, a mutant against
its wild type — and the same track legitimately appears in more than one comparison. That
cannot be inferred from filenames, so it is declared in a spec file per dataset and the
generator only renders it.

Start from a generated skeleton, one container per track:

```bash
python3 bin/make_trackhub.py --init \
    --hub-root . --genome sacCer3 --dataset mystudy2026 \
    --bigwig-dir /path/to/pipeline/results/bigwig \
    --spec docs/containers_mystudy2026.json
```

Then edit that file: merge tracks into comparison containers, give each a `shortLabel`
and `longLabel`, set `graphType` (`bar` for a time course, `points` otherwise), set
`viewLimits`, and choose colours. Render it:

```bash
python3 bin/make_trackhub.py \
    --hub-root . --genome sacCer3 --dataset mystudy2026 \
    --spec docs/containers_mystudy2026.json \
    --bigwig-dir /path/to/pipeline/results/bigwig \
    --dataset-json /path/to/dataset.json
```

Every container is rendered as a `multiWig` `solidOverlay`, because overlaying the
profiles being compared is what makes a difference visible — stacked rows leave the eye
to do a subtraction it cannot do.

## Call sets

A set of origin calls is not a profile and does not belong in a comparison container —
there is nothing to overlay it with. Declare it in the spec's top-level `tracks` array
instead, and the generator renders it as a standalone `bigBed` stanza:

```json
"tracks": [
  {"id": "trep_peaks", "file": "syncseq.trep.peaks.bb", "type": "bigBed 8 +",
   "shortLabel": "Trep peaks", "longLabel": "...", "spectrum": true,
   "scoreMin": 100, "visibility": "pack", "priority": 7}
]
```

`spectrum` shades each feature by its score column, so a weak call looks weak. Ship the
plain BED beside the `.bb` — it is what a reader downloads — plus the `.as` that declares
any fields beyond the standard ones, and record how the calls were made in a doc of their
own. [docs/trep_peak_calls.md](trep_peak_calls.md) is the worked example.

### Settings that matter more than they look

`maxHeightPixels` is `max:default:min` and the **middle** value is the height actually
drawn. A default of 32 renders a signal track at a third of its usable size, which is the
usual reason a track looks cramped. The hub uses `100:100:8`.

`graphTypeDefault` is `bar` or `points`, nothing else. Filled bars suit a dense time
course where the envelope is the message; points suit a comparison of two or three
profiles, where fills would hide each other.

`negateValues on` flips the sign of the data. The Trep track uses it so that *early*
replication plots high, consistent with the copy-number tracks where high means more
replicated — otherwise an origin that is a peak everywhere else becomes a trough in Trep.
Note that UCSC draws negative values in `altColor` rather than `color`, so a negated track
needs both set.

The generator copies the bigWigs in, writes `<genome>/<dataset>/trackDb.txt`, and rewrites
`<genome>/trackDb.txt` as one `include` per dataset directory present. Re-running it for an
existing dataset refreshes that dataset only.

### One path rule to know

UCSC resolves the two kinds of relative path in a hub differently, and the difference is
easy to get wrong:

- an `include` path is relative to **the file containing it**, so `sacCer3/trackDb.txt`
  carrying `include muller2014/trackDb.txt` reads `sacCer3/muller2014/trackDb.txt`, as you
  would expect;
- a `bigDataUrl` is relative to **the trackDb declared in `genomes.txt`** — that is
  `sacCer3/trackDb.txt` — regardless of which included file the stanza came from.

So a track living in `sacCer3/muller2014/` needs `bigDataUrl muller2014/<file>.bw`, not
`bigDataUrl <file>.bw`. Getting this wrong produces a hub whose files are all present and
all individually reachable, and which still fails in the browser with
`couldn't open .../sacCer3/<file>.bw`. `bin/make_trackhub.py` applies the prefix and
`bin/validate_hub.py` resolves paths the same way UCSC does, so the mistake is caught
before the hub is served rather than after.

## 3. How tracks are grouped

Two levels, and a dataset only has to think about the second.

**Sections** are UCSC supertracks — the folders a reader sees in the track controls. They
are declared once for the whole hub in `docs/sections.json` and follow the `data_type`
values in `docs/oridb_hub_candidate_datasets.tsv`, so the browser's layout and the dataset
census stay in step: `timing`, `initiation`, `footprint`, `forkdir`, `pausing`, `rerep`,
`ssdna`. A section with no dataset in it is not rendered.

**A dataset says only which section each container belongs to.** Put a default in the spec
and override per container where a dataset spans more than one kind of measurement:

```json
{
  "dataset": "mystudy2026",
  "section": "timing",
  "containers": [
    { "id": "timecourse", "...": "..." },
    { "id": "hu", "section": "initiation", "...": "..." }
  ]
}
```

The generator writes `parent oridb_<section>` into each top-level stanza and regenerates
`sacCer3/sections.trackDb.txt`. Naming a section that `docs/sections.json` does not define
is an error, not a silently orphaned track. Adding a dataset still never edits another
dataset's configuration.

**Every section is closed by default** (`superTrack on`). The only section open on first
load is the overview, which holds one canonical track per measurement type, also declared
in `docs/sections.json`. That is deliberate: it keeps the number of tracks the browser
draws independent of how many datasets the hub carries. If a new dataset is the best
available reference for its measurement type, propose swapping it into the overview rather
than opening its section.

`bin/validate_hub.py` fails if more than `--max-default-tracks` (10) rows would be drawn on
first load, and prints the per-section counts so you can see what a reader meets when they
open a folder. A dataset that arrives with everything switched on will not pass CI.

## 4. Choose view limits

Set `viewLimits` per container in the spec. The conventions in use are `0.9:2.1` for
copy-number profiles, `0.9:1.5` for MFA-seq, `-85:-25` with `negateValues` for Trep in
minutes, and `-1:1` for fork directionality.

Set these from what the measurement means, not from the data's extremes. A track scaled to
its own min and max looks dramatic and compares badly with the track above it; a track on a
scale the reader can interpret — unreplicated to fully replicated — compares directly with
every other dataset in the hub.

## 5. Validate and commit

```bash
python3 bin/validate_hub.py --hub-root . --chrom-sizes sacCer3.chrom.sizes
git add sacCer3/<dataset> sacCer3/trackDb.txt sacCer3/sections.trackDb.txt \
        docs/containers_<dataset>.json
git commit -m "Add <dataset>"
```

## Adding a new assembly

Add a `genome`/`trackDb` pair to `genomes.txt` and create `<genome>/` with its dataset
directories. For *S. pombe*, note that UCSC has no native assembly: GenArk hosts
GCA_000002945.3 (labelled ASM294v3), whose chromosome lengths match ASM294v2, so data need
no lift — only sequence-name translation from `chr1`/`chr2`/`chr3` to
`CU329670.1`/`CU329671.1`/`CU329672.1`.

## Size

Git is a poor home for large binaries. The current tracks are ~90 kB each because they are
1 kb bins over a 12 Mb genome, which is fine. If a future dataset brings tracks of hundreds
of megabytes, serve those from object storage and keep only the `trackDb.txt` and
`dataset.json` here, pointing `bigDataUrl` at the absolute URL.
