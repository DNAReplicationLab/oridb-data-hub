# Adding a dataset

A dataset arrives as a set of bigWig files on an assembly the hub already serves, plus a
`dataset.json` describing where it came from.

## 1. Name it

Use a short, stable identifier: first author and year is usually right (`muller2014`), or
the source accession if there is no paper. It becomes the directory name and the prefix of
every container track, so it should not change later.

## 2. Generate the stanza

```bash
python3 bin/make_trackhub.py \
    --hub-root . \
    --genome sacCer3 \
    --dataset muller2014 \
    --bigwig-dir /path/to/pipeline/results/bigwig \
    --descriptions descriptions.tsv \
    --dataset-json /path/to/dataset.json
```

`descriptions.tsv` is two columns — track name, then the long label a reader sees. Without
it the track name is used, which is serviceable but terse.

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

Grouping is by experimental variant, inferred from the track filename: `sortseq`,
`syncseq`, `mfaseq`, `hu`, `trep`. A variant with more than one `syncseq` track becomes an
overlaid `multiWig` container; everything else becomes a composite with one row per track.

If a new dataset uses a variant the generator does not know, add it to `VARIANT` in
`bin/make_trackhub.py` with a label and a one-line description of what was measured. Prefer
that to hand-editing a generated `trackDb.txt` — the file carries a header saying it is
generated, and the next regeneration will discard manual edits.

## 4. Choose view limits

Pass `--viewlimits` as `variant=lo:hi` pairs. The default is `1:2` for copy-number tracks,
which is the relative copy number scale, and `25:85` for Trep in minutes.

Set these from what the measurement means, not from the data's extremes. A track scaled to
its own min and max looks dramatic and compares badly with the track above it; a track on a
scale the reader can interpret — unreplicated to fully replicated — compares directly with
every other dataset in the hub.

## 5. Validate and commit

```bash
python3 bin/validate_hub.py --hub-root . --chrom-sizes sacCer3.chrom.sizes
git add sacCer3/<dataset> sacCer3/trackDb.txt
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
