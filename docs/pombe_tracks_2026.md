# The first *S. pombe* tracks in the OriDB hub: Pu-seq and Mcm2-ChEC

Added 23 November 2026. Two deposited datasets, converted and served on ASM294v2; no
re-analysis, no remapping of reads. Every number below came out of the conversion itself
and is also recorded in the two `dataset.json` files.

Until now the hub served nothing on ASM294v2 except OriDB's own 741 pombe origin records.
These two datasets give the pombe side its first measurements, and between them the two
axes a tier rule needs: **firing** (Pu-seq) and **licensing** (Mcm2-ChEC).

## Why these two

| axis | dataset | what it measures | why it, and not another |
|---|---|---|---|
| firing, timing, fork direction | Daigaku et al. 2015, GSE62108 | polymerase usage | the reference pombe polymerase-usage dataset; later Pu-seq papers normalise against it, and its four derived profiles are published ready to serve |
| licensing | Foss et al. 2021, GSE150800 | MCM2-7 occupancy by ChEC | the only pombe MCM occupancy dataset found on ASM294v2 with a free-MNase control deposited beside it |

Two assays that fail differently: Pu-seq sees an origin only if it fires often enough to
bend the polymerase-usage profile, ChEC sees MCM loading whether or not the origin fires.
That is the point of requiring corroboration between them.

## 1. Pu-seq (Daigaku 2015)

**Source.** Daigaku Y *et al.* A global profile of replication polymerase usage. *Nat
Struct Mol Biol* 2015; PMID 25664722; doi 10.1038/nsmb.2962. GEO GSE62108, the four
series-level derived files. The 25 per-sample strand-resolved wigs in `GSE62108_RAW.tar`
are **not** served: they are the input the four profiles were derived from.

**Assembly.** The submitter states ASM294v2 and the coordinates agree with it: the highest
interval on each chromosome falls within one bin of the chromosome end.

| track | file | bins | bin size | value range |
|---|---|---|---|---|
| Trep | `puseq.trep.bw` | 12,556 | 1 kb | 24.16 – 56.39 minutes |
| Leftward forks | `puseq.leftward_fork.bw` | 12,568 | 1 kb | 0.063 – 0.939 |
| Termination | `puseq.termination.bw` | 41,860 | 300 bp | 0 – 13.37 % |
| Origin calls | `puseq.origins.bb` | 1,145 calls | 300 bp | 0.067 – 83.81 % efficiency |

Three things were done to the files, and nothing else:

1. **Sequences renamed.** `chr1`, `chr2`, `chr3` → `I`, `II`, `III`, which is what PomBase,
   this hub and OriDB's pombe records use.
2. **Edges clipped.** The last bin of a chromosome is written at the full 1 kb span even
   where the sequence ends inside it. Four intervals were clipped to the chromosome end
   (one in Trep, three in the leftward-fork profile) and **one was dropped**: a bin at
   `I:5,579,501–5,580,500` starts 368 bp past the end of chromosome I, so there is nothing
   to clip it to.
3. **The termination profile was re-binned.** It is published as a 1 kb window stepped
   every 300 bp, so consecutive intervals overlap by 700 bp — which bigWig forbids. Each
   value is now written to a 300 bp bin centred on its own window. No value changed and
   the sampling grid is preserved; the window it was computed over (1 kb) is recorded in
   `dataset.json`.

**Origin efficiency is served as a call set, not a profile.** The published bedGraph is
1,145 discrete 300 bp origin calls rather than a continuous signal, and six positions are
listed twice with different efficiencies. A bigWig cannot hold overlapping intervals, and
merging the duplicates would discard a published value, so the file became a bigBed in
which every record survives. The score column scales efficiency to 0–1000 across the
published range and the raw percentage is kept as its own field.

**This is the call set the pombe tier rule will register.** Those six duplicated positions
have to be resolved by the registry step, not here.

## 2. Mcm2-ChEC (Foss 2021)

**Source.** Foss EJ *et al.* Chromosomal Mcm2-7 distribution and the genome replication
program in species from yeast to humans. *PLoS Genet* 2021; PMID 34473702;
doi 10.1371/journal.pgen.1009714. GEO GSE150800 — a five-species series, of which the
pombe arm is three files of **per-base read depths**, 12.57 million rows each. There is no
bedGraph for pombe in the archive, only for the mouse and human arms.

Eight tracks, 30.4 MB in total:

| track | source column | mean depth per base |
|---|---|---|
| `mcm2chec.hu_rep1.all.bw` | HU rep1, all fragment sizes | 46.8 |
| `mcm2chec.g1_rep1.all.bw` | G1 rep1, all fragment sizes | 54.1 |
| `mcm2chec.hu_rep1.51_100bp.bw` | HU rep1, 51–100 bp | 11.5 |
| `mcm2chec.g1_rep1.51_100bp.bw` | G1 rep1, 51–100 bp | 9.4 |
| `mcm2chec.hu_rep1.151_200bp.bw` | HU rep1, 151–200 bp | 13.1 |
| `mcm2chec.g1_rep1.151_200bp.bw` | G1 rep1, 151–200 bp | 12.6 |
| `mcm2chec.hu_rep2.all.bw` | HU rep2 | 146.5 |
| `free_mnase.rep1.bw` | free MNase rep1 | 57.8 |

Two decisions, both reversible and both recorded:

1. **Binned to 25 bp.** At one base per interval a single track is 38.2 MB of bigWig and
   the eight would be 300 MB; at 25 bp each is about 3.7 MB. 25 bp is still an order of
   magnitude finer than the 300 bp – 1 kb Pu-seq bins beside it. ChEC resolves finer than
   this, and `convert.py --bin 1` reproduces the single-base version for anyone who needs
   it — it just should not be what GitHub Pages serves to every reader.
2. **Scaled to each library's own mean depth**, so 1 is average coverage. The libraries
   differ up to 15-fold in depth (9.4 to 146.5 reads per base), and the containers that
   overlay HU against G1, or rep1 against rep2, are unreadable otherwise. Every raw mean
   depth is in `dataset.json`, so the published counts can be recovered by multiplying
   back.

**The free-MNase control is served beside the samples**, in its own container. ChEC signal
is accessibility plus tethered cleavage, and a hub that served only the samples would
invite the control to be skipped.

**Display note.** `viewLimits` are capped at 10 on purpose. ChEC peaks reach 300 times mean
coverage, so a scale that fits the tallest peak flattens everything else;
`bin/validate_hub.py` warns about the clipping, and for these ten tracks the warning is
expected rather than a defect to fix.

## 3. What the hub now draws by default on ASM294v2

The overview section gained four pombe entries — Trep, fork direction, Pu-seq origin calls
and Mcm2 occupancy in HU — which join OriDB's own origin records. Five tracks on first
load, against a budget of ten.

## 4. What this does not yet give

The licensing axis has a **profile but no call set**: MCM sites have to be called from the
ChEC signal before an origin can be said to be licensed, and that calling rule is a
scientific choice with one HU and one G1 replicate behind it. Until it exists, a pombe site
can be corroborated by Pu-seq alone, which is one axis and therefore not enough for
*likely* under the tier rule.

## 5. Reproducing

```bash
python3 pipelines/puseq-daigaku2015/convert.py  --hub-root . --raw-dir /tmp/puseq
python3 pipelines/mcm2chec-foss2021/convert.py  --hub-root . --raw-dir /tmp/foss
for genome in $(awk '/^genome /{print $2}' genomes.txt); do
  for dir in "$genome"/*/; do
    dataset=$(basename "$dir")
    spec="docs/containers_${dataset}_${genome}.json"
    [ -f "$spec" ] || spec="docs/containers_${dataset}.json"
    python3 bin/make_trackhub.py --hub-root . --genome "$genome" --dataset "$dataset" --spec "$spec"
  done
done
python3 bin/validate_hub.py --hub-root .
```

Both pipelines download what they need from GEO and cache it; both need `bedGraphToBigWig`
(and the Pu-seq one `bedToBigBed`) on `PATH`. The conversion is deterministic, which is why
CI can re-run the generator and fail on a diff.
