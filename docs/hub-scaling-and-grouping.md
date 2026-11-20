# Keeping the hub usable as it grows

How the OriDB data hub should be structured so that adding datasets does not make it
unusable in the UCSC Genome Browser. Written against UCSC's own hub documentation
([trackDb reference](https://genome.ucsc.edu/goldenPath/help/trackDb/trackDbHub.html),
[grouping quick start](https://genome.ucsc.edu/goldenPath/help/hubQuickStartGroups.html),
[track hub user guide](https://genome.ucsc.edu/goldenPath/help/hgTrackHubHelp.html)).

## What actually overloads, and what does not

Hub data is fetched by byte range on demand — UCSC retrieves only the part of a bigWig or
bigBed covering the window in view, and caches it. Total repository size is therefore not the
constraint, and neither is the number of datasets *stored*. Three things do degrade:

1. **Tracks drawn at once.** Every track visible in the default view is a fetch, a decode and a
   rendered row. The hub currently serves 11 top-level tracks from 2 datasets, all declared
   `visibility full` or `pack`. At that rate the 78 priority-1 and priority-2 datasets in the
   candidate catalogue would come to several hundred tracks, all on by default. That is the
   overload, and it arrives long before any storage limit.
2. **Track-controls navigability.** A flat list of several hundred stanzas under one hub heading
   is unusable even when nothing is drawn.
3. **trackDb assembly.** UCSC checks the timestamp of every included file. Its own guidance is
   that include statements "may cause slower load times on high-latency connections", with
   concatenation suggested instead. The hub uses one `include` per dataset, so this cost grows
   linearly with the catalogue.

## The grouping mechanisms UCSC provides

UCSC distinguishes three containers by how closely related the tracks are:

| Container | Use it when | Already used here |
|---|---|---|
| `container multiWig` | all children are signal tracks that make sense overlaid on one scale | yes — sync-seq timepoints, Ok-seq strands |
| `compositeTrack on` | children are configured together or have internal subgrouping (strain, assay, timepoint) | no |
| `superTrack on` | children share no configuration and you simply want them shown together — "a set of composite tracks or a mix of bigBed, bigWig and bam" | no |

A supertrack is the answer to "slightly separated sections": it renders as its own labelled,
collapsible folder in the track controls, with its own position in the list. Supertracks can
contain composites and multiWigs, but not the reverse. The declaring stanza behaves like a
light switch that is off unless you add `show`: `superTrack on show` displays by default,
`superTrack on` alone does not. Children declare `parent <superTrackName>` and set their own
`visibility`.

For collections that outgrow even that, UCSC documents the **faceted composite**
(`compositeTrack faceted` plus `metaDataUrl`), intended for datasets with too many tracks to
coerce into a subgroup structure and which mix signals and peaks under one heading. That is a
fair description of where this hub ends up, and the catalogue TSV is the natural source for the
metadata table. It is the right second step, not the first.

## Proposed layout for this hub

Two levels, mapping onto the structure the dataset catalogue already uses:

- **Level 1 — one supertrack per data type.** Replication timing; Initiation mapping; Licensing
  and replisome footprints; Fork direction and polymerase usage; Fork pausing and termination;
  Re-replication; ssDNA and replication stress. These are exactly the values of the `data_type`
  column in `oridb_hub_candidate_datasets.tsv`, so the browser's sections and the catalogue stay
  in step by construction, and a new dataset's section follows from its catalogue row.
- **Level 2 — the existing per-dataset containers become children.** Each dataset's current
  multiWig or bigBed stanzas gain one `parent <superTrack>` line. Nothing else about the
  per-dataset directory changes.

```
# sacCer3/superTracks.txt  (new, generated)
track oridb_timing
superTrack on show
shortLabel Replication timing
longLabel Replication timing and copy-number profiles through S phase
priority 1

track oridb_footprint
superTrack on
shortLabel Licensing footprints
longLabel ORC and MCM2-7 occupancy, origin licensing
priority 3
...

# sacCer3/muller2014/trackDb.txt  (one added line per top-level stanza)
track muller2014_timecourse
parent oridb_timing
container multiWig
...
```

This preserves the rule the repo was built on — adding a dataset never edits another dataset's
config — because the only shared file touched is the generated supertrack list.

## Default visibility is the real control

Grouping makes the hub navigable; visibility policy is what stops it being slow. Proposed rules:

1. **At most one supertrack carries `show`.** Everything else is `superTrack on`, so a first-time
   visitor gets one section drawn and the rest one click away.
2. **A target for the default view: no more than 8–10 drawn tracks, whatever the catalogue size.**
   This is the number to hold constant as datasets accumulate; it should be asserted in
   `bin/validate_hub.py` so a dataset cannot be added with everything switched on.
3. **Within a dataset, one or two tracks default to `full`**, the rest `hide` or `dense`. The
   current hub has 11 of 11 visible, which is the habit to break now rather than at 50 datasets.
4. **`onlyVisibility <mode>`** caps tracks that only make sense in one display mode, without
   discarding the user's per-subtrack settings.
5. **`priority`** orders sections; keep reference datasets above perturbation panels.

A dataset such as GSE249817 (258 ORC ChEC-seq samples across 40-plus mutants) should enter as a
single composite with subgroups for genotype and replicate and `dimensions dimX=genotype`,
with only the wild-type subtrack on — not as 258 stanzas.

## Other limits worth knowing

- `bigBedMaxItems` is a server-side cap on features loaded from a bigBed; `maxItems` and
  `maxWindowToDraw` are the per-track equivalents for item tracks.
- UCSC re-checks hub file timestamps every 300 s, so a push can take five minutes to appear;
  `udcTimeout=5` in the browser URL cuts that to five seconds while testing.
- `hubCheck` via the Hub Development tab validates settings but, with `-noTracks`, does not fetch
  the data files — `bin/validate_hub.py` remains the check that the bigDataUrls resolve.

## Suggested order of work

1. Add supertrack generation to `bin/make_trackhub.py` and emit `sacCer3/superTracks.txt`;
   add the `parent` line to each dataset's top-level stanzas. Small, mechanical, and it can be
   done before any new dataset arrives.
2. Set the default-visibility policy above, and add the drawn-track budget as a check in
   `bin/validate_hub.py`.
3. Decide whether `make_trackhub.py` should emit one concatenated `sacCer3/trackDb.txt` instead
   of per-dataset includes. The per-dataset files stay the authoring unit either way; this only
   changes what is served, and it removes the per-include timestamp checks UCSC warns about.
4. Revisit faceted composites once the hub carries more than roughly twenty datasets, using the
   catalogue TSV as the metadata table.
