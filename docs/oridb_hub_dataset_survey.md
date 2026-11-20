# Candidate datasets for the OriDB v3 data hub — *Saccharomyces cerevisiae*

A census of published, publicly deposited genome-wide replication datasets for *S. cerevisiae*, assembled from NCBI GEO and from ENA/SRA deposits that have no GEO record. **164 datasets are in scope**, split into two groups that are served differently:

- **155 track-hub candidates** — continuous signal or interval data that a genome browser displays directly. These are the subject of most of this document.
- **9 single-molecule datasets** — per-read data that UCSC and comparable browsers do not render usefully. These are catalogued separately in their own section and are destined for the group's own single-molecule visualisation tool, not the track hub. Two of the nine also yield a conventional aggregate profile, so they appear in both routes.

### Already in the hub

Two datasets are already served by the hub, giving 11 tracks in total. They are marked **[in hub]** in the listings below and carry the dataset directory name in the `in_hub` column of the TSV.

| Hub dataset | Accessions | Tracks | Contribution |
|---|---|---|---|
| `muller2014` | GSE42243, GSE48212 | 7 | Sort-seq, sync-seq and MFA-seq timing, reprocessed on sacCer3 |
| `mcguffee2013` | GSE40696 | 4 | Ok-seq fork directionality, initiation signal and initiation-site calls, re-analysed on sacCer3 |

The machine-readable version is `oridb_hub_candidate_datasets.tsv`, one row per dataset. The `category` column holds `hub track` or `single-molecule`, and `serving_route` records which viewer a dataset is destined for.

## Scope

Included: assays that measure the replication process itself — replication timing and copy-number profiling, initiation site and origin-efficiency mapping, licensing and replisome protein occupancy, fork direction and polymerase usage, fork position/pausing/termination, single-molecule replication, ssDNA during S phase, and re-replication. All eras and all strains and conditions are included, with a priority flag separating reference-quality wild-type data from perturbation and legacy-platform datasets.

Excluded: nucleosome and histone-modification mapping, chromatin maturation and histone transfer, cohesin/condensin and Hi-C, meiotic recombination, transcription, and repair-focused break mapping — even where a replication term appears in the title. 12 series that the searches returned and that were examined and rejected on these grounds are listed at the end of the TSV with priority 0 and a reason.

## How the survey was run

1. A repository-first sweep of GEO (`db=gds`, *S. cerevisiae*, series only) across 30 assay-specific query terms returned 522 distinct series; every one was triaged on title, design and assay type rather than on the term that found it.
2. Eight deliberately broader terms were then run as a completeness test. They returned 166 further series, of which five proved genuinely in scope — so the assay-term sweep was close to saturating for this scope.
3. Each candidate's full GEO SOFT record was retrieved, giving per-sample processing text, platform, instrument and the complete supplementary-file inventory.
4. ENA/SRA was searched independently and every hit cross-checked against the GEO set via BioProject. This is what surfaced the nanopore single-molecule studies, which are deposited in ENA and have no GEO record at all — a literature- or GEO-only survey would have missed them.
5. Citations come from the PMID GEO itself records where one exists; where GEO records none, the publication was found by title search and the match is flagged in the `citation_link` column (`GEO-linked`, `title-matched`, `probable match`, `none found`). DOIs were resolved from PubMed, not composed.

## Track-hub candidates

| Data type | P1 | P2 | P3 | Total |
|---|---|---|---|---|
| replication timing | 4 | 17 | 16 | 37 |
| initiation mapping | 1 | 14 | 17 | 32 |
| protein footprinting | 7 | 10 | 21 | 38 |
| fork direction / polymerase usage | 6 | 7 | 9 | 22 |
| fork position, pausing, termination | 0 | 8 | 4 | 12 |
| re-replication | 0 | 3 | 7 | 10 |
| ssDNA / replication stress | 0 | 1 | 3 | 4 |
| **Total** | 18 | 60 | 77 | 155 |

Priority 1 means a dataset usable as a reference track for its assay class in unperturbed cells. Priority 2 means a comparative or perturbation dataset with wild-type controls, or reference data needing reprocessing or a coordinate lift. Priority 3 means peripheral: legacy array platform, a single locus, a non-standard reference genome, or a study whose replication content is incidental.

## The state the data is actually in

Counts below are for the 155 track-hub candidates.

| Current state of the data | Datasets |
|---|---|
| continuous signal track published | 67 |
| tabular processed data published | 38 |
| array probe-level only | 24 |
| interval/call set published | 16 |
| raw reads only | 10 |


| Coordinate system | Datasets |
|---|---|
| sacCer3 / R64 | 86 |
| microarray probe-level | 36 |
| sacCer2 | 7 |
| not stated in the archive record | 6 |
| custom / non-S288C reference | 6 |
| sacCer1 | 5 |
| not stated in GEO | 4 |
| multi-species, build not resolved | 1 |
| SGD coordinates, version not stated | 1 |
| per-species references, version not stated | 1 |
| SGD reference, version not stated | 1 |
| S288C reference, version not stated | 1 |

Three points matter for ingestion planning.

First, **sacCer3 is now the norm but not the rule**: 86 of 155 track-hub candidates state sacCer3/R64 in their own processing description, while 7 state sacCer2 and 5 sacCer1. The sacCer2 group includes the two Okazaki-fragment datasets that define the cerevisiae fork-direction axis, so a lift is unavoidable there.

Second, **36 datasets are microarray probe-level only**. Their coordinates are whatever the array design used, which is typically an SGD build from the submission year rather than a named assembly; serving them means re-deriving positions from the platform record, not lifting a file.

Third, **10 track-hub candidates publish no processed data at all** — the archive holds reads and nothing else.

A caution that applies throughout: the `assembly` column records what the authors *stated*, flagged in `assembly_source`. It is not an independent verification. Where a dataset becomes load-bearing, the build should be confirmed against the file itself before the track is published — this project has already seen one case (GSE40696) where the published signal tracks had to be assigned to sacCer2 empirically by cross-correlation, because the deposit does not say so unambiguously.


## Priority 1 — the reference set

The datasets I would treat as the backbone of the hub: each is the best available genome-wide measurement of its kind in unperturbed cells.


### replication timing

- **GSE36045** — Müller CA 2012, *Genome Res*, [10.1101/gr.139477.112](https://doi.org/10.1101/gr.139477.112). sort-seq across 4 Saccharomyces species + hybrid. 10 samples. *State:* continuous signal track published; per-species references, version not stated. Müller & Nieduszynski 2012; cerevisiae arm directly comparable to Muller 2014.
- **GSE42243** **[in hub]** — Müller CA 2014, *Nucleic Acids Res*, [10.1093/nar/gkt878](https://doi.org/10.1093/nar/gkt878). sort-seq (FACS S/G2 copy number), SOLiD + Illumina. 4 samples. *State:* continuous signal track published; sacCer1 (2003 build). Müller 2014; already in the hub as the muller2014 dataset, reprocessed onto sacCer3 by this project.
- **GSE48212** **[in hub]** — Müller CA 2014, *Nucleic Acids Res*, [10.1093/nar/gkt878](https://doi.org/10.1093/nar/gkt878). sync-seq / MFA-seq / sort-seq time course. 14 samples. *State:* continuous signal track published; sacCer1 (2003 build). Müller 2014 companion series; already in the hub, reprocessed onto sacCer3 by this project.
- **GSE48561** — Hawkins M 2013, *Cell Rep*, [10.1016/j.celrep.2013.10.014](https://doi.org/10.1016/j.celrep.2013.10.014). sort-seq time course (6 S-phase samples, 2 strains). 16 samples. *State:* continuous signal track published; sacCer1 (2003 build). Hawkins 2013 high-resolution profiles + stochastic model.

### initiation mapping

- **GSE71050** — Peace JM 2016, *Genome Res*, [10.1101/gr.196857.115](https://doi.org/10.1101/gr.196857.115). quantitative BrdU-IP-seq (QBU), Fkh1/Fkh2. 64 samples. *State:* tabular processed data published; sacCer3 / R64. Peace 2016; quantitative, absolute origin firing measure.

### protein footprinting

- **GSE150800** — Foss EJ 2021, *PLoS Genet*, [10.1371/journal.pgen.1009714](https://doi.org/10.1371/journal.pgen.1009714). Mcm2-7 MNase (ChEC) across species incl. S. cerevisiae. 22 samples. *State:* continuous signal track published; sacCer3 / R64. Foss 2021; cross-species MCM distribution.
- **GSE16926** — Eaton ML 2010, *Genes Dev*, [10.1101/gad.1913210](https://doi.org/10.1101/gad.1913210). ORC ChIP-seq + nucleosome mapping at origins. 5 samples. *State:* interval/call set published; SGD coordinates, version not stated (2010 submission). Eaton 2010; the ACS/nucleosome reference and an OriDB v2 continuity anchor.
- **GSE173527** — Reuter LM 2024, *Nat Commun*, [10.1038/s41467-024-51538-9](https://doi.org/10.1038/s41467-024-51538-9). ChIP-exo 5.0 of Mcm4 and Orc2, G1 and G2. 8 samples. *State:* continuous signal track published; sacCer3 / R64. Near base-pair footprints.
- **GSE173529** — Weekes C 2026, *Nat Commun*, [10.1038/s41467-026-75695-1](https://doi.org/10.1038/s41467-026-75695-1). SuperSeries: MCM2-7 loading and initial DNA melting at near-bp resolution. 22 samples. *State:* continuous signal track published; sacCer3 / R64. ChIP-exo 5.0 + permanganate; umbrella of GSE173525/173527/173528.
- **GSE240779** — Reuter LM 2024, *Nat Commun*, [10.1038/s41467-024-51538-9](https://doi.org/10.1038/s41467-024-51538-9). ORC and MCM ChIP-seq/ChIP-exo, loading-dependent ORC release. 44 samples. *State:* continuous signal track published; sacCer3 / R64. Reuter 2024; genome-wide ORC and MCM2-7 calls (published table is conservative).
- **GSE242131** — Foss EJ 2024, *Elife*, [10.7554/eLife.88087](https://doi.org/10.7554/eLife.88087). MCM2/4/6 ChEC-seq + MNase-seq; 17,618 MCM binding sites. 13 samples. *State:* tabular processed data published; sacCer3 / R64. Foss 2024; publishes a single-base call set, sacCer3 - highest-resolution origin call set available.
- **GSE249817** — Chappleboim M 2024, *Nucleic Acids Res*, [10.1093/nar/gkae249](https://doi.org/10.1093/nar/gkae249). ORC ChEC-seq, wild-type and >40 ORC mutants. 258 samples. *State:* continuous signal track published; sacCer3 / R64. Chappleboim 2024; 258 samples, R64/sacCer3.

### fork direction / polymerase usage

- **GSE125855** — Zhou ZX 2019, *Nat Commun*, [10.1038/s41467-019-11995-z](https://doi.org/10.1038/s41467-019-11995-z). HydEn-seq polymerase usage, polymerase variants. 58 samples. *State:* continuous signal track published; multi-species, build not resolved. Zhou 2019; extended HydEn-seq panel.
- **GSE134225** — Sriramachandran AM 2020, *Mol Cell*, [10.1016/j.molcel.2020.03.027](https://doi.org/10.1016/j.molcel.2020.03.027). GLOE-seq SuperSeries, SSBs and replication patterns. 96 samples. *State:* continuous signal track published; sacCer3 / R64. Sriramachandran 2020; strand-specific SSB/replication maps.
- **GSE154811** — Kara N 2021, *PLoS Biol*, [10.1371/journal.pbio.3000886](https://doi.org/10.1371/journal.pbio.3000886). TrAEL-seq, WT and mutants (+ hESC controls). 46 samples. *State:* tabular processed data published; sacCer3 / R64. Kara 2021; fork direction and position from 3' ends.
- **GSE33786** — Smith DJ 2012, *Nature*, [10.1038/nature10895](https://doi.org/10.1038/nature10895). Okazaki fragment sequencing, WT and pol32. 4 samples. *State:* interval/call set published; sacCer2 (2008 build). Smith & Whitehouse 2012; the original cerevisiae Okazaki dataset.
- **GSE40696** **[in hub]** — McGuffee SR 2013, *Mol Cell*, [10.1016/j.molcel.2013.03.004](https://doi.org/10.1016/j.molcel.2013.03.004). Ok-seq (Okazaki fragment sequencing), OEM initiation/termination. 5 samples. *State:* continuous signal track published; sacCer2 (2008 build) (custom reference for some samples). McGuffee 2013; the cerevisiae Ok-seq reference. Published .sgr signal is sacCer2.
- **GSE62181** — Clausen AR 2015, *Nat Struct Mol Biol*, [10.1038/nsmb.2957](https://doi.org/10.1038/nsmb.2957). HydEn-seq ribonucleotide mapping, 14 strains. 29 samples. *State:* continuous signal track published; not stated in GEO. Clausen 2015; polymerase usage reference.

## Priority 2 — comparative and perturbation datasets


### replication timing (17)

| Accession | Study | Method | Current state | Notes |
|---|---|---|---|---|
| GSE135178 | Batrakou DG 2020, *Nat Protoc*, [10.1038/s41596-019-0287-7](https://doi.org/10.1038/s41596-019-0287-7) | sort-seq protocol reference (S-phase sorted, 1 kb bins) | interval/call set published; sacCer3 / R64 | Batrakou 2020 methods paper; small demonstration dataset |
| GSE17120 | Koren A 2010, *Genome Res*, [10.1101/gr.102764.109](https://doi.org/10.1101/gr.102764.109) | array CGH (S vs G1), mrc1 scaling of the timing program | tabular processed data published; custom / non-S288C reference | Agilent 4x44K; legacy platform |
| GSE199447 | Santos MM 2022, *Genome Biol*, [10.1186/s13059-022-02788-7](https://doi.org/10.1186/s13059-022-02788-7) | copy-number replication profiles, global early-replication strain | tabular processed data published; sacCer3 / R64 | SSDDCS strain; whole-genome timing perturbation |
| GSE205068 | Kwan EX 2023, *Cell Rep*, [10.1016/j.celrep.2023.112161](https://doi.org/10.1016/j.celrep.2023.112161) | heavy/light density-transfer replication timing incl. rDNA | tabular processed data published; microarray probe-level; coordinates set by platform GPL10930 | Isotopic density transfer; rDNA timing coordinates anaphase |
| GSE212974 | Weiβ M 2023, *Cell Rep*, [10.1016/j.celrep.2023.112045](https://doi.org/10.1016/j.celrep.2023.112045) | copy-number origin activity after targeted Ask1 degradation/recruitment | continuous signal track published; sacCer3 / R64 | Single-copy locus proteomics of early vs late origins |
| GSE32002 | Gispan A 2014, *BMC Biol*, [10.1186/s12915-014-0079-z](https://doi.org/10.1186/s12915-014-0079-z) | array CGH (S vs G1), 18 strains, several stresses | tabular processed data published; microarray probe-level; coordinates set by platform GPL10930/GPL14566 | Large legacy timing matrix (Agilent CGH); coordinates follow array design |
| GSE33677 | Knott SR 2012, *Cell*, [10.1016/j.cell.2011.12.012](https://doi.org/10.1016/j.cell.2011.12.012) | BrdU-IP-seq + RNA-seq, forkhead mutants | tabular processed data published; sacCer3 / R64 | Knott 2012 sequencing arm |
| GSE33702 | Knott SR 2012, *Cell*, [10.1016/j.cell.2011.12.012](https://doi.org/10.1016/j.cell.2011.12.012) | BrdU-IP-chip time course + HU, forkhead mutants | array probe-level only; microarray probe-level; coordinates set by platform GPL14887 | Knott 2012; NimbleGen HX12 CGH array |
| GSE52987 | Murakami H 2014, *Cell*, [10.1016/j.cell.2014.06.028](https://doi.org/10.1016/j.cell.2014.06.028) | copy-number replication profiles from total DNA, 8 time points x 6 strains | tabular processed data published; sacCer2 (2008 build) | ARS+/ars-delta, DDK overexpression, tof1-delta |
| GSE55130 | Peace JM 2014, *PLoS One*, [10.1371/journal.pone.0098501](https://doi.org/10.1371/journal.pone.0098501) | BrdU-IP-chip, rif1 and mec1-100 time course | array probe-level only; microarray probe-level; coordinates set by platform GPL14887 | Peace 2014; NimbleGen array arm |
| GSE55155 | Peace JM 2014, *PLoS One*, [10.1371/journal.pone.0098501](https://doi.org/10.1371/journal.pone.0098501) | BrdU-IP-seq in HU, WT and rif1-delta | tabular processed data published; sacCer3 / R64 | Peace 2014 sequencing arm |
| GSE89337 | Müller CA 2017, *J Cell Biol*, [10.1083/jcb.201701061](https://doi.org/10.1083/jcb.201701061) | sort-seq, 6 non-Saccharomyces yeasts + a cerevisiae reference | continuous signal track published; sacCer3 / R64 (custom reference for some samples) | Agier/Fischer 2018; cerevisiae sample usable as a control profile |
| GSE89611 | Shen Y 2017, *Science*, [10.1126/science.aaf4791](https://doi.org/10.1126/science.aaf4791) | sort-seq of a synthetic chromosome II strain | continuous signal track published; not stated in GEO | synII (Sc2.0); timing of a synthetic chromosome vs native |
| GSE90151 | Foss EJ 2017, *Proc Natl Acad Sci U S A*, [10.1073/pnas.1614781114](https://doi.org/10.1073/pnas.1614781114) | flow-sorted G1/S/early-G2/late-G2 sequencing | tabular processed data published; sacCer3 / R64 | Sir2 replication gaps; repetitive vs unique sequence balance |
| GSE94795 | Ostrow AZ 2017, *Proc Natl Acad Sci U S A*, [10.1073/pnas.1612422114](https://doi.org/10.1073/pnas.1612422114) | BrdU-IP-seq + RNA-seq + 4C, forkhead dimerization | interval/call set published; sacCer3 / R64 | Ostrow 2017 sequencing arm |
| GSE97953 | Hafner L 2018, *Cell Rep*, [10.1016/j.celrep.2018.03.113](https://doi.org/10.1016/j.celrep.2018.03.113) | Pol2 (leading-strand) ChIP-seq timing + copy number | continuous signal track published; sacCer3 / R64 | Rif1 and chromosome-internal origins |
| SRP264744 | Frenkel N 2021, *Genome Res*, [10.1101/gr.266510.120](https://doi.org/10.1101/gr.266510.120) | replication speed profiling | raw reads only; not stated in the archive record | Frenkel 2021 Genome Res; large run set, histone-acetylation effect on fork speed |

### initiation mapping (14)

| Accession | Study | Method | Current state | Notes |
|---|---|---|---|---|
| GSE103073 | Fang D 2017, *Genes Dev*, [10.1101/gad.306571.117](https://doi.org/10.1101/gad.306571.117) | copy-number origin activity in HU, Dbf4-forkhead recruitment | continuous signal track published; sacCer3 / R64 | Fang 2017; HU-arrested origin activity |
| GSE111058 | Soudet J 2018, *Genome Res*, [10.1101/gr.239582.118](https://doi.org/10.1101/gr.239582.118) | BrdU-seq + MNase-seq, transcription and the initiation programme | continuous signal track published; sacCer3 / R64 | Nrd1/Ysh1 anchor-away |
| GSE136602 | Forey R 2020, *Mol Cell*, [10.1016/j.molcel.2020.02.021](https://doi.org/10.1016/j.molcel.2020.02.021) | BrdU-IP-seq / ChIP-seq, origin activity and checkpoint | continuous signal track published; sacCer3 / R64 | Sequencing arm; larger and on sacCer3 |
| GSE155831 | Lee PH 2021, *Nucleic Acids Res*, [10.1093/nar/gkaa1148](https://doi.org/10.1093/nar/gkaa1148) | BrdU-IP-seq, quiescent (G0) vs G1 release | continuous signal track published; sacCer3 / R64 | Origin utilisation after quiescence re-entry |
| GSE166733 | Joshi I 2022, *PLoS One*, [10.1371/journal.pone.0263569](https://doi.org/10.1371/journal.pone.0263569) | BrdU/HU origin activation compared in two laboratory strains | tabular processed data published; microarray probe-level; coordinates set by platform GPL10930 | Strain-to-strain origin activation differences (W303 vs BY) |
| GSE173528 | Weekes C 2026, *Nat Commun*, [10.1038/s41467-026-75695-1](https://doi.org/10.1038/s41467-026-75695-1) | permanganate ssDNA PIP-seq, Mcm4-Flag | continuous signal track published; sacCer3 / R64 | DNA melting at origins |
| GSE186490 | *no publication identified* | copy-number + Mcm4/Sld2 ChIP through S phase | continuous signal track published; sacCer3 / R64 | PP2A-Cdc55/Pph3 dephosphorylation of the pre-IC |
| GSE199448 | Santos MM 2022, *Genome Biol*, [10.1186/s13059-022-02788-7](https://doi.org/10.1186/s13059-022-02788-7) | copy number after targeted origin deletions (ARS1008/1009, ARS816/818) | tabular processed data published; sacCer3 / R64 | Origin-deletion series; direct test of origin contribution |
| GSE200349 | He Y 2022, *Proc Natl Acad Sci U S A*, [10.1073/pnas.2212134119](https://doi.org/10.1073/pnas.2212134119) | BrdU-IP-seq, Rpd3 and Fkh1 origin stimulation | interval/call set published; sacCer3 / R64 | sacCer3; origin activity with rDNA context |
| GSE21014 | Crabbé L 2010, *Nat Struct Mol Biol*, [10.1038/nsmb.1932](https://doi.org/10.1038/nsmb.1932) | BrdU-IP-chip, WT + 20 mutants, Affymetrix tiling 1.0R | array probe-level only; microarray probe-level; coordinates set by platform GPL7250 | Crabbe 2010; large mutant panel of origin usage |
| GSE215190 | Petrie MV 2022, *Genes (Basel)*, [10.3390/genes13122202](https://doi.org/10.3390/genes13122202) | quantitative BrdU-IP-seq (QBU) + Dbf4 ChIP-seq | interval/call set published; SGD reference, version not stated | Ctf19-activated centromeric origins |
| GSE262693 | Wu Y 2024, *Nat Commun*, [10.1038/s41467-024-52408-0](https://doi.org/10.1038/s41467-024-52408-0) | syn-seq (synchronisation-sequencing) copy number, Orc2-NTD mutants | interval/call set published; sacCer3 / R64 | ORC intrinsically disordered region and licensing |
| GSE57617 | Yoshida K 2014, *Mol Cell*, [10.1016/j.molcel.2014.04.032](https://doi.org/10.1016/j.molcel.2014.04.032) | BrdU-IP-seq, wt/sir2/rpd3, genome + rDNA | continuous signal track published; sacCer3 / R64 | sacCer3-aligned; rDNA-aware processing |
| GSE61743 | Bogenschutz NL 2014, *PLoS One*, [10.1371/journal.pone.0114545](https://doi.org/10.1371/journal.pone.0114545) | BrdU/ChIP on an origin-depleted chromosome arm | interval/call set published; microarray probe-level; coordinates set by platform GPL10725 | Seven origins deleted on chrVI-L; non-canonical initiation |

### protein footprinting (10)

| Accession | Study | Method | Current state | Notes |
|---|---|---|---|---|
| GSE130273 | Foss EJ 2019, *PLoS Genet*, [10.1371/journal.pgen.1008138](https://doi.org/10.1371/journal.pgen.1008138) | MCM2-ChEC + MNase + replication profiling at the rDNA | tabular processed data published; sacCer3 / R64 | Sir2 and Mcm2-7 displacement in the rDNA |
| GSE147927 | Badjatia N 2021, *Cell Rep*, [10.1016/j.celrep.2020.108640](https://doi.org/10.1016/j.celrep.2020.108640) | ChIP-exo compendium of hundreds of yeast proteins | continuous signal track published; sacCer3 / R64 | Rossi 2021; includes ORC and MCM subunits among 1251 samples |
| GSE149163 | Lee CSK 2021, *Nat Commun*, [10.1038/s41467-020-20277-y](https://doi.org/10.1038/s41467-020-20277-y) | ORC/MCM ChIP-seq + BrdU-IP-seq, humanised ORC | continuous signal track published; sacCer2 (2008 build) | Chimeric ORC strains with WT controls |
| GSE173525 | Weekes C 2026, *Nat Commun*, [10.1038/s41467-026-75695-1](https://doi.org/10.1038/s41467-026-75695-1) | ChIP-exo 5.0 of Mcm4-Flag, G1 vs 30 min release | continuous signal track published; sacCer3 / R64 | Mcm2/5 gate mutant |
| GSE240472 | Reuter LM 2024, *Nat Commun*, [10.1038/s41467-024-51538-9](https://doi.org/10.1038/s41467-024-51538-9) | ORC/MCM ChIP-seq, cdc6-1 and cdc46-1, pre-RC formation | continuous signal track published; sacCer3 / R64 | Companion to the ChIP-exo series |
| GSE285768 | Lichauco C 2025, *Elife*, [10.7554/eLife.97438](https://doi.org/10.7554/eLife.97438) | MCM ChEC/ChIP + MNase + timing, Sir2 and Fun30, rDNA | continuous signal track published; sacCer3 / R64 | 168 samples; rDNA replication timing via MCM positioning |
| GSE52614 | Yu C 2014, *Mol Cell*, [10.1016/j.molcel.2014.09.017](https://doi.org/10.1016/j.molcel.2014.09.017) | BrdU-IP + replisome protein ChIP-seq in HU | continuous signal track published; sacCer3 / R64 | Genome-wide map of proteins at forks |
| GSE69065 | Gros J 2015, *Mol Cell*, [10.1016/j.molcel.2015.10.022](https://doi.org/10.1016/j.molcel.2015.10.022) | Mcm2-7 ChIP-seq, post-licensing sliding | continuous signal track published; sacCer3 / R64 | Belsky 2015; MCM positions distinct from ORC sites |
| SRP242629 | Hoggard T 2020, *Proc Natl Acad Sci U S A*, [10.1073/pnas.2004664117](https://doi.org/10.1073/pnas.2004664117) | ORC and MCM ChIP-seq + whole-genome sequencing | raw reads only; not stated in the archive record | Hoggard 2020 PNAS; licensing-efficiency comparison, deposited to SRA without a GEO series |
| SRP489912 | Hoggard T 2024, *PLoS Genet*, [10.1371/journal.pgen.1011366](https://doi.org/10.1371/journal.pgen.1011366) | ORC/MCM and chromatin profiling | raw reads only; not stated in the archive record | Hoggard 2024 PLoS Genet; SRA-only deposit |

### fork direction / polymerase usage (7)

| Accession | Study | Method | Current state | Notes |
|---|---|---|---|---|
| GSE101698 | Garbacz MA 2018, *Nat Commun*, [10.1038/s41467-018-03270-4](https://doi.org/10.1038/s41467-018-03270-4) | ribonucleotide mapping, Pol delta at leading-strand initiation | continuous signal track published; not stated in GEO | Steric-gate variants |
| GSE115897 | Porcella SY 2020, *PLoS Genet*, [10.1371/journal.pgen.1008755](https://doi.org/10.1371/journal.pgen.1008755) | Okazaki fragment sequencing, Pol alpha recruitment | interval/call set published; sacCer3 / R64 | Ligase-depletion enrichment; 70 samples |
| GSE118078 | Kahli M 2019, *Nucleic Acids Res*, [10.1093/nar/gky1242](https://doi.org/10.1093/nar/gky1242) | Okazaki fragment sequencing, nuclease redundancy | tabular processed data published; sacCer3 / R64 | Fragment termini and maturation |
| GSE141884 | *no publication identified* | Okazaki fragments + ribonucleotides, Pol3 degron | tabular processed data published; sacCer3 / R64 | Lagging-strand processing and fork dynamics |
| GSE42929 | Yanga W 2013, *FEBS Lett*, [10.1016/j.febslet.2013.06.014](https://doi.org/10.1016/j.febslet.2013.06.014) | Okazaki fragment sequencing, S288C | interval/call set published; sacCer3 / R64 | Single-sample deposit of Okazaki fragments |
| GSE64521 | Reijns MAM 2015, *Nature*, [10.1038/nature14183](https://doi.org/10.1038/nature14183) | ribonucleotide mapping of replicative polymerase contributions | interval/call set published; sacCer3 / R64 | Polymerase usage in the fully replicated genome |
| GSE89721 | Gan H 2017, *Mol Cell*, [10.1016/j.molcel.2017.09.018](https://doi.org/10.1016/j.molcel.2017.09.018) | eSPAN / BrdU-IP-seq, Rad53 and strand coupling | continuous signal track published; sacCer3 / R64 | Strand-specific nascent DNA |

### fork position, pausing, termination (8)

| Accession | Study | Method | Current state | Notes |
|---|---|---|---|---|
| GSE139860 | *no publication identified* | Okazaki fragment sequencing, fork arrest at tDNAs | interval/call set published; sacCer3 / R64 | Depletion panel; fork pausing quantification |
| GSE140034 | Shyian M 2020, *Genes Dev*, [10.1101/gad.331868.119](https://doi.org/10.1101/gad.331868.119) | Top1 and replisome (Cdc45, Mcm4, Pol2) ChIP-seq | continuous signal track published; sacCer3 / R64 | Fork pausing complex and topoisomerases |
| GSE144321 | Westhorpe R 2020, *Nucleic Acids Res*, [10.1093/nar/gkaa963](https://doi.org/10.1093/nar/gkaa963) | replisome ChIP-seq, Tof1/Mrc1 uncoupling | continuous signal track published; sacCer3 / R64 | Helicase-polymerase uncoupling; replisome positions |
| GSE16218 | Azvolinsky A 2009, *Mol Cell*, [10.1016/j.molcel.2009.05.022](https://doi.org/10.1016/j.molcel.2009.05.022) | Pol2 and Rrm3 ChIP-chip, fork pause sites | tabular processed data published; microarray probe-level; coordinates set by platform GPL8192 | Azvolinsky 2009; classic pause-site map (PCR array) |
| GSE19061 | Fachinetti D 2010, *Mol Cell*, [10.1016/j.molcel.2010.07.024](https://doi.org/10.1016/j.molcel.2010.07.024) | BrdU-IP-chip + ChIP-chip, Top2-mediated termination | array probe-level only; microarray probe-level; coordinates set by platform GPL7250 | One of the few explicit termination-site datasets |
| GSE214930 | Choudhary R 2023, *Cell Rep*, [10.1016/j.celrep.2023.112747](https://doi.org/10.1016/j.celrep.2023.112747) | SuperSeries: Sen1/Rrm3 and replication termination | continuous signal track published; sacCer3 / R64 | ChIP-seq and ChIP-chip arms |
| GSE235881 | Hadjicharalambous A 2023, *Wellcome Open Res*, [10.12688/wellcomeopenres.19617.1](https://doi.org/10.12688/wellcomeopenres.19617.1) | TrAEL-seq, fork locations in HU and recovery | tabular processed data published; sacCer3 / R64 | Replisome mutants; fork recovery time course |
| GSE71973 | Osmundson JS 2017, *Nat Struct Mol Biol*, [10.1038/nsmb.3342](https://doi.org/10.1038/nsmb.3342) | Okazaki fragment sequencing, Pif1-family helicases at tRNA genes | continuous signal track published; sacCer3 / R64 | Fork arrest quantified by Watson:Crick ratio |

### re-replication (3)

| Accession | Study | Method | Current state | Notes |
|---|---|---|---|---|
| GSE143572 | Menzel J 2020, *Nucleic Acids Res*, [10.1093/nar/gkaa197](https://doi.org/10.1093/nar/gkaa197) | Rerep-seq, enrichment of re-replicated DNA | continuous signal track published; sacCer3 / R64 | Method paper for re-replication mapping |
| GSE4181 | Green BM 2006, *Mol Biol Cell*, [10.1091/mbc.e05-11-1043](https://doi.org/10.1091/mbc.e05-11-1043) | array CGH of re-replication | raw reads only; microarray probe-level; coordinates set by platform GPL3412 | Green 2006; first genome-wide re-replication map |
| GSE4487 | Tanny RE 2006, *Mol Biol Cell*, [10.1091/mbc.e05-11-1037](https://doi.org/10.1091/mbc.e05-11-1037) | array CGH of re-replication, inhibitory controls | raw reads only; microarray probe-level; coordinates set by platform GPL3499/GPL884 | Green & Li 2006 |

### ssDNA / replication stress (1)

| Accession | Study | Method | Current state | Notes |
|---|---|---|---|---|
| GSE4099 | Feng W 2006, *Nat Cell Biol*, [10.1038/ncb1358](https://doi.org/10.1038/ncb1358) | ssDNA microarray in HU, cerevisiae and pombe | raw reads only; microarray probe-level; coordinates set by platform GPL3388/GPL3398 | Feng 2006; classic ssDNA origin-firing map (legacy array) |

## Priority 3 — peripheral datasets

77 further track-hub candidates are catalogued at priority 3 and listed in full in the TSV:

| Data type | Datasets |
|---|---|
| replication timing | 16 |
| initiation mapping | 17 |
| protein footprinting | 21 |
| fork direction / polymerase usage | 9 |
| fork position, pausing, termination | 4 |
| re-replication | 7 |
| ssDNA / replication stress | 3 |

The common reasons for a priority-3 call are: a legacy microarray platform whose coordinates must be re-derived from the array design; a dataset confined to one locus or one chromosome; alignment to a strain-specific or *de novo* assembly rather than sacCer3; or a SuperSeries record that only points at component series already catalogued.


## Single-molecule datasets — a separate serving route

These 9 datasets are deliberately held apart from the track-hub list. Their content is per-read: individual molecules carrying BrdU incorporation tracts, fork positions, fork directions and velocities, or per-molecule chromatin states. UCSC and comparable browsers have no useful rendering for this — collapsing it to a coverage track discards exactly the information the assay was done to obtain. They are intended for the group's own single-molecule visualisation tool.

Two of the nine (ERP161317, ERP172809) also produce a conventional genome-wide profile as their published product, so they can feed both routes: the aggregate profile as a hub track, the per-molecule data to the viewer. Their `serving_route` value records this.

This is also the group with the heaviest data cost: none of the ENA deposits publishes processed output, so the per-molecule calls must be regenerated with the authors' own software (DNAscent, RepNano/NanoForkSpeed) from raw signal.

| Accession | Priority | Study | Method | Serving route | Current state | Notes |
|---|---|---|---|---|---|---|
| ERP120013 | 1 | Hennion M 2020, *Genome Biol*, [10.1186/s13059-020-02013-3](https://doi.org/10.1186/s13059-020-02013-3) | nanopore BrdU single-molecule (FORK-seq) | single-molecule viewer | raw reads only; not stated in the archive record | Hennion 2020; per-molecule fork direction and initiation/termination. Processed RFD tracks are in the paper's supplement, not the archive |
| ERP134869 | 1 | Theulot B 2022, *Nat Commun*, [10.1038/s41467-022-31012-0](https://doi.org/10.1038/s41467-022-31012-0) | nanopore BrdU single-molecule (NanoForkSpeed) | single-molecule viewer | raw reads only; not stated in the archive record | Theulot 2022; fork velocity distributions genome-wide |
| ERP161317 | 1 | Theulot B 2025, *Nat Commun*, [10.1038/s41467-024-55520-3](https://doi.org/10.1038/s41467-024-55520-3) | nanopore BrdU single-molecule (Nanotiming) | both: aggregate profile to hub, per-molecule to viewer | raw reads only; not stated in the archive record | Theulot 2025; timing including subtelomeres that short-read profiles cannot resolve |
| GSE121941 | 1 | Müller CA 2019, *Nat Methods*, [10.1038/s41592-019-0394-y](https://doi.org/10.1038/s41592-019-0394-y) | nanopore BrdU detection on ultra-long reads (D-NAscent) | single-molecule viewer | interval/call set published; sacCer3 / R64 | Muller 2019; the founding DNAscent dataset |
| ERP172809 | 2 | Pellet J 2025, *Nucleic Acids Res*, [10.1093/nar/gkaf754](https://doi.org/10.1093/nar/gkaf754) | nanopore BrdU single-molecule | both: aggregate profile to hub, per-molecule to viewer | raw reads only; not stated in the archive record | Pellet 2025; single-chromosome strain vs 16-chromosome reference |
| GSE194025 | 2 | Claussin C 2022, *Mol Cell*, [10.1016/j.molcel.2022.02.010](https://doi.org/10.1016/j.molcel.2022.02.010) | nanopore replisome mapping (Mcm4-MNase cleavage) | single-molecule viewer | interval/call set published; custom / non-S288C reference | Single-molecule replisome progression |
| GSE220757 | 2 | Chanou A 2023, *Nucleic Acids Res*, [10.1093/nar/gkad1022](https://doi.org/10.1093/nar/gkad1022) | MATAC-seq single-molecule methyl-accessibility at origins | single-molecule viewer | tabular processed data published; sacCer3 / R64 | Chromatin states linked to origin efficiency |
| ERP110889 | 3 | *no publication identified* | nanopore BrdU single-molecule | single-molecule viewer | raw reads only; not stated in the archive record | IBENS deposit including E. coli controls; no separate publication identified - probably the precursor deposit to FORK-seq |
| GSE212101 | 3 | Meng X 2023, *Genes Dev*, [10.1101/gad.350054.122](https://doi.org/10.1101/gad.350054.122) | nanopore replisome mapping, Pol epsilon domain mutants | single-molecule viewer | interval/call set published; custom / non-S288C reference | Aligned to a Canu custom assembly |

## Datasets held outside GEO

Eleven in-scope deposits exist in ENA/SRA with no GEO series. Five of them are nanopore studies listed in the single-molecule section above (ERP110889, ERP120013, ERP134869, ERP161317, ERP172809); the remaining 6 are track-hub candidates:

| Accession | Priority | Study | Method | Runs | Notes |
|---|---|---|---|---|---|
| SRP242629 | 2 | Hoggard T 2020, *Proc Natl Acad Sci U S A*, [10.1073/pnas.2004664117](https://doi.org/10.1073/pnas.2004664117) | ORC and MCM ChIP-seq + whole-genome sequencing | 17 | Hoggard 2020 PNAS; licensing-efficiency comparison, deposited to SRA without a GEO series |
| SRP264744 | 2 | Frenkel N 2021, *Genome Res*, [10.1101/gr.266510.120](https://doi.org/10.1101/gr.266510.120) | replication speed profiling | 448 | Frenkel 2021 Genome Res; large run set, histone-acetylation effect on fork speed |
| SRP489912 | 2 | Hoggard T 2024, *PLoS Genet*, [10.1371/journal.pgen.1011366](https://doi.org/10.1371/journal.pgen.1011366) | ORC/MCM and chromatin profiling | 20 | Hoggard 2024 PLoS Genet; SRA-only deposit |
| SRP040498 | 3 | *no publication identified* | MNase-seq for MCM distribution | 2 | W303; no publication identified from the archive record alone |
| SRP253620 | 3 | Xu P 2021, *Nucleic Acids Res*, [10.1093/nar/gkab801](https://doi.org/10.1093/nar/gkab801) | ribonucleotide mapping around ARS elements | 34 | Probable match to Xu 2021 NAR - confirm before citing |
| SRP469961 | 3 | Bugallo A 2024, *Nucleic Acids Res*, [10.1093/nar/gkae192](https://doi.org/10.1093/nar/gkae192) | Okazaki/nascent strand sequencing | 69 | Bugallo 2024 NAR |

## Genome-wide studies with no archive deposit

Several studies that OriDB already cites predate routine data deposition. Their data exists only as journal supplementary material, and their origin calls are already represented in the OriDB tables. Each was checked for a GEO link directly rather than assumed; none has one.

| Study | DOI | Content | Data availability |
|---|---|---|---|
| Raghuraman MK 2001, *Science* | [10.1126/science.294.5540.115](https://doi.org/10.1126/science.294.5540.115) | Density-transfer microarray replication timing — the original genome-wide timing profile | journal supplement only |
| Wyrick JJ 2001, *Science* | [10.1126/science.1066101](https://doi.org/10.1126/science.1066101) | ORC and MCM ChIP-chip — the first genome-wide licensing map | journal supplement only |
| Yabuki N 2002, *Genes Cells* | [10.1046/j.1365-2443.2002.00559.x](https://doi.org/10.1046/j.1365-2443.2002.00559.x) | Replication profile and early-firing origins | journal supplement only |
| Xu W 2006, *BMC Genomics* | [10.1186/1471-2164-7-276](https://doi.org/10.1186/1471-2164-7-276) | ORC and Mcm2 binding on tiling arrays | journal supplement only |
| Alvino GM 2007, *Mol Cell Biol* | [10.1128/MCB.00719-07](https://doi.org/10.1128/MCB.00719-07) | Replication timing in hydroxyurea | journal supplement only |

## Suggested ingestion order

Ordered by value per unit of work, taking the hub's existing Müller 2014 dataset as the template and the agreed policy of using author-processed data first.

1. **GSE242131** (Foss 2024, MCM ChEC-seq). 17,618 single-base MCM binding sites on sacCer3, already a call set — this becomes a bigBed with no reprocessing and is the highest-resolution origin evidence available.
2. **GSE240779 / GSE173527** (Reuter 2024, ORC and MCM ChIP-exo). sacCer3 signal tracks; gives an independent licensing footprint beside Foss, with ORC and MCM separated.
3. **GSE249817** (Chappleboim 2024, ORC ChEC-seq). sacCer3; wild-type plus a large ORC mutant panel, so one ingestion yields both a reference track and a comparative resource.
4. **GSE40696** (McGuffee 2013, Ok-seq). The fork-direction axis. Requires a sacCer2→sacCer3 lift of the published .sgr tracks, and the OEM derivation is reproducible from them.
5. **GSE48561 and GSE36045** (Hawkins 2013; Müller 2012, sort-seq). Timing data from the same assay family as the dataset already in the hub, so the display configuration carries over directly.
6. **GSE16926** (Eaton 2010, ORC ChIP-seq and nucleosome positions). The continuity anchor with OriDB v2's ACS element layer; coordinates need confirming before publication.

Items 1–3 are all sacCer3 with processed files already published, so the first three steps should be close to configuration-only work.

The single-molecule datasets are deliberately absent from this order: they are not blocked on hub work, and their sequencing follows whatever schedule the separate viewer sets. Within that group, GSE121941, ERP120013 and ERP134869 are the three to start from.


## Caveats

- **Assembly claims are author statements, not verifications.** 109 datasets state a build in their GEO processing text; 55 do not, and are marked accordingly.
- **8 in-scope datasets have no publication I could identify**, and 27 are linked to a publication by title match rather than by a PMID recorded in GEO. The `citation_link` column marks which is which; the one `probable match` (GSE115897) should be confirmed before it is cited.
- **The file inventory counts what the archive lists**, including per-sample files. Some series bundle processed files inside a `*_RAW.tar`, in which case the readiness column reports what is individually listed.
- **SuperSeries are catalogued alongside their component series.** They are flagged in the notes and should not be ingested separately.
- **Coverage is complete for GEO and for ENA/SRA study records under the terms used.** Datasets published only as journal supplementary tables, or held in lab-specific repositories, are not systematically represented beyond the legacy studies listed above.
