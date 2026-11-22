# Candidate *Schizosaccharomyces pombe* replication datasets for the OriDB data hub

Survey date: 21 November 2026. Compiled by querying NCBI GEO (`db=gds`), NCBI SRA and BioProject
(E-utilities), EBI BioStudies/ArrayExpress and PubMed. Every accession below was retrieved from a live
archive record in this task and its GEO SOFT header read; no accession is quoted from memory. No data
files were downloaded, so statements about file *contents* (as opposed to file names, formats and the
archive's own `Genome_build` declarations) are explicitly flagged where they are inferences.

## 1. What the pombe replication data landscape looks like

The fission yeast replication literature splits cleanly into two eras, and the archive record reflects
that split. The microarray era (2005-2012) produced the datasets that define the field's origin
nomenclature - ssDNA mapping under hydroxyurea (Feng 2006), genome-wide origin characterisation
(Heichinger 2006), pre-RC/ORC-MCM tiling arrays (Hayashi 2007), two-colour copy-number timing (Eshaghi
2007), the origin-firing programme (Wu and Nurse 2009) and the Masai laboratory's BrdU-IP tiling arrays
on Mrc1 and Rif1 (2011-2012). Those datasets are deposited as raw CEL, GPR or Agilent feature-extraction
text, carry no stated genome build, and in several cases their probes are indexed by ORF name rather
than by coordinate, so they are coordinate-bearing only to the extent that the GEO platform record
supplies an annotation table. The sequencing era (2015-2025) is dominated by polymerase-usage methods,
which are a genuine pombe speciality: Pu-seq from the Carr laboratory at Sussex (Daigaku 2015 and at
least four later applications), HydEn-seq from the Kunkel laboratory (Zhou 2019), and more recently
strand-resolved methods borrowed from the chromatin field - SCAR-seq (Charlton 2024) and GLOE-seq
(Martienssen laboratory, deposited 2024, released 2025). Alongside these sit short-read BrdU-IP-seq
timing profiles (Daigaku 2017; Kanoh 2015; Hasegawa and Kanoh, unpublished) and ChIP-seq or ChEC-seq
maps of the licensing machinery (Masuda 2020; Foss 2021). Assembly practice in the sequencing era is
near-uniform and favourable to the hub: the archive states ASM294v2 - or a point release of it
(ASM294v2.19, ASM294v2.23, "assembly 2.23") - for almost every sequencing dataset examined, with two
exceptions noted in the table. Two categories of the hub's controlled vocabulary are empty for pombe:
I found no re-replication dataset, and no single-molecule (nanopore or PacBio) replication dataset. For
the latter I scanned all 143 public Oxford Nanopore SRA runs with organism *S. pombe* and found only
genome-assembly and RNA work; the one run whose strain genotype resembles a Pu-seq strain
(`rnh201 cdc6-L591G`, SRR7193431 / SRR7193498) belongs to study SRP148618, a *Schizosaccharomyces*
genome sequencing and assembly project from Edinburgh, not to a replication experiment.

Two datasets that the literature would lead one to expect are **not** listed below because I could not
verify an archive accession for them. Heichinger et al. 2006 (EMBO J, PMID 17053780,
doi 10.1038/sj.emboj.7601390) has no GEO or ArrayExpress accession that surfaced in BioStudies,
ArrayExpress or GEO searches, and the publisher full text available to me was abstract-only, so I could
not read its data-availability statement. Xu et al. 2012 (Genome Biol, PMID 22531001,
doi 10.1186/gb-2012-13-4-r27) states that filtered reads were submitted to the Sequence Read Archive
with the accessions listed in its Additional file 2, Table S8; that table is not part of the PMC full
text I retrieved, `elink` returns no pubmed-to-SRA link for the paper, and I therefore could not
establish the accession. Both would be worth chasing by hand - the Xu profiles were, per that paper,
submitted to OriDB itself.

## 2. Candidate table

Fifteen candidates, one row each. `n_samples` is the archive's sample count; where a series is
multi-species the pombe-only count is given as well. `processed_files` lists what the archive actually
holds, taken from the GEO supplementary-file inventory.

| accession | repository | data_type | method | title | year | first_author | journal | pmid | doi | platform_type | instrument | n_samples | assembly | assembly_source | processed_files | readiness | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GSE62108 | GEO | fork direction / polymerase usage | Pu-seq: ribonucleotide tracking in pol-delta-L591G and pol-epsilon-M630F, rnh201 background | A global profile of replication polymerase usage | 2015 | Daigaku | Nat Struct Mol Biol | 25664722 | 10.1038/nsmb.2962 | short-read sequencing | Illumina HiSeq 2000 | 11 (all pombe) | ASM294v2 | GEO sample `Genome_build` (GSM1519714) | series-level `PU-seq_Trep.wig`, `PU-seq_leftward_moving_fork.wig`, `PU-seq_terminaton-events.wig`, `origin-effici.bedgraph`; 25 per-sample strand-resolved wig files | high - four derived genome-wide profiles already in wig/bedGraph | The reference pombe polymerase-usage dataset; later Pu-seq papers normalise against it. Chromosome naming inside the wig files was not inspected (no files downloaded). Charlton 2024 deposited `GSE62108_origin_20plus_chromosome_liftover.bed` in GSE241062, which indicates another group needed to re-map these origin coordinates. |
| GSE70033 | GEO | replication timing | BrdU-IP-seq in synchronised cells, converted to local replication rate at 25/50/75 % genome replicated | The role of PCNA ubiquitylation during unperturbed genome replication in S. pombe | 2017 | Daigaku | PLoS Genet | 28481910 | 10.1371/journal.pgen.1006789 | short-read sequencing | NextSeq 500 | 14 (all pombe) | ASM294v2.23 | GEO sample `Genome_build` (GSM1715906) | 6 series-level wig files (`BrdU-IP-seq-cdc25-wt-gr25pc/50pc/75pc` and matching mutant) | high - genome-wide wig, wild-type series included | Same laboratory and binning convention (300 bp) as GSE62108, so the timing and fork-direction tracks will register against each other without extra work. |
| GSE125335 | GEO | protein footprinting | Orc4 and Mcm2 ChIP-seq in G1, fold-enrichment profiles produced with DROMPA2 | Genome-wide map of DNA replication origin in Schizosaccharomyces pombe | 2020 | Masuda | Open Biol | 32692956 | 10.1098/rsob.200052 | short-read sequencing | Illumina HiSeq 2500 | 6 (all pombe) | ASM294v2.19 | GEO sample `Genome_build` (GSM3569251) | per-sample wig, one file per chromosome (`..._Orc4_G1_chr1.wig`, `_chr2`, `_chr3`), plus one bedGraph | high - wig, but split per chromosome and using `chr1/chr2/chr3` naming | Gives the licensing axis: ORC and MCM occupancy against which OriDB's 741 pombe origin records can be read. Three files per sample must be concatenated and renamed to I/II/III. |
| GSE65293 | GEO | initiation mapping | BrdU-IP-seq under HU plus Rif1 and Mcm4 ChIP-seq; Rif1 binding-site mutants | Rif1 binds to G quadruplexes and suppresses replication over long distances | 2015 | Kanoh | Nat Struct Mol Biol | 26436827 | 10.1038/nsmb.3102 | short-read sequencing | Illumina MiSeq; Illumina Genome Analyzer IIx | 15 (all pombe) | ASM294v2 | GEO sample `Genome_build` (GSM1591762) | 9 series-level bedGraph files (BrdU-IP HU wild-type, rif1-delta, two Rif1 binding-site mutants; Rif1 ChIP in wild-type, rap1-delta, taz1-delta; Mcm4 ChIP) | high - named bedGraph per condition, no tar unpacking needed | Carries both an early-origin firing profile and the trans-acting regulator that shapes it, which is the single most-requested comparison in pombe timing. |
| GSE93178 | GEO | initiation mapping | Mcm4 ChIP-seq in G1 and BrdU-IP-seq under HU in hsk1 and mrc1 backgrounds | Checkpoint-independent regulation of origin firing by Mrc1 through interaction with Hsk1 kinase | 2017 | Matsumoto | Mol Cell Biol | 28069740 | 10.1128/MCB.00355-16 | short-read sequencing | Illumina Genome Analyzer IIx; Illumina MiSeq | 7 (all pombe) | ASM294v2 | GEO sample `Genome_build` (GSM2445931) | 5 per-sample bedGraph files inside `GSE93178_RAW.tar` | medium-high - bedGraph, but inside a RAW tar | Same pipeline and build as GSE65293; a natural second track in a Masai-laboratory initiation section. |
| GSE94812 | GEO | replication timing | BrdU-DNA-seq of newly replicated DNA under HU in telomere-binding-protein deletions | Genome-wide replication timing in the mutants deleted for telomere-binding proteins in S. pombe | 2017 (submitted); released 2019 | Hasegawa (GEO contributor; no linked publication) | not established | none linked in GEO | not established | short-read sequencing | Illumina MiSeq | 8 (all pombe) | ASM294v2 | GEO sample `Genome_build` (GSM2483550) | 48 per-sample wig files, one per sequence, named `..._I.wig`, `_II.wig`, `_III.wig`, `_MT.wig`, `_AB325691.wig` | medium-high - wig already using PomBase chromosome names I/II/III | Linked to BioProject PRJNA374497 and SRA study SRP099408. No publication is linked in the GEO record, so provenance for a public hub track is weaker than for the rows above; worth asking Kanoh before serving. |
| GSE98462 | GEO | initiation mapping | S-phase versus unreplicated DNA ratio on Agilent arrays (origin efficiency), with ssDNA and Rad52 ChIP-chip sub-series | Genome-wide profiles of replication initiation, single-stranded DNA accumulation, and Rad52 binding in cells exposed to replication stress (SuperSeries) | 2018 | Gómez-Escoda | Genome Res | 29934426 | 10.1101/gr.224527.117 | two-colour microarray (Agilent 4x44K, GPL16235) | not applicable | 72 (all pombe) | not stated by the submitter; platform probes are annotated `chrI:162-221`-style | GPL16235 platform annotation table (`CHROMOSOMAL_LOCATION`) | 72 per-sample text files with a baseline-corrected origin-efficiency `Value` column | medium - coordinates exist but only via a probe-to-coordinate join | A second laboratory's origin-efficiency measurement, with a perturbation panel. Probe spacing on a 44K array over a 12.6 Mb genome is coarse (roughly 1 probe per 290 bp nominally, but the array is a ChIP-on-chip design, so effective spacing was not verified). Sub-series: GSE98444, GSE98445, GSE98446, GSE107097, GSE107099, GSE113559. |
| GSE125855 | GEO | fork direction / polymerase usage | HydEn-seq: alkaline cleavage at genomic ribonucleotides, 5'-end mapping, in rnh201 polymerase-mutant strains | Mapping polymerase usage by tracking ribonucleotide incorporation (HydEn-seq) during DNA replication | 2019 | Zhou | Nat Commun | 31488849 | 10.1038/s41467-019-11995-z | short-read sequencing | Illumina HiSeq 4000 | 58 total, 24 pombe | S. pombe assembly 2.23 (cerevisiae samples on a custom L03 assembly) | GEO sample `Genome_build` (GSM3583288) | 116 per-sample bigWig files, forward and reverse strand, normalised 5'-end counts | high on format (bigWig), medium on curation (multi-species series; pombe files must be selected by sample) | Methodologically independent of Pu-seq and therefore the best available cross-check on it; the Kunkel laboratory's leading/lagging assignment differs in detail from the Sussex one, which is itself worth showing in a browser. |
| GSE278836 | GEO | fork direction / polymerase usage | GLOE-seq: strand-specific mapping of DNA nicks and breaks in wild-type, dcr1, rnh and triple mutants | Transcription-replication conflict resolution by nuclear RNA interference [GLOE-seq] | 2024 (submitted); released 2025 | Cheng (GEO contributor; no linked publication) | not established | none linked in GEO | not established | short-read sequencing | NextSeq 500 | 8 (all pombe) | ASM294v2 | GEO sample `Assembly` field (GSM8555672) | 8 series-level bigWig files: `wt`, `dcr`, `rnh`, `tri`, each forward and reverse | highest on format - bigWig on ASM294v2, nothing to convert | Sub-series of GSE278850 (Martienssen laboratory). No publication linked at the time of retrieval, so the interpretation of the tracks is not yet citable; the files themselves are hub-ready today. |
| GSE241062 | GEO | fork direction / polymerase usage | SCAR-seq: strand-specific nascent-chromatin sequencing, giving histone partition and replication fork directionality | The fork protection complex promotes parental histone recycling and epigenetic maintenance | 2024 | Charlton | Cell | 39094569 | 10.1016/j.cell.2024.07.017 | short-read sequencing | NextSeq 2000 | 70 total, 58 pombe | "S.pombe h90 contig" | GEO sample `Assembly` field (GSM7714584) | 62 per-sample bigWig, including explicit `..._RFD.bw` fork-directionality tracks; plus `GSE241062_GSE62108_origin_20plus_chromosome_liftover.bed` | medium - bigWig and RFD ready, but on a non-standard assembly | The only pombe RFD track derived from a chromatin method rather than a polymerase method. The stated assembly is a laboratory h90 contig set, not ASM294v2, so the bigWigs cannot be trusted on ASM294v2 without checking the contig-to-chromosome relationship, and may require realignment from reads. |
| GSE153731 | GEO | fork position, pausing, termination | Pu-seq around an inducible RTS1 replication fork barrier and a 10x rDNA barrier, recombination-dependent replication | Replication dynamics of recombination-dependent replication forks | 2021 | Naiman | Nat Commun | 33568651 | 10.1038/s41467-021-21198-0 | short-read sequencing | NextSeq 500 | 36 (all pombe) | PomBase reference sequence with inserted RTS1 and 10x rDNA barrier constructs | GEO sample `Genome_build` (GSM4650274) | 72 per-sample CSV files, 300 bp binned strand-specific counts | medium - coordinate-bearing but on a modified genome, and CSV not bedGraph | The best available pombe fork-pausing and termination dataset. Because the reference carries engineered inserts, any track served on stock ASM294v2 must either drop the insert region or be served with an explicit caveat. |
| GSE6523 | GEO | protein footprinting | ORC and MCM (pre-RC) ChIP-chip on Affymetrix tiling arrays | Genome-wide localization of ORC and MCM binding sites and identification of replication origins in S. pombe | 2007 | Hayashi | EMBO J | 17304213 | 10.1038/sj.emboj.7601585 | tiling microarray (Affymetrix S. pombe Tiling 1.0FR, GPL3847) | not applicable | 10 (all pombe) | not stated by the submitter; platform probes carry RefSeq accessions `NC_003424.2` etc. with START/END | GPL3847 platform annotation table | raw CEL only (10 files) | low - raw CEL, requires full reprocessing | Historically the defining pre-RC map for pombe and the source of many OriDB origin calls. Ingestion means CEL normalisation plus a probe-coordinate join; the probe coordinates are on RefSeq chromosome accessions (CU329670/CU329671/CU329672 equivalents), which maps onto ASM294v2 one-to-one. |
| GSE6977 | GEO | replication timing | Two-colour array copy-number timing through S phase, HU block and release | Genome-wide profiling of DNA replication timing in pombe | 2007 | Eshaghi | PLoS One | 17684567 | 10.1371/journal.pone.0000722 | spotted long-oligo microarray (10K, GPL1932) | not applicable | 70 (all pombe) | not stated by the submitter; platform table has no coordinate columns, probes indexed by ORF (`SPBC1348.01`) | GPL1932 platform annotation table | raw GenePix `.gpr` only (70 files) | low - raw two-colour scans, and probes carry no coordinates | The classic pombe timing profile, cited widely, but reconstructing a browser track means re-deriving probe positions from PomBase systematic IDs and then refitting the sigmoid timing model. Worth listing as the historical reference, not worth ingesting early. |
| GSE4099 | GEO | ssDNA / replication stress | ssDNA mapping in hydroxyurea (S/G1 Cy3:Cy5 ratio) on a spotted array, in both budding and fission yeast | Single-stranded DNA formation during S phase in the presence of hydroxyurea in S. cerevisiae and S. pombe | 2006 | Feng | Nat Cell Biol | 16429127 | 10.1038/ncb1358 | spotted microarray (Eurogentec AR-ScPGS, GPL3398) | not applicable | 10 total, 2 pombe | not stated; platform table supplies `chr` (1/2/3) and `coord (kb)` | GPL3398 platform annotation table | none - no supplementary files; the ratio values sit in the GEO sample data tables | low - values must be pulled from sample tables and joined to kb-resolution probe coordinates | The only genome-wide pombe ssDNA dataset found, and so the only candidate for the ssDNA/replication-stress section. Coordinates are given in kilobases on a 2005-era sequence; agreement with ASM294v2 was not verified. |
| GSE150800 | GEO | protein footprinting | MCM2-ChEC (MNase tethered to Mcm2), G1 and HU, compared across five species | Chromosomal Mcm2-7 distribution and the genome replication program in species from yeast to humans | 2021 | Foss | PLoS Genet | 34473702 | 10.1371/journal.pgen.1009714 | short-read sequencing | Illumina HiSeq 2500 | 22 total, 4 pombe | Schizosaccharomyces_pombe_ASM294v2 | GEO sample `Genome_build` (GSM5536939) | pombe-specific processed text files of read depths (`GSE150800_pombe_MCM2ChEC_HU_G1_read_depths...`, `processed_pombe_MCM2ChEC_HU_rep2...`, `processed_pombe_free_MNase_rep1...`); bedGraph exists for the mouse and human arms only | medium - correct assembly, but pombe arm delivered as text read-depth tables, not bedGraph | Only four pombe samples, and the free-MNase control is the one that matters for interpreting the signal. Valuable as a licensing cross-check on GSE125335 using an orthogonal (ChEC rather than ChIP) chemistry. |

### Datasets examined and deliberately not tabled

Four further Pu-seq applications from the same laboratory are available and would make a coherent
comparative panel later, but each is a perturbation of the GSE62108 measurement rather than a new axis:
GSE113747 (Pai 2019, J Cell Sci, PMID 30674555, bedGraph origin positions and efficiencies),
GSE165503 (Zach 2021, Wellcome Open Res, PMID 33796794, ASM294v2, CSV counts),
GSE192342 (Rtf2 and RTS1, no publication linked, CSV counts on a modified reference) and
GSE247371 (SUMO protease at arrested forks, CSV). Three array-era series from the Masai laboratory
(GSE28182, Hayano 2011, PMID 21518960; GSE34369, Hayano 2012, PMID 22279046; GSE28879, de Castro 2012,
PMID 21989386) sit on Affymetrix platform GPL7715, whose GEO record carries no annotation table at all,
so probe coordinates would have to come from outside GEO; GSE28879 is the exception worth revisiting
because it ships per-sample WIG as well as CEL. GSE14864 (Wu and Nurse 2009, Cell, PMID 19269364)
provides per-chromosome normalised text alongside CEL, also on GPL7715. GSE88714 (Perrot 2018,
PLoS Genet, PMID 29466359) and GSE210821 (Lanteri 2022, Genes, PMID 35886027) use the same
coordinate-bearing Agilent platform as GSE98462 and are the natural companions to it.

## 3. Ranking, with the reasoning

The ranking below is by what a dataset adds to a hub that currently serves nothing on ASM294v2, divided
by what ingestion costs. Cost is stated as the concrete operation required, not as a score.

1. **GSE62108, Pu-seq (Daigaku 2015).** Adds the axis pombe is known for and the axis OriDB's origin
   records most need: which strand each polymerase copies, and hence fork direction, origin efficiency
   and termination position, genome-wide. No other candidate delivers four distinct derived profiles
   from one ingest. Cost: wig and bedGraph to bigWig conversion, plus establishing the chromosome naming
   used inside the files and mapping it to I/II/III. Because another group has already published a
   liftover of this dataset's origin list, the naming question must be answered before conversion, not
   after.
2. **GSE70033, BrdU-IP-seq replication rate (Daigaku 2017).** Adds the timing axis - the measurement a
   browser user expects to see under an origin - and does so with wild-type data in a form that needs no
   model refitting, because the submitter has already converted BrdU-IP signal into local replication
   rate at three fixed points in S phase. Cost: wig to bigWig, same 300 bp binning and same laboratory
   as the Pu-seq above, so no reconciliation work. The declared build is ASM294v2.23 rather than plain
   ASM294v2; point releases of ASM294v2 do not change chromosome lengths, but that should be confirmed
   against PomBase before the two are shown on one axis.
3. **GSE125335, Orc4 and Mcm2 ChIP-seq (Masuda 2020).** Adds licensing: where ORC sits and where MCM is
   loaded, which is the physical basis for the origin coordinates OriDB already serves and is not
   recoverable from any timing or polymerase-usage track. Cost: concatenating three per-chromosome wig
   files per sample and renaming `chr1/chr2/chr3` to `I/II/III`; low, but not zero.
4. **GSE65293, Rif1 and BrdU-IP (Kanoh 2015).** Adds the regulator layer - Rif1 occupancy next to the
   firing profile it suppresses - with nine ready-named bedGraph files on ASM294v2. Ranked below the
   first three only because its BrdU-IP profiles are HU-arrested, so they report early-origin firing
   rather than the unperturbed programme.
5. **GSE278836, GLOE-seq (Cheng, released 2025).** Adds nick-level strand asymmetry and is the only
   candidate that is already bigWig on ASM294v2, so ingestion cost is effectively zero. Ranked fifth
   because no publication is linked yet: the hub would be serving interpretation-free tracks.
6. **GSE125855, HydEn-seq (Zhou 2019).** Adds an independent chemistry for the same quantity as Pu-seq,
   which is exactly the comparison a browser makes easy and a figure does not. Cost: selecting the 24
   pombe samples out of 58 and checking that "assembly 2.23" means ASM294v2.23; the files are bigWig
   already.
7. **GSE93178, Mcm4 and BrdU (Matsumoto 2017).** Adds the Hsk1/Mrc1 axis on the same build and pipeline
   as GSE65293. Cost: unpacking a RAW tar, then bedGraph to bigWig.
8. **GSE98462, array origin efficiency (Gómez-Escoda 2018).** Adds a second laboratory's origin-efficiency
   measurement and a perturbation panel, and is the only array-era candidate whose platform record
   carries explicit `chrI:162-221` coordinates. Cost: a probe-to-coordinate join across 72 text files,
   then binning to bedGraph; moderate, and the resulting resolution is array-limited.
9. **GSE94812, BrdU-DNA-seq in telomere-protein mutants (Hasegawa, unpublished).** Adds subtelomeric
   timing, already in wig with PomBase chromosome names. Cost: near zero. Ranked here only because there
   is no linked publication to cite.
10. **GSE150800, MCM2-ChEC (Foss 2021).** Adds an orthogonal licensing chemistry, but only four pombe
    samples and the pombe arm is delivered as read-depth text rather than bedGraph. Cost: parsing the
    submitter's text format and binning it.
11. **GSE153731, Pu-seq at RTS1 (Naiman 2021).** Adds fork pausing and termination, a category otherwise
    unrepresented. Cost: the reference carries engineered barrier inserts, so either the insert interval
    is dropped or the track is served with a prominent caveat; plus CSV to bedGraph conversion of 72
    files.
12. **GSE241062, SCAR-seq (Charlton 2024).** Adds RFD from a chromatin method and 58 pombe samples of
    perturbation data. Cost is the highest of any sequencing-era candidate: the stated assembly is an
    in-house h90 contig set, so the bigWigs would have to be validated against ASM294v2 or the reads
    realigned.
13. **GSE6523, pre-RC tiling arrays (Hayashi 2007).** Adds the historical pre-RC map underlying much of
    OriDB's existing pombe content. Cost: CEL normalisation plus a probe join; the probe coordinates are
    on RefSeq chromosome accessions, which is the cleanest of the array-era coordinate situations.
14. **GSE4099, ssDNA under HU (Feng 2006).** Adds the only ssDNA/replication-stress measurement found
    for pombe. Cost: the values are not in supplementary files at all but in the GEO sample data tables,
    and the probe coordinates are in kilobases on a 2005 sequence.
15. **GSE6977, array timing (Eshaghi 2007).** Adds the classic timing profile, but its probes are
    indexed by ORF with no coordinates anywhere in the platform record, and only raw `.gpr` scans are
    deposited. Reconstructing a track means re-deriving probe positions from PomBase and refitting the
    published timing model - the highest cost on this list for a quantity that items 2 and 9 already
    supply at higher resolution.

## 4. Recommendation: ingest these three first

**GSE62108 (Pu-seq), GSE70033 (BrdU-IP-seq replication rate) and GSE125335 (Orc4/Mcm2 ChIP-seq).**

The reasoning is coverage, not just quality. These three are orthogonal: one measures which polymerase
copies which strand (fork direction, termination, origin efficiency), one measures when each locus is
copied, and one measures where the licensing machinery sits before any of it happens. Together they put
a quantitative, three-layer context under every one of the 741 pombe origin records the site already
serves, which is the specific gap this exercise exists to close. All three state ASM294v2 or a point
release of it in the archive record, so none requires a cross-assembly lift. All three deposit
coordinate-bearing processed files (wig and bedGraph), so none requires reprocessing raw reads with the
authors' software. They come from two laboratories and three assay families, so the first pombe section
of the hub will not read as one group's view. The combined conversion work is wig/bedGraph to bigWig,
one chromosome-name mapping to establish (GSE62108), and one per-chromosome file concatenation and
rename (GSE125335).

If a fourth is wanted at the same time, take **GSE278836** - its eight bigWigs are already on ASM294v2
and need no conversion at all - on the understanding that it is served as data without a citable
interpretation until the Martienssen paper appears.

## 5. What I verified and what I did not

Verified by live retrieval in this task: every accession above resolves in GEO; the title, sample count,
platform, instrument and supplementary-file inventory of each series as reported by the GEO esummary and
SOFT records; the `Genome_build` or `Assembly` statement in a representative sample of each series; the
presence or absence of coordinate columns in platform records GPL16235, GPL3847, GPL1932 and GPL3398;
and the PubMed metadata (first author, year, journal, PMID, DOI) for every linked publication.

Not verified: the contents of any data file, since nothing was downloaded - so chromosome naming inside
wig and bigWig files, bin sizes beyond what the submitter states, and the actual coordinate span of each
track remain unconfirmed. GPL7715 (Affymetrix Sp20b_M) returned no annotation table from GEO, so my
statement that the Masai-laboratory tiling arrays carry no coordinates in the archive is based on that
absence, not on a positive check elsewhere. The archive accessions for Heichinger 2006 and Xu 2012 could
not be established, as described in section 1. The claim that no pombe single-molecule replication
dataset exists is based on a title scan of all 143 public Oxford Nanopore SRA runs for the organism plus
GEO and PubMed searches; a dataset deposited under a non-pombe organism label, or in a repository
outside GEO/SRA/ENA/ArrayExpress, would not have been found.
