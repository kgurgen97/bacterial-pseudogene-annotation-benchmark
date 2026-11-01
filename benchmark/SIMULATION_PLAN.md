# SIMULATION_PLAN v1 — frozen 2026-10-30, before any genome is mutated

**This file is frozen.** Its SHA-256 is asserted at runtime by `benchmark/simulate.py` and
`benchmark/evaluate.py`, which refuse to run if it has changed. Nothing in it may be edited because
of a result. If a decision recorded here turns out to be wrong, the remedy is a new, separately
labelled analysis with its own plan file — never an edit to this one. This is the same discipline the
parent study's seven frozen plans follow.

Written after `SIMULATION_NOVELTY_AUDIT.md` (commit `f16d535`) and before any simulator code exists.

## 0. The question, and the one that earns the analysis its place

Four callers disagree about which bacterial genes are pseudogenes. The parent study established that
they disagree, that the disagreement survives phylogeny, and that selection corroborates
Pseudofinder's unique calls (+0.035 per genome, p = 1 × 10⁻¹¹) but not DFAST's larger set
(+0.006, p = 0.06, and 15/26 genomes in the reduced strata). It cannot say **why**, because natural
pseudogenes carry no label.

**Primary scientific question (P).** Do the four callers detect *different classes of
gene-inactivating lesion*, and can that difference account for the DFAST/Pseudofinder divergence?

**Secondary questions (S1–S4).** Per-caller sensitivity by lesion type and position; false-positive
rate on genes left untouched; whether either generalises across genera and reduction strata; and
computational cost.

## 1. Experimental unit, and why

**The experimental unit is the parental genome.** Genes within one genome share an assembly, a GC
content, a database representation and a translation table, so they are not independent replicates.
The parent study was forced to retrofit this (`CLUSTERING_RESULTS.md`: a pooled test at p = 3 × 10⁻⁵⁴
became p = 0.015 under a genome-clustered sign test). This plan adopts it from the start: every
headline number is computed per genome first, then aggregated across genomes.

Genus is the second clustering level, used where ≥ 2 genera contribute.

## 2. Genome selection, fixed now

From RefSeq, `assembly_level = "Complete Genome"`, with GenBank/RefSeq CDS annotation present.
Inclusion, all required:

1. **Single replicon** preferred; multi-replicon genomes allowed only if every replicon is retained
   and contig headers are preserved verbatim by every caller (§ 5). Plasmid-only contigs are kept in
   the input but excluded from the evaluation universe.
2. `transl_table = 11` for every CDS. Genomes using table 4 or 25 are **excluded** at this stage
   rather than special-cased; the exclusion is reported.
3. ≥ 1 500 annotated CDS, so the eligible set is large enough to lesion a small fraction.
4. Stratum assignment available from the parent study's `results/corpus_traits.tsv` or the expansion
   manifests, so results map onto § 3.5 of the manuscript. Target composition across the four strata
   S1 extreme obligate, S2 obligate intracellular, S3 recent restriction, S4 free-living breadth,
   with ≥ 3 distinct genera per stratum at Stage C.
5. No genome that the parent study's dN/dS or locus analyses used as a *partner* reference, to keep
   the benchmark from sharing an orthologue set with the analysis it is meant to inform.

Selection is by a fixed seeded draw within strata, recorded in `manifests/benchmark_genomes.tsv`
with accession, organism, genus, stratum, replicon count, CDS count and assembly accession.
**Genomes are never reselected after seeing a caller's output.**

## 3. The high-confidence intact gene set

Per parental genome, a CDS enters the eligible set only if **all** hold:

1. No `/pseudo` or `/pseudogene` qualifier in the source annotation.
2. Length is a multiple of 3, and ≥ 90 codons — long enough to carry distinct early, middle and late
   lesion positions.
3. Begins with a table-11 start codon (`ATG`, `GTG`, `TTG`) and ends with a stop codon.
4. No in-frame stop codon before the terminal one.
5. No `/exception`, `/transl_except`, `/ribosomal_slippage` or `/artificial_location` qualifier.
6. **Does not overlap any other annotated feature**, of any type, and has ≥ 30 bp of clear sequence
   on both flanks. Lesioning a gene that overlaps a neighbour makes truth unassignable for the
   neighbour.
7. ≥ 500 bp from any replicon end.
8. **Not a near-duplicate of another CDS in the same genome**: no other protein in the genome aligns
   at ≥ 90 % identity over ≥ 90 % of the shorter length (DIAMOND, self-search). A caller can rescue a
   broken gene by homology to its intact paralogue, which makes a miss uninterpretable.

Criterion 8 is the one most likely to be forgotten and most likely to be raised by a reviewer.

Every exclusion is counted by reason in `results/eligibility.tsv`. If fewer than 400 genes survive in
a genome, the genome is dropped and the drop reported.

## 4. Lesion menu, and the truth tiers

One lesion per gene. Lesioned genes are disjoint. Each parental genome yields **R replicate modified
assemblies**, each with its own seed; lesions are re-drawn per replicate, so a gene may be lesioned
in one replicate and untouched in another, which supplies matched negatives.

**Lesion density is fixed at 2.5 % of eligible genes per replicate, capped at 120 genes.** This is
deliberately below the pseudogene fraction of real reduced genomes, so global statistics such as
coding density move little and the callers are not pushed into a regime where they are responding to
the genome rather than to the gene. It is also the main difference in kind from Jundzill's
genome-wide 0.5–2 % nucleotide deletion.

### Disabling arms

| code | lesion | positions |
|---|---|---|
| `M1_stop_sub` | one substitution creating an in-frame `TAA`/`TAG`/`TGA` | early / middle / late |
| `M2_fs_ins1` | one nucleotide inserted | early / middle / late |
| `M3_fs_del1` | one nucleotide deleted | early / middle / late |
| `M4_fs_indel2` | two nucleotides inserted or deleted | early / middle / late |
| `M5_multi` | two or three independent disabling lesions in one gene, types drawn from M1–M4 | mixed |
| `M6_start_loss` | start codon mutated to a non-start codon | start |
| `M7_trunc_term` | terminal truncation: in-frame stop placed in the last 10 % of codons | late |
| `M8_disrupt_int` | internal segment of 30–150 bp deleted, length not a multiple of 3 | middle |

Positions are defined on codon index: **early** 5–20 %, **middle** 40–60 %, **late** 80–95 % of the
CDS. `M7` is distinct from `M1 late` in that it targets the final decile and is scored as a
truncation, not an internal stop — the task requires terminal truncation and internal disruption to
be separate categories, and they are.

### Control arms

| code | control | assumption |
|---|---|---|
| `C0_parental` | unmodified parental assembly | the baseline every caller is paired against |
| `C1_synonymous` | substitution verified to leave the translated protein byte-identical | negative |
| `C2_intergenic` | edit in an intergenic region overlapping no annotated feature, ≥ 30 bp from any | negative |
| `C3_inframe_indel` | 3 bp inserted or deleted | **ambiguous, not negative** |
| `C4_tech_error` | genome-wide random single-nucleotide deletions at 0.05 %, no targeted gene lesions | technical, scored separately |

`C3` is labelled ambiguous because an in-frame indel may or may not disable a protein; it is never
counted as either a true positive or a false positive in the primary analysis. `C4` is the
sequencing/assembly-artefact scenario that Cooley & Wright make mandatory; it is reported in its own
section and **never pooled** with the biological arms.

### Truth tiers, assigned at generation time

- **`high`** — `M1` early/middle, `M2`, `M3`, `M4` early/middle, `M5`, `M6` *only where no alternative
  table-11 start codon exists in frame within the first 10 % of codons*, verified per gene.
- **`medium`** — `M1`/`M4` late, `M7`, `M8`.
- **`ambiguous`** — `C3`; `M6` where a downstream alternative start exists.
- **`negative`** — `C0`, `C1`, `C2`, and every eligible gene in a replicate that received no lesion.
- **`technical`** — `C4`.

**The primary analysis uses `high` only.** `medium` is a pre-registered secondary. `ambiguous` is
reported descriptively and never enters a confusion matrix. Not every mutated gene is a true
pseudogene and this plan does not pretend otherwise.

### Per-lesion record

For every lesion: parental accession, replicate id, seed, gene locus tag, protein id, contig,
parental start/end/strand, codon index, lesion code, position stratum, original and modified local
sequence (21 bp window), the full coordinate delta, expected effect, truth tier, and for `M6` whether
an alternative start was found. Written to `results/ground_truth.tsv` plus a Parquet copy.

**Coordinate liftover.** Every insertion and deletion shifts downstream coordinates. A per-replicate
liftover map (parental → modified, and the inverse) is emitted as `results/liftover/<replicate>.tsv`
and is the only sanctioned way to translate coordinates. No analysis may assume parental coordinates
hold in a modified assembly.

### Generation-time validation, blocking

`simulate.py` exits non-zero without writing results unless, for every lesion:

1. The modified FASTA at the recorded modified coordinate contains exactly the intended bases.
2. Re-translating the modified gene in its original frame produces the expected consequence — a
   premature stop for `M1`/`M7`, a frame change for `M2`/`M3`/`M4`/`M8`, a lost start for `M6`.
3. `C1` substitutions leave the translated protein identical.
4. `C2` edits fall in no annotated feature.
5. The liftover map round-trips: parental → modified → parental is the identity on all unlesioned
   positions.
6. No two lesions in one replicate fall in the same gene or within 30 bp of each other.
7. Contig count, contig order and contig names are unchanged from the parental FASTA.
8. The parental FASTA is byte-identical to the downloaded original.

## 5. Running the callers

All four are specified. Execution is gated on availability, not on convenience.

| caller | version | invocation requirements |
|---|---|---|
| Bakta | 1.12.2 | `--keep-contig-headers`, `--translation-table 11`, fixed database version recorded |
| DFAST | 1.5.3 | `--use_original_name t`, `--gcode 11` |
| Pseudofinder | clone `ac5051a` | `annotate` with a fixed protein database; the parent study's two configurations are not mixed |
| PGAP | `2026-06-18.build8602` | **requires a container runtime; none is installed on this machine** |

`--keep-contig-headers` and `--use_original_name` are not optional conveniences. They make each
caller's output contig names identical to the input, which — combined with knowing the exact
coordinate of every lesion created — means truth is matched **inside each caller's own output**, with
no cross-caller contig reconciliation anywhere. The parent study resolved 26 % of 886 176 interval
lookups by file order (`LOCUS_RESULTS.md`); that uncertainty is designed out here rather than
inherited, and any run in which a caller's output contig names do not match the input is a failed run,
not a run to be reconciled.

**Every caller is run on both the parental and the modified assembly, by us.** Comparing a mutated
assembly against the unchanged RefSeq annotation of its parent is not a benchmark of PGAP or of
anything else, and will not be done under any label. This is why PGAP cannot be substituted for by
RefSeq: it must be rerun or reported as not run.

**PGAP contingency.** PGAP's own documentation offers no native build and requires Docker, Podman or
Singularity/Apptainer; its hardware needs ("8 cores, 16G RAM, and 100GB of storage", ~25 min for a
5.7 Mb genome) are well inside this machine's 20 cores / 23 GB / 772 GB. The blocker is therefore one
system package plus a ~30 GB supplemental-data download — a change outside this repository, not to be
made without the owner's explicit go-ahead. Until it is made, the benchmark **executes three callers
and delivers the PGAP arm as a reproducible execution plan** in `SIMULATION_METHODS.md`, with its
absence stated in the results and in any manuscript text. A three-caller result is reported as a
three-caller result.

Versions, databases, parameters, exact command lines, wall clock, peak RSS and CPU time are recorded
per run in `results/tool_runs.tsv`, and the environment is dumped to `logs/benchmark_env.log`.

## 6. Output normalisation — classes kept apart

A documented parser maps each caller's output for a gene locus to exactly one of:

`explicit_pseudogene` · `disrupted` · `partial` · `intact_cds` · `absent` · `ambiguous`

Rules fixed now:

- `explicit_pseudogene` requires an explicit pseudogene assertion: a `/pseudo` or `/pseudogene`
  qualifier, feature type `pseudogene`, or Pseudofinder's pseudogene output.
- **DFAST `partial` is `partial`, not `explicit_pseudogene`.** It is scored under endpoint B and never
  under endpoint A.
- **`absent` is `absent`, not a pseudogene call.** A missing CDS is scored under endpoints B and C and
  never counted as an explicit pseudogene call.
- `disrupted` covers a gene split into ≥ 2 CDS features, or a feature annotated as frameshifted or
  internally stopped without an explicit pseudogene assertion.
- Where a caller's output supports two classes for one locus, the class is `ambiguous` and the gene is
  excluded from that caller's primary analysis, counted and reported.

The parser is driven by a fixture file of hand-written caller output covering every class, and its
tests are part of the Stage A gate.

## 7. The three endpoints

Scored per truth gene, so one disrupted gene split into several annotations can never count as
several true positives. The evaluation universe is the set of truth genes, not the set of predictions.

- **Endpoint A — explicit pseudogene classification.** Positive iff the overlapping call is
  `explicit_pseudogene`.
- **Endpoint B — disruption or abnormal structure detected.** Positive iff
  `explicit_pseudogene`, `disrupted`, `partial`, or `absent`.
- **Endpoint C — gene recovered as any annotated feature.** Positive iff any feature overlaps the
  lesioned gene's lifted interval.

**All endpoints are paired against `C0`.** A gene counts as a true positive only if the caller calls
it positive in the modified replicate and **did not** call it positive in the parental genome. A gene
counts as a false positive only if the caller calls it positive in the modified replicate, it carried
no lesion, and it was not already positive in the parental. This removes each caller's baseline call
rate from every number and is the only way the comparison is fair between a conservative caller and a
liberal one. Unpaired figures are reported as a secondary.

Overlap for matching is ≥ 50 % reciprocal with the lesioned gene's lifted interval, the same
criterion as `LOCUS_PLAN.md`, computed in the modified assembly's own coordinates.

## 8. Metrics

Per caller, per arm, per genome: TP, FP, FN, TN; sensitivity; specificity; precision; F1; balanced
accuracy; and false-positive rate on untouched eligible genes. Precision and F1 are reported **only**
where the negative universe is the full untouched eligible set, since both depend on prevalence and
prevalence here is set by design at 2.5 %; the plan states now that **sensitivity and FPR are the
interpretable pair** and precision is contextual.

Breakdowns, pre-registered: by lesion code; by position stratum; by truth tier; by genus; by
reduction stratum; and computational cost per genome.

**Uncertainty.** Percentile bootstrap over genomes, 2 000 resamples, seed 20261030, resampling whole
genomes with all their genes. Where ≥ 2 genera contribute, a genus-clustered bootstrap is reported
alongside. No confidence interval is computed from pooled genes.

**Caller comparison.** Matched, because every caller sees the same genes and the same lesions:
McNemar's exact test on paired per-gene calls as the pooled view, and — as the primary — a sign test
on per-genome sensitivity differences, requiring ≥ 5 non-tied genomes, with exact ties dropped and the
number dropped always reported. This mirrors the parent study's corrected practice.

**Multiple testing.** Benjamini–Hochberg within the primary family, which is the caller-pair
comparisons under endpoint A on `high`-tier lesions. The per-lesion-type and per-stratum breakdowns
are a separate exploratory family, BH-corrected within themselves and labelled exploratory.

**Failed runs.** If any caller fails on any replicate of a genome, that genome is dropped for **all**
callers — complete cases only — and every failure is listed in `SIMULATION_RESULTS.md` with its error.
Dropping a genome for one caller only would bias the comparison toward whichever caller is most
fragile. A caller that produces output but no call for a gene yields `absent`, which is a result, not
a failure.

## 9. Pre-registered predictions, with the bars that decide them

Each is **SUPPORTED**, **UNSUPPORTED**, or **NOT TESTABLE**. Bars are fixed here.

| # | prediction | bar |
|---|---|---|
| **B1** | **gate.** Every executed caller detects ≥ 50 % of `high`-tier `M1_stop_sub` lesions under endpoint B | if no caller clears it, the pipeline or parser is broken, not the callers; the run is void and nothing downstream is reported |
| **B2** | Pseudofinder's endpoint-A sensitivity exceeds DFAST's on frameshift lesions (`M2`+`M3`+`M4`) | sign test over genomes, p < 0.05 after BH |
| **B3** | **the novelty claim.** Within at least one caller, endpoint-A sensitivity differs across lesion codes | the widest per-caller gap between two lesion codes exceeds 20 percentage points with non-overlapping bootstrap intervals |
| **B4** | Every caller's FPR on untouched eligible genes, paired against `C0`, is < 2 % | upper bootstrap bound < 2 % |
| **B5** | **the explanation.** DFAST is more sensitive than Pseudofinder on at least one lesion code | if DFAST is more sensitive on none, the "nested plus noise" reading of `DFAST_RESULTS.md` is favoured over "different lesion classes", and that is reported as the answer even though it is the less interesting one |

B1 gates everything. B5 is written so that both outcomes are publishable and neither can be reported
as a vindication of the parent study's preferred reading.

**`C4` is not used to score any caller's accuracy.** It is reported as a separate observation:
how many unexpected pseudogene calls each caller makes per megabase under technical error alone.

## 10. Stages and gates

No stage begins until the previous one's gate is met, in writing.

**Stage A — validation. 3 genomes, deliberately diverse, small known lesion sets.**
Gate, all required: every one of the eight generation-time validations in § 4 passes; FASTA validity
confirmed; liftover round-trips; contig names identical in every caller's output; the parser
reproduces a hand-written fixture for all six classes; every confusion-matrix cell verified against a
hand-computed toy example; unit and integration tests pass. No performance number is reported from
Stage A.

**Stage B — pilot. 10–20 genomes, ≥ 4 genera, all arms.**
Gate: caller failure rate < 10 % of runs; measured runtime and storage extrapolate within budget to
Stage C; B1 met; and at least one arm shows a detection rate strictly between 0 and 100 % for some
caller, since an all-or-nothing pilot cannot discriminate between callers at any scale. Every failure
is documented. **Neither genomes nor lesions are reselected on the basis of which caller does well.**

**Stage C — scaled benchmark. 50–100 genomes**, final n and the replicate count R chosen from the
pilot's observed between-genome variance and the compute budget, not from a significance target. The
choice and its arithmetic are written into `SIMULATION_METHODS.md` before Stage C runs.

Execution is resumable and checkpointed per (genome, replicate, caller). Disk, memory, process count
and output size are monitored; the budget is a hard cap, and exceeding it stops the run rather than
silently dropping a caller or a condition.

## 11. What the benchmark may and may not conclude

- Synthetic lesion detection is a property of **software**. It is evidence about what a caller can
  see, not about whether a caller's calls on natural genomes are biologically correct. The two are
  reported in separate sections and the synthetic result is never described as biological accuracy.
- A caller with low sensitivity here is not thereby shown to be wrong about natural pseudogenes, and a
  caller with high sensitivity is not thereby validated. Natural pseudogenes are old, multi-cause and
  unlabelled; these lesions are recent, single-cause and known.
- Consistency with the dN/dS results is assessed as **concordance or discordance between two
  independent lines of evidence**, not as one validating the other.
- Negative, inconclusive and failed results are reported in `SIMULATION_RESULTS.md` with the same
  prominence as positive ones.

## 12. Deliverables

`SIMULATION_NOVELTY_AUDIT.md` (done, `f16d535`) · this plan · `SIMULATION_METHODS.md` ·
`SIMULATION_RESULTS.md` · `SIMULATION_LIMITATIONS.md` · `manifests/benchmark_genomes.tsv` ·
`results/ground_truth.tsv` and `.parquet` · `results/liftover/` · `results/eligibility.tsv` ·
`results/tool_runs.tsv` · `results/predictions.tsv` · `results/confusion.tsv` ·
`results/summary_genome.tsv`, `results/summary_genus.tsv` · `figures/` · `workflow/` ·
`tests/` · `logs/` · and a scientific review answering whether this strengthens the manuscript.

Raw assemblies and caller output stay out of Git; manifests, checksums and retrieval instructions are
committed.
