# Does the pseudogene caller change what the literature concludes?

**Written 2026-10-21, before any four-caller association was computed.** Bakta was 22/363 and DFAST
342/363 at the time of writing, so no result existed for either. Nothing here may be changed because of
an outcome.

## 1. Why the design had to change

The first version of this project asked whether swapping the pseudogene caller changes the verdicts of
three of my own unpublished studies (H10, H11, H12). It does not — but that is weak evidence, because
those three are not published, two of them are nulls, and a null is hard to break. A reviewer would ask,
correctly, why those three and why they should stand for the field.

So I went looking for published cross-genome associations with pseudogene fraction to re-test instead.
**There are almost none.** The literature on bacterial pseudogene fraction consists of:

- methods papers (Lerat & Ochman's Ψ-Φ 2004; Pseudofinder 2022, which benchmarks on **two** genomes);
- single-lineage case studies (*Shigella*, *Yersinia pestis*, *Serratia symbiotica*, *M. leprae*);
- reviews asserting qualitatively that host-restricted and intracellular bacteria carry elevated
  pseudogene fractions;
- a 2024 paper arguing many PGAP-flagged pseudogenes are not pseudogenes at all, and tracing the
  inflation partly to assembler choice.

That absence is itself worth reporting. It also means "re-test N published regressions" is not available.
What *is* available, and is stronger, is to take the field's **canonical claims** — the ones every review
makes and on which the interpretation of pseudogene counts rests — and ask whether they are detected the
same way by four different callers on one corpus.

## 2. Corpus (frozen) and its limits, stated up front

363 complete bacterial genomes, already annotated by four callers:
**PGAP** (RefSeq), **Pseudofinder**, **Bakta** 1.12.2 (DB 6.0 full), **DFAST** 1.5.3.

Counting rules were fixed earlier, per caller, and are not revisited here: PGAP
`pseudo/(pseudo+CDS)`; Pseudofinder `pseudos/(pseudos+intact)`; Bakta the `pseudogenes:` line of its own
summary; DFAST the disrupted classes only (frameshift + internal_stop_codon + both), never `partial`.

**The corpus is two phyla, not bacteria in general**: Pseudomonadota 259 (of which 137 are
Gammaproteobacteria from the H11 set) and Actinomycetota 104. Genome size spans 0.88–12.01 Mb, which is
the range the reduction claim needs. Every analysis below is run **within phylum** as well as pooled, and
the paper will not claim generality beyond these two phyla.

## 3. The claims under test, as the literature makes them

**C1 — genome reduction.** Pseudogene accumulation accompanies genome reduction, so pseudogene fraction
should be higher in smaller genomes. Predicted direction: **negative** correlation between pseudogene
fraction and genome size. This is the mechanistic core of every review's interpretation.
*Label-free and available for all 363 genomes.*

**C2 — host restriction.** Host-associated and intracellular bacteria carry elevated pseudogene fractions
relative to free-living ones. Predicted direction: **host > free-living**.
*Requires a lifestyle label; see section 4.*

C1 is the primary because it needs no external label and uses the whole corpus. C2 is secondary because
its labels are the weak link, and that weakness is a property of available metadata, not of the callers.

## 4. Lifestyle labels for C2 — and why the primary set is small

Labels come from Madin et al. 2020 (*Scientific Data*, CC-BY), field `isolation_source`, joined to the
corpus in a fixed cascade: NCBI `tax_id`, then species name, then genus.

- `host*` → host-associated; `soil* / water* / sediment* / rock*` → free-living; everything else excluded.
- **Primary C2 set: tax_id and species joins only — 74 genomes (32 host, 42 free-living)**, balanced
  across both phyla (Actinomycetota 19/20, Pseudomonadota 13/22).
- A genus-modal label is **not** a species' lifestyle: one genus can hold free-living saprophytes and
  obligate pathogens. The genus-fallback set (211 genomes, 96 host / 115 free-living) is therefore
  reported as **exploratory only** and never as the primary.
- `doubling_h`, which would have given an effective-population-size proxy, covers only 43/363 and is
  **not used**.

## 5. Statistics

- **C1**: Spearman correlation of pseudogene fraction against genome size, and partial Spearman given
  GC content, pooled and within each phylum.
- **C2**: Mann–Whitney U on pseudogene fraction, host vs free-living, pooled and within phylum.
- Every test is run **independently for each of the four callers** on the identical genome set.
- α = 0.05, two-sided throughout. No multiplicity correction across callers, because the callers are not
  competing hypotheses — each is a separate analysis of the same claim, and correcting across them would
  obscure precisely the disagreement being measured. Stated here so it is not mistaken for an oversight.

**A mechanical dependency that must be reported, not hidden.** Pseudogene fraction has CDS count in its
denominator, and CDS count scales with genome size, so C1 is not free of arithmetic coupling. All three
quantities are therefore reported for every caller: the fraction, the raw pseudogene count, and
pseudogenes per Mb. If the callers disagree about C1, the report states which of the three measures the
disagreement appears in.

## 6. Primary endpoint — declared now

> For each claim, two numbers:
> **(a) how many of the four callers return a significant association in the predicted direction**, and
> **(b) the ratio of the largest to the smallest absolute effect size across the four callers.**

Pre-committed readings:

| outcome | what it means |
|---|---|
| 4 of 4 agree in direction and significance | the claim is robust to caller choice; (b) then quantifies how much the *magnitude* a paper would report depends on its caller |
| callers split on significance | **the headline**: whether a published claim is detectable depends on which caller the authors happened to use |
| callers split on direction | severe; the measure cannot support cross-study comparison at all |
| no caller finds the association | the canonical claim does not hold in this corpus, which is a finding about these two phyla and not about the callers |

The last row matters: this design can conclude that the field's canonical claim fails here, and that
outcome gets written up as readily as the others.

## 7. What this cannot do

- It does not re-test any specific published paper's analysis, because the papers do not report
  cross-genome regressions to re-test. It tests the claims those papers assert.
- Two phyla only.
- One configuration per caller. Pseudofinder with a different reference, or Bakta with the light
  database, would give different absolute counts.
- C2's labels are species- or genus-level ecology, not measured lifestyle for these strains.

## 8. Deliverables

`analysis_literature_retest.py`, `results/literature_retest.{tsv,json}`, and a results section that
reports the primary endpoint first. The script is to run on whatever callers are complete and be rerun
when Bakta finishes; the staged run and the final run are both recorded.
