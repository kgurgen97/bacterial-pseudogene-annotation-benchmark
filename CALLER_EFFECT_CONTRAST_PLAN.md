# Frozen plan: do the configurations differ in effect magnitude, or only in verdict?

**Frozen before the contrasts were computed.** `analysis_caller_effect_contrast.py` asserts this
file's SHA-256 at runtime. A frozen plan is never edited because of a result.

---

## 1 The gap this closes

`MANUSCRIPT_Q1_FINAL_CANDIDATE.md § 9.3` reports that the four configurations give a **3.68-fold
spread** in effect size on the canonical genome-reduction association, and that two report a
significant relationship while two do not (PGAP p = 0.016 and Pseudofinder p = 7.0 × 10⁻⁶ against
Bakta p = 0.065 and DFAST p = 0.228).

**A significant result under one configuration and a non-significant result under another is not
evidence that the two effects differ.** Two correlations straddling an arbitrary threshold can be
statistically indistinguishable from each other. The manuscript currently draws its practical
conclusion from exactly that comparison, which is a known statistical error and a fair target for
desk rejection.

This plan specifies the **direct paired contrast** that the claim actually requires.

## 2 Data

`results/literature_retest.tsv`, 363 genomes, already computed; no annotation is re-run and no genome
is downloaded. Genus is parsed from the organism string in `results/corpus523.tsv` by the same rule
`analysis_dnds.py` uses. A genome is eligible when all four configurations have a non-missing measure
and a genus can be resolved.

## 3 Estimands

For each metric *m* ∈ {**fraction** (primary), count, per Mb} and each ordered pair of configurations
(A, B):

**Δρ(A, B | m) = ρ_A − ρ_B**, where ρ is the Spearman correlation between the per-genome pseudogene
measure and genome size, **computed on the same genomes**.

Secondary: the same contrast on the **partial** Spearman given GC content.

The per-configuration ρ values themselves are already published (§ 9.3) and are reported again here
only as context.

## 4 Inference

- **Cluster bootstrap over genera**, because genomes within a genus are not independent. Genera are
  resampled with replacement; within a resample both ρ's are recomputed on the *same* genomes, so the
  contrast is paired at every replicate.
- **B = 5 000** replicates, seed **20261101**, one fresh generator per metric.
- Report Δρ, the 2.5th and 97.5th percentiles, and a two-sided bootstrap p-value
  `2 × min(P(Δ* ≤ 0), P(Δ* ≥ 0))`, floored at `1/B`.
- **Multiplicity**: 6 unordered pairs per metric. Holm–Bonferroni within metric. Both unadjusted and
  adjusted p-values are reported.
- Eligible genomes and genera are reported for every contrast.

## 5 Pseudofinder's dependence on PGAP

Pseudofinder post-processes PGAP-derived annotation, so the pair (PGAP, Pseudofinder) is **not two
independent measurements** and its contrast cannot be read as a comparison of independent pipelines.
Contrasts are therefore reported in two blocks:

- **Independent block**: PGAP–Bakta, PGAP–DFAST, Bakta–DFAST.
- **Dependent block**: every pair involving Pseudofinder, reported with the dependence stated inline
  and excluded from the headline reading.

The headline reading uses the **independent block only**.

## 6 Magnitude versus verdict — reported separately and never merged

Two tables:

1. **Verdict table**: per configuration, ρ, p and whether p < 0.05. This is what the manuscript
   currently reports.
2. **Contrast table**: Δρ with CI and bootstrap p. This is what the claim requires.

A pair is labelled **"verdict differs, magnitudes not distinguishable"** when the two configurations
fall on opposite sides of p = 0.05 while the Δρ interval contains zero. That label is declared now so
it cannot be avoided later.

## 7 Declared readings — fixed before the result is seen

| observation | what the manuscript must say |
|---|---|
| **≥ 1 independent pair** has a Δρ interval excluding zero after Holm adjustment | configurations differ in effect **magnitude**, for those pairs, and the claim stands in a form naming them |
| **no independent pair** survives adjustment, but ≥ 1 does unadjusted | the evidence is **suggestive only**; report unadjusted and adjusted side by side and do not claim a magnitude difference |
| **no independent pair** has an interval excluding zero even unadjusted | **retreat**: the manuscript must say the **verdict** differs while the **magnitudes are not distinguishable**, and the title, abstract, § 9.3 and the Conclusions must be corrected accordingly |

The third row is the outcome that costs the paper its strongest-sounding sentence. It is written down
now precisely so that it cannot be renegotiated after the fact.

## 8 What this cannot do

- It cannot establish which configuration is correct: there is no ground truth here.
- It cannot separate a genuine difference in detection from a difference in definition.
- A null contrast does not mean the configurations agree; with 363 genomes clustered in genera the
  power to resolve small Δρ is limited, and that limit is reported rather than interpreted as
  equivalence.
- No equivalence margin is pre-specified, so a non-significant contrast is **absence of evidence**.

## 9 Scope

This plan authorises the contrasts above on existing data. It authorises no new annotation, no new
dN/dS, no RBH re-run and no downloads.
