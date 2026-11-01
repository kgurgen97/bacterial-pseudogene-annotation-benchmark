# Do the configurations differ in effect magnitude, or only in verdict?

Executes `CALLER_EFFECT_CONTRAST_PLAN.md` (sha256 `0d1d6458…`), frozen and committed at `f0aacc7`
**before** the contrasts were computed. Results: `results/caller_effect_contrast.json`,
`results/caller_effect_contrast/*.tsv`. No annotation re-run, no genome downloaded, no dN/dS
recomputed.

---

## 1 The error being corrected

`MANUSCRIPT_Q1_FINAL_CANDIDATE.md § 9.3` concluded that annotation choice changes a downstream result
because PGAP (p = 0.016) and Pseudofinder (p = 7.0 × 10⁻⁶) reported a significant genome-reduction
association while Bakta (p = 0.065) and DFAST (p = 0.228) did not.

**Comparing verdicts is not comparing effects.** Two correlations on opposite sides of p = 0.05 can be
statistically indistinguishable from each other. The claim requires the **paired difference** of the
correlations, with an interval.

## 2 Design

363 genomes, **362 eligible**, **185 genera**. For each metric and each pair of configurations,
Δρ = ρ_A − ρ_B on the **same genomes**, with a **genus cluster bootstrap** (B = 5 000, seed 20261101)
that recomputes both correlations inside every replicate so the contrast is paired. Holm adjustment
across the six pairs within each metric and block. Two-sided bootstrap p floored at 1/B.

**Pseudofinder is quarantined.** It post-processes PGAP output, so pairs involving it are not two
independent measurements. The headline reading uses the **independent block only**: PGAP–Bakta,
PGAP–DFAST, Bakta–DFAST.

> **Deviation D1 from the frozen plan, § 2.** The plan said to parse genus from
> `results/corpus523.tsv`'s organism string. That column is populated for the 160 *expansion* genomes
> only; for these 363 it is empty, so the plan as written yields **zero** eligible genomes. Genus was
> taken instead from `results/corpus_traits.tsv`, which covers exactly these 363 and carries a parsed
> genus column, with the corpus523 parse kept as a fallback. **The frozen plan was not edited.** This
> is the same lineage-source mismatch `PHYLO_RESULTS.md` documented for its own plan.

## 3 Primary metric — pseudogene **fraction**

### 3.1 Verdict table — what the manuscript previously reported

| configuration | ρ | p | verdict | partial ρ given GC |
|---|---|---|---|---|
| PGAP | +0.1260 | 0.0163 | **significant** | +0.1118 |
| Pseudofinder | +0.2334 | 7.0 × 10⁻⁶ | **significant** | +0.1538 |
| Bakta | +0.0968 | 0.0654 | n.s. | +0.1258 |
| DFAST | +0.0634 | 0.228 | n.s. | −0.0488 |

### 3.2 Contrast table — what the claim requires

**Independent pairs (headline):**

| pair | Δρ | 95 % CI | p | Holm p | label |
|---|---|---|---|---|---|
| PGAP − Bakta | +0.0292 | [−0.0726, +0.1304] | 0.570 | 1.000 | **verdict differs, magnitudes not distinguishable** |
| PGAP − DFAST | +0.0626 | [−0.0129, +0.1344] | 0.104 | 0.312 | **verdict differs, magnitudes not distinguishable** |
| Bakta − DFAST | +0.0334 | [−0.0799, +0.1427] | 0.582 | 1.000 | no difference detected |

**Pairs involving Pseudofinder — not independent of PGAP, excluded from the headline:**

| pair | Δρ | 95 % CI | Holm p |
|---|---|---|---|
| PGAP − Pseudofinder | −0.1074 | [−0.1878, −0.0272] | 0.020 |
| Pseudofinder − Bakta | +0.1366 | [+0.0173, +0.2572] | 0.028 |
| Pseudofinder − DFAST | +0.1700 | [+0.0839, +0.2542] | 0.0006 |

> ### Declared reading: **RETREAT**
> On the primary metric, **no independent pair has an interval excluding zero**, adjusted or not.
> Per § 7 of the frozen plan, the manuscript must say that the **verdict** differs while the
> **magnitudes are not distinguishable**, and the title, abstract, § 9.3 and the Conclusions must be
> corrected. They have been.

## 4 Secondary metrics — and the finding that replaces the overclaim

### 4.1 Raw count: independent configurations **do** differ

| pair | Δρ | 95 % CI | Holm p | label |
|---|---|---|---|---|
| **PGAP − Bakta** | **+0.2084** | **[+0.1225, +0.2967]** | **0.0006** | **magnitudes differ** |
| **Bakta − DFAST** | **−0.1879** | **[−0.2792, −0.0975]** | **0.0006** | **magnitudes differ** |
| PGAP − DFAST | +0.0205 | [−0.0315, +0.0724] | 0.432 | no difference detected |

All four configurations are individually significant on this metric (ρ +0.29 to +0.88), so here the
verdicts agree while two of three **magnitudes genuinely differ**.

### 4.2 Per Mb: like fraction, nothing distinguishable

| pair | Δρ | 95 % CI | Holm p |
|---|---|---|---|
| PGAP − Bakta | +0.0376 | [−0.0629, +0.1380] | 0.920 |
| PGAP − DFAST | +0.0605 | [−0.0131, +0.1331] | 0.344 |
| Bakta − DFAST | +0.0229 | [−0.0898, +0.1337] | 0.920 |

### 4.3 The pattern across metrics is the real result

| metric | normalisation | independent pairs with distinguishable magnitudes |
|---|---|---|
| **count** | none | **2 of 3** |
| fraction | ÷ that configuration's own CDS count | 0 of 3 |
| per Mb | ÷ genome length | 0 of 3 |

**Whether the configurations appear to disagree about how pseudogene burden scales with genome size
depends on the metric.** On raw counts they demonstrably differ — unsurprisingly, since they emit
between 6 844 and 130 748 calls corpus-wide. Dividing by each configuration's own CDS count, or by
genome length, absorbs most of that, and the configurations then become statistically
indistinguishable in how burden *scales*.

This is weaker than "annotation choice changes the conclusion" and more useful: it says where the
pipeline-dependence lives, and where it does not.

## 5 Magnitude against verdict, the two confusions kept apart

| situation | pairs | reading |
|---|---|---|
| verdicts differ, magnitudes do **not** | PGAP–Bakta and PGAP–DFAST on fraction and per Mb | the familiar error; **cannot** support a claim that effects differ |
| verdicts agree, magnitudes **do** differ | PGAP–Bakta and Bakta–DFAST on count | a real difference that a verdict comparison would have **missed entirely** |
| neither | Bakta–DFAST on fraction and per Mb | no evidence of difference, which is not evidence of equivalence |

The second row is the point: comparing significance verdicts is not merely over-permissive, it is also
**under-powered in the other direction** — it misses genuine differences when both configurations
happen to clear the threshold.

## 6 Limitations

1. **No equivalence margin was pre-specified**, so a non-significant contrast is **absence of
   evidence**, not evidence of equivalence. With 362 genomes in 185 genera, power to resolve a
   Δρ of ~0.03 is low, and the fraction-metric nulls are consistent with real differences of that size.
2. **Genus clustering is coarse.** 185 genera over 362 genomes means most genera contribute one
   genome; the cluster bootstrap is close to a genome bootstrap here, and gains little.
3. **Spearman ρ differences are not a natural effect size for comparison**; they are bounded and
   non-additive. The contrast is interpretable in direction and rough magnitude, not as a ratio.
4. **The metric comparison in § 4.3 was not pre-specified as a hypothesis**; the three metrics were
   pre-specified, but the pattern across them is a post-hoc observation and is labelled as such.
5. **Pseudofinder's contrasts cannot be read as pipeline comparisons** at all.
6. This says nothing about **which** configuration is right.

## 7 What changed in the manuscript

| location | before | after |
|---|---|---|
| title | "determines whether a canonical genome-reduction result is significant" | names the metric dependence instead |
| abstract | "two configurations report a significant relationship and two report none" framed as the consequence | states that on fraction the magnitudes are **not** distinguishable, and that the difference is on raw counts |
| § 9.3 | verdict comparison as the practical result | both tables, the retreat, and the metric pattern |
| § 12.1 | "the reported result is pipeline-dependent" | narrowed to counts, with normalisation named as what removes it |
| Conclusions | item 3 | rewritten |
