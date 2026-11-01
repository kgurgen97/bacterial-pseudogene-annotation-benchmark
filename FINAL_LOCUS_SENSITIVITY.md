# G3 after correction: what it survives, and the one place it weakens

2026-10-31. **Post hoc throughout.** Nothing here is registered and nothing here is a verdict. G3's
registered verdict comes from the frozen `LOCUS_PLAN.md` analysis; this asks how far the number moves
when defensible choices are made differently. `LOCUS_PLAN.md` is untouched.

## The corrected result being tested

| | |
|---|---|
| scored loci | 103 224 |
| discordant | 4 886 (was 7 565 before the classifier fix) |
| concordant-intact | 98 338 |
| **G3** | **+0.0307** per genome over **92 genomes**, 78 positive, sign test **p = 6.16 × 10⁻¹²** |

## The battery

> ### Corrected after adversarial review, 2026-10-31
>
> Three statements in the original version of this section were wrong, and all three are fixed below.
>
> 1. **"Fifteen checks including a genus-clustered bootstrap" was false.** The count of 15 comes from
>    a filter that keeps only entries carrying a top-level `p`, and S1 (the bootstraps) and S2 (the
>    leave-one-genus-out refits) do not. Neither was ever among the 15. Verified:
>    `results/review_response.json` `R4`.
> 2. **Two of the 15 are degenerate** — `S4_min_3` and `S5_without_sole_PGAP` reproduce the base
>    statistic digit for digit (the latter drops 30 loci). So the honest count is **13 non-degenerate
>    checks, plus the bootstraps and 29 refits reported separately.**
> 3. **The S1 interval was for the wrong estimator.** It is a CI for the *pooled* gap, not for the
>    reported median of per-genome gaps — `analysis_locus_sensitivity.py` prints that caveat itself,
>    but this section then used it as though it belonged to +0.0307. The correct interval is now
>    computed and is **tighter**, not looser (see S1 below).

Thirteen non-degenerate variants. **All stay positive and significant at 0.05**, and so do the
bootstraps and all 29 leave-one-genus-out refits.

### S1 — clustering respected

Two different estimators, kept apart. The **reported** statistic is the median of per-genome gaps;
the pooled gap is a different quantity.

| estimator | unit resampled | 95 % CI | |
|---|---|---|---|
| **median of per-genome gaps — the reported statistic** | **genus** | **[+0.0206, +0.0401]** | **excludes 0** |
| pooled gap | genome | [+0.0169, +0.0553] | excludes 0 |
| pooled gap | genus | [+0.0167, +0.0524] | excludes 0 |

Clustering by genus barely widens either interval, and the CI for the statistic actually quoted is
the tightest of the three. This is the check most likely to have killed G3 and it does not.

**But the headline p must change.** A genus-level sign test — one value per genus, which is the
honest unit given up to 6 genomes per genus and 6 genera with ≥ 5 — gives **+0.0302, 25 of 28 genera
positive, p = 2.7 × 10⁻⁵**. The per-genome p of 6.16 × 10⁻¹² overstates the evidence by roughly seven
orders of magnitude, which is the very error `DNDS_RESULTS.md` criticised § 3.7 for. **p ≈ 3 × 10⁻⁵ is
the defensible figure.**

### S2 — no single genus carries it

Every one of the 29 leave-one-genus-out refits stays positive and significant. The worst four lose 5 %
of the effect:

| dropped | genomes | gap | p |
|---|---|---|---|
| *Aliarcobacter* | 87 | +0.0292 | 8.6 × 10⁻¹¹ |
| *Bordetella* | 89 | +0.0292 | 3.0 × 10⁻¹¹ |
| *Chlamydia* | 87 | +0.0292 | 8.6 × 10⁻¹¹ |
| *Mycoplasma* | 86 | +0.0298 | 1.5 × 10⁻¹⁰ |

### S3 — a stricter baseline makes it slightly stronger

Requiring **three** callers to agree a gene is intact, rather than two, drops 1 106 baseline loci and
gives **+0.0319**, 104 % of base, same p. A more conservative baseline does not erode the gap.

### S4 — it strengthens as genomes are required to carry more loci

| minimum per class per genome | genomes | gap | p |
|---|---|---|---|
| 3 (registered) | 92 | +0.0307 | 6.2 × 10⁻¹² |
| 5 | 85 | +0.0310 | 8.1 × 10⁻¹² |
| 10 | 72 | +0.0307 | 5.8 × 10⁻¹² |
| 20 | 43 | +0.0335 | 2.2 × 10⁻¹⁰ |
| 50 | 21 | **+0.0393**, 21/21 positive | 9.5 × 10⁻⁷ |

Genomes with the most data show the largest gap, which is the direction a real effect should move.

### S5 — no one caller drives it

**The original version of this check did not test what it claimed**, as the review pointed out. It
removed only the loci where a caller is the *sole* pseudogene caller, which for Pseudofinder leaves
2 361 loci where it is still implicated. The real test removes **every** locus with that caller on the
pseudogene side:

| caller removed from the pseudogene side | loci dropped | genomes | gap | p |
|---|---|---|---|---|
| **Pseudofinder** (84.1 % of discordant loci) | **4 107** | 51 | **+0.0384** | 1.8 × 10⁻⁸ |
| PGAP (45.6 %) | 2 226 | 84 | +0.0386 | 3.3 × 10⁻¹³ |
| DFAST (40.7 %) | 1 987 | 78 | +0.0305 | 9.7 × 10⁻¹² |
| Bakta (2.9 %) | 143 | 90 | +0.0336 | 7.8 × 10⁻¹⁴ |

**G3 survives removing Pseudofinder from 84 % of the discordant set entirely**, and does so at a
*larger* effect than base. The conclusion is the one originally claimed, but the published check did
not establish it and the real one does — and it also disposes of the circularity worry, since
Pseudofinder is the caller derived from PGAP.

For completeness the sole-caller version, which is a weaker test:

| sole-caller loci removed | loci dropped | genomes | gap | p |
|---|---|---|---|---|
| Pseudofinder-only | 1 746 | 77 | +0.0287 | 1.4 × 10⁻⁸ |
| DFAST-only | 540 | 88 | +0.0302 | 1.2 × 10⁻⁹ |
| Bakta-only | 48 | 90 | +0.0306 | 7.8 × 10⁻¹⁴ |
| PGAP-only | 30 | 92 | +0.0307 | 6.2 × 10⁻¹² (= base; degenerate) |

### S6 — the one place it weakens

Requiring **two** callers to agree a locus is a pseudogene cuts the discordant set from 4 886 to 2 522
and the gap to **+0.0221, 72 % of base, p = 1.5 × 10⁻⁴**. Still positive and significant, but this is
the weakest arm and the honest reading is that part of G3's magnitude comes from single-caller calls.
Those are also where a caller's own error rate is highest, so this should be stated in the manuscript
rather than left to a reader to discover.

### S7 — not a property of the reduced genomes

Dropping any one reduction stratum leaves the gap between +0.0292 and +0.0329, all significant,
including dropping S1 extreme obligate entirely.

## Separately: the 26 % order-based contig mapping

`analysis_g3_contig_sensitivity.py`, reported in full there. The order fall-back **cannot choose
wrongly in a genome with one contig**, because `order[index]` is then the only contig there is.

| arm | genomes | gap | p |
|---|---|---|---|
| all | 92 | +0.0307 | 6.2 × 10⁻¹² |
| **single-contig only** (hazard impossible) | 68 | **+0.0265** | 1.3 × 10⁻⁸ |
| multi-contig only | 24 | +0.0363 | 2.8 × 10⁻⁴ |

Both arms hold independently with consistent effect sizes. And a fact that bounds the hazard further:
**not one of the 99 genomes has two contigs of equal length**, so the uniqueness test can only fail
when a caller's region length is *missing*, never because two contigs are indistinguishable by
length.

## Interpretation that must travel with G3

**Pseudofinder is not an independent caller.** It is a post-processor and was given PGAP's own
annotation as input, so PGAP and Pseudofinder cannot disagree about a gene PGAP called pseudogene
(`benchmark/PF_INPUT_DEPENDENCE.md`: in 155 of 160 genomes PGAP's pseudogenes carry a CDS with
`/pseudo`, which `annotate.py` converts straight to a pseudogene; verified on real output, 267 of 267
in `pf_intact.gff` = 0). The L2 asymmetry — 1 742 Pseudofinder-only against PGAP-intact, 0 the reverse
— is **structurally guaranteed and is not biological agreement.** Any text implying otherwise is
wrong, and this qualification belongs wherever the nesting is mentioned.

G3 itself does not depend on that pairing: it contrasts discordant loci against concordant-intact
loci regardless of which callers disagree, and S5 shows no caller carries it.

## What G3 does not establish

- It is an **evolutionary plausibility check, not a measure of caller precision or recall.** Elevated
  ω at discordant loci says selection has relaxed where callers disagree; it does not say which caller
  is right about any gene.
- dN/dS remains **one-sided**: high ω corroborates a pseudogene call, low ω does not refute one.
  **Superseded after review:** this overstates what ω supports. The FRAME2 control shows a frameshift
  **alone** raises measured ω by +0.0315, so an elevated value can come from the lesion's effect on the
  alignment rather than from relaxed selection. ω is an **association and plausibility measure with
  substantial methodological limitations**, not confirmation of a caller's accuracy. Corrected wording
  in `MANUSCRIPT_CANDIDATE_RBH_V2.md § D.1`.
- Overlap is not orthology. A locus is a 50 % reciprocal-overlap component, not a verified pair.
- The orthologue is a DIAMOND best hit, not a confirmed orthologue.

## S8 — the overlap threshold, the last outstanding check, now done

The 50 % reciprocal-overlap criterion was the one sensitivity check left open. Rebuilding the whole
analysis per threshold was infeasible — the validated run took 4.7 hours because it computes ω for
103 224 loci. But the threshold changes the **partition**, not the ω of a given span, so
`analysis_overlap_sensitivity.py` rebuilds only the partition, using `analysis_locus.py`'s own
union-find with `RECIPROCAL` overridden, and reuses the ω already computed for each locus it can match
back. Thirty-four seconds per threshold. **All outputs are isolated under
`results/overlap_sensitivity/`; `results/loci.tsv` and `results/locus.json` were not touched**, and
that is asserted after every run.

**Faithfulness control first.** At the base threshold of 0.50 the rebuild must recover the validated
locus set. It recovers **103 224 of 103 224 — 100.0 %** — and reproduces G3 to the digit
(+0.0307, 78/92, p = 6.16 × 10⁻¹²). So the reproduction is exact and the other rows mean something.

| reciprocal overlap | loci | discordant | base loci preserved | G3 on matched loci | p |
|---|---|---|---|---|---|
| 0.30 | 222 823 | 29 762 | 99.4 % | +0.0303 | 1.1 × 10⁻¹¹ |
| 0.40 | 226 649 | 31 709 | 99.9 % | +0.0303 | 1.1 × 10⁻¹¹ |
| **0.50 (base, control)** | 228 326 | 32 053 | **100.0 %** | **+0.0307** | 6.2 × 10⁻¹² |
| 0.60 | 229 049 | 31 338 | 99.6 % | +0.0286 | 3.3 × 10⁻¹¹ |
| 0.70 | 229 330 | 30 572 | 99.2 % | +0.0271 | 3.2 × 10⁻¹² |

**This sweep is conditional on the loci already scored at 0.50; it is not a complete recomputation.**
Every row above is G3 over the subset of the validated locus set that the rebuilt partition can match
back, with each locus keeping the ω computed in the validated run. Loci that are **new** at a
threshold have no ω and are dropped, and no ω is recomputed in the reading frame the new partition
would imply. The sweep therefore answers a narrower question than a full rebuild would: whether
re-partitioning the *same scored loci* moves the effect. It does not answer what G3 would be on a
locus set fully rebuilt and rescored at each threshold.

**On that conditional basis the threshold does not decide G3.** Across 0.30 to 0.70 the effect moves
between +0.0271 and +0.0307 — a span of 0.0036, about 12 % of the effect — and every threshold is
significant at better than 10⁻¹⁰, with 99.2–99.9 % of the base loci preserved. The 120 000–127 000
"new" candidate loci are nearly constant across thresholds, which says they are mostly the components
the validated run dropped for want of an orthologue hit — ω attrition rather than a threshold effect —
but that is an inference from their near-constancy, not a measurement of their ω.

## Files

`analysis_locus_sensitivity.py` → `results/locus_sensitivity.json` ·
`analysis_g3_contig_sensitivity.py` → `results/g3_contig_sensitivity.json` · input `results/loci.tsv`
