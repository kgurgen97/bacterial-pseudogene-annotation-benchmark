# Plan: the forward-ambiguity criterion, corrected

**Frozen before the analysis ran.** `analysis_forward_ambiguity.py` asserts this file's SHA-256 at
runtime and refuses to proceed if it differs. A frozen plan is never edited because of a result.

## 1 The defect in the existing criterion

`analysis_rbh_audited.py` declared a scientific criterion "forward ambiguity ≤ 10 %" and scored it
**pooled**: 347 of 6 360 loci, **5.46 %**, PASS. The class breakdown makes that pass meaningless:

| class | forward-ambiguous (margin < 1 %) |
|---|---|
| concordant intact | **1.9 %** |
| **discordant** | **13.9 %** |

Discordant loci are **seven times** more likely to have a near-tied forward runner-up, so the
orthologue assignment is least certain exactly where the analysis places its claim. The pooled
criterion passed only because baseline loci outnumber discordant loci four to one.
`RBH_ORTHOLOGY_VALIDATION.md § 8.2` recorded that a class-stratified version would fail, and § 13.7
listed it as the open item expansion would not fix.

## 2 Why the replacement is not a rate threshold

The obvious correction — the same threshold applied per class — would fail at 13.9 % against 10 %. But
a rate threshold is a **proxy**, and its relationship to the quantity anyone cares about is unknown: a
10 % ambiguity rate that does not move G3 matters less than a 2 % rate that does.

As with retention, **ω is known for every forward-ambiguous locus.** Ambiguity is a property recorded
per locus, not a reason the locus was dropped from the data. So the consequence of ambiguity is
directly measurable, and a proxy threshold is unnecessary. The criterion is therefore declared
**consequence-based**, with the class-stratified rates reported alongside as descriptive statistics
and the proxy's verdict reported as the failure it is.

**This is only legitimate because the direct measurement exists.** Replacing a failing proxy with a
measurable consequence is otherwise indistinguishable from choosing the criterion that passes, and
both verdicts are therefore reported side by side wherever this is summarised.

## 3 Declared quantities

Gap definition, combination across genomes and interval estimation are inherited unchanged from the
project: per genome, median ω over discordant loci minus median ω over concordant-intact loci,
≥ 3 loci per side, median of per-genome gaps, genome- and genus-clustered bootstraps. ω is never
recomputed.

| id | arm | purpose |
|---|---|---|
| **A0** | all scored loci | the **truth**, +0.0268 |
| **A1** | forward-unambiguous loci only | what dropping ambiguous loci does |
| **A2** | forward-ambiguous loci only | what the ambiguous loci themselves look like |
| **A3** | IPW on unambiguous loci, modelling ambiguity per genome and class from log₁₀ span, forward bitscore and log₁₀ candidate count | recovery of A0 from unambiguous loci alone |

A3 uses the machinery frozen under `RBH_STRATIFIED_PLAN.md`: logistic fit by Newton–Raphson, weights
`1/p̂` stabilised to mean 1 within the cell, capped at 10, effective sample size reported, and a cell
skipped if either side has too few loci to model.

## 4 Declared criterion

The corrected criterion has three parts, all declared now:

1. **Consequence**: |A1 − A0| ≤ **0.005**. Dropping every forward-ambiguous locus must not move G3 by
   more than that.
2. **No genome flips**: no genome's gap may change sign between A0 and A1.
3. **Rates are reported per class, per genome and per genus**, with no threshold attached, and the
   proxy criterion's class-stratified verdict stated explicitly.

Part 3 carries no pass or fail. Parts 1 and 2 do.

## 5 Declared mechanism test

Why are discordant loci seven times more ambiguous? Two explanations with different consequences, and
the analysis must separate them rather than assert one:

- **H1 paralogy.** Discordant loci are enriched in multi-copy families — transposases and insertion
  sequences are both pseudogene hotspots and large families — so several partner proteins are genuinely
  near-equivalent. The ambiguity is real and biological, and no aligner setting removes it.
- **H2 weak signal.** Discordant loci are shorter and more degraded, so their bitscores are smaller and
  the top two land closer together by chance. The ambiguity is a power problem.

Declared tests:

- **Adjustment.** A logistic model of ambiguity on log₁₀ span, forward bitscore and log₁₀ candidate
  count, with a class indicator, fitted pooled. The class coefficient **before and after** adjustment
  is reported. If adjustment removes it, H2; if it survives, H1 or something else.
- **Homology of the runner-up.** For ambiguous loci, the best and runner-up partner subjects are
  aligned to each other by DIAMOND `blastp`. If the runner-up is a strong hit of the best, the two are
  paralogues and the ambiguity is H1. The fraction is reported per class, with the bitscore of the
  best-versus-runner-up alignment relative to the best's self-score.
- **Candidate counts.** The median number of forward candidates passing e-value 10⁻⁵, by class.

### 5.1 What each mechanism would mean, declared in advance

| finding | conclusion |
|---|---|
| class effect vanishes after adjustment, runner-up rarely homologous | **H2**: a power problem. Ambiguity is noise, and arm A1 is the fix |
| class effect survives, runner-up usually homologous | **H1**: real paralogy among pseudogene-prone families. Not fixable by thresholds, and a limitation to state in the manuscript rather than a defect to repair |
| class effect survives, runner-up rarely homologous | neither; something unmodelled, and the honest report is that the asymmetry is unexplained |

## 6 Gate

Technical gates inherited: ω bit-identical to its input, input checksum stable, this plan's hash
asserted. The forward-ambiguity criterion passes only if parts 1 and 2 of § 4 both hold.

Nothing in this plan authorises expansion to 92 RBH targets, and nothing in it makes RBH biological
corroboration of G3.

## 7 Scope

12-target pilot, 6 360 loci, input `results/rbh3/per_locus_rbh_audited.tsv`.
