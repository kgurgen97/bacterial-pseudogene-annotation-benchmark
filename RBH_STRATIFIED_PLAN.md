# Plan: genome-stratified handling of differential RBH retention

**Frozen before the estimators were run.** `analysis_rbh_stratified.py` asserts this file's SHA-256 at
runtime and refuses to proceed if it differs. This plan is never edited because of a result; if a
result demands a different design, the remedy is a new, separately labelled plan.

## 1 The problem this addresses

`RBH_ORTHOLOGY_VALIDATION.md` reports that reciprocal-best-hit confirmation retains **95.0 %** of
concordant-intact loci and **60.4 %** of discordant loci — a **34.5 percentage point** imbalance — and
that the imbalance is concentrated: *Shigella* retains 33.8 % of its discordant loci at a 51.7 point
imbalance while supplying 43 % of all discordant loci, and *Sodalis* adds 382 more at 13.2 points.
Three of four scientific quality criteria failed, and retention imbalance was named the top-leverage
blocker.

G3 is **already** a median over per-genome gaps, so pooling across genomes is not the defect. The
defect is **within-genome**: in a genome that retains a third of its discordant loci, the retained gap
is computed over a non-random third, and § 4.2 of that document shows the excluded discordant loci are
systematically the shorter, higher-identity, lower-ω ones.

## 2 The property that makes this checkable rather than assumable

**ω is known for every locus, including every excluded one.** RBH exclusion removes loci from an *arm*,
not from the data. The quantity a selection correction is supposed to recover — the full-universe gap
E1 = **+0.0268** — is therefore already measured.

So this is not a missing-data problem that must be argued about. The correction can be **validated
against the truth it is trying to recover.** That is the design below, and it is the reason this plan
exists rather than a longer argument about assumptions.

## 3 Declared estimators

All four use the project's existing gap definition — per genome, median ω over discordant loci minus
median ω over concordant-intact loci, requiring ≥ 3 loci per side — combined across genomes by the
median of per-genome gaps, with genome-level and genus-level bootstrap intervals. ω is never
recomputed.

| id | estimator | purpose |
|---|---|---|
| **S0** | E1, all scored loci | the **truth** the corrections must recover |
| **S1** | E2 naive, RBH-retained loci, unweighted | the biased estimate, for comparison |
| **S2** | **IPW**: inverse probability of retention, fitted per genome and per locus class on log₁₀ span and alignment coverage, applied to retained loci only | the principled correction; its bias against S0 is the headline |
| **S3** | covariate-matched, 1:1 within genome, span and coverage only, caliper 0.10 | removes within-genome covariate imbalance among retained loci |
| **S4** | restricted to genomes whose retention imbalance is ≤ 10 points | the clean but low-power comparison |

### 3.1 S2 in detail

For each genome and each class separately, fit a logistic model of retention on
`[1, log10(span), forward_coverage]` by Newton–Raphson, with a ridge term of 10⁻⁶ on the slopes for
numerical stability. Weight each retained locus by `1 / p̂`, stabilised by dividing by the mean weight
within its (genome, class) cell so that weights average 1. Compute weighted medians per class, then the
gap, then the median across genomes.

Declared in advance:

- **Covariates are log₁₀ span and alignment coverage only.** Percent identity is excluded because it is
  close to a deterministic function of the same divergence ω measures; conditioning on it conditions on
  a mediator. This is the choice § 6.2 of the validation document argued for.
- **Weight truncation**: any stabilised weight above **10** is clipped to 10, and the clipped fraction
  is reported. Clipping is applied before the weighted median, never chosen after seeing the gap.
- **A cell is skipped** if either class has fewer than 10 retained loci or fewer than 3 excluded loci,
  because a retention model cannot be fitted without variation in the outcome.
- **Effective sample size** `(Σw)² / Σw²` is reported per cell; a cell whose ESS falls below 30 % of its
  retained count is flagged as weight-dominated.

### 3.2 Disclosure about S3

S3's configuration — span and coverage, caliper 0.10 — is the one cell of the ten-cell grid in
§ 6.3 of the validation document that met the balance criterion |SMD| ≤ 0.10, a criterion declared
before that grid was run. **It was nonetheless identified in a grid whose G3 values were also visible,
so S3 is a specification frozen after seeing those values, not a clean pre-registration.** It is
reported as a supporting arm, never as the headline, and this paragraph travels with it.

## 4 What each outcome would mean

Declared before running, so that no result can be re-labelled afterwards:

| observation | conclusion |
|---|---|
| **S2 recovers S0** within ±0.005 | retention is explainable by span and coverage; the retained estimate's bias is correctable, and differential retention is **not** a threat to G3 |
| **S2 lands between S1 and S0** | retention depends partly on the measured covariates and partly on something unmeasured — most plausibly ω itself. Partial correction only; the residual is the honest uncertainty |
| **S2 misses S0 by more than ±0.010, or overshoots** | retention is **not** explainable by span and coverage. The missingness is not at random on these covariates, IPW is invalid here, and the retained estimate cannot be repaired by weighting. This is a negative result and would be reported as one |
| **S1, S2, S3, S4 agree in sign and lie within ±0.010 of one another** | the sign-stability criterion passes |
| **S4 disagrees in sign with S1–S3** | the low-imbalance genomes behave differently from the high-imbalance ones, and the pooled retained estimate is driven by the genomes that drop the most |

A correction that merely moves the estimate in the comfortable direction is **not** success. The
criterion is agreement with S0, which is known, and overshooting counts as failure.

## 5 The diagnostic that interprets all of them

Per genome, the retention imbalance and the gap are both known. Their **Spearman correlation across the
12 genomes** is reported. A strong positive correlation would mean the genomes that discard the most
discordant loci are also the ones reporting the largest gaps, which would make the pooled retained
estimate a selection artefact. A null correlation is evidence against that, and is reported as such
without being called proof.

## 6 Gates

Technical gates are inherited unchanged from `analysis_rbh_audited.py`: ω bit-identical to its input,
input checksum stable across the run, and this plan's hash asserted. The scientific criterion this plan
exists to settle is **retention imbalance**, and it is declared to pass only if S2 recovers S0 within
±0.005 **and** no cell is flagged weight-dominated for more than 25 % of genomes.

Nothing here licenses the inference that verified orthology corroborates the biology. The selection
warning in § 7 of `RBH_ORTHOLOGY_VALIDATION.md` stands in full.

## 7 Scope

12-target pilot, 6 360 loci, input `results/rbh3/per_locus_rbh_audited.tsv`. Expansion to 92 targets
remains out of scope and is not authorised by this plan.
