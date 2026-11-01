# Plan: RBH orthology at full corpus scale, 92 targets

**Frozen before the expanded run.** Hash asserted by `analysis_rbh_expansion_check.py`. A frozen plan
is never edited because of a result.

## 1 Why expand, and what the pilot could not answer

`RBH_ORTHOLOGY_VALIDATION.md` § 15.7 lists the items still open after both blockers were measured. Two
of the three are consequences of pilot scale and nothing else:

- **One genome per genus.** All 12 pilot genomes sat in 12 distinct genera, so genus clustering was
  vacuous and every genus-clustered interval in § 1–15 is really a genome-level bootstrap. The corpus
  has **28 genera with up to 6 genomes each**, so only expansion can test within-genus structure.
- **Per-genome instability.** S2 (IPW) showed no per-genome improvement over the naive estimate — 4 of
  9 genomes, sign test p = 1.0 — and S3's per-genome values were the least stable of the four arms.
  Both are small-sample problems that more genomes address directly.
- **A3 underpowered.** Ambiguity-weighted recovery was fittable in 3 of 12 genomes.

The two blockers themselves are **not** reasons to expand: they were measured and found small
(retention +0.0022, ambiguity −0.0030). Expansion is to test whether those measurements hold at scale
and to obtain genus-level inference, not to improve them.

## 2 Scope and inputs

| | |
|---|---|
| input | `results/reference_sensitivity/per_locus_all.tsv`, **36 535 loci, 92 (query, partner) pairs, 28 genera, 4 strata** |
| control already established | its control column reproduced `results/loci.tsv` for **36 535 of 36 535** ω `[results/reference_sensitivity_all.json]` |
| pilot evidence | `results/rbh3/*` and `results/rbh_stratified/*`, **preserved**: every expanded run writes to its own directory |

## 3 Analyses, unchanged

The three analyses run exactly as frozen, on more data. Their own plans are adopted here by hash and
are not modified:

| analysis | plan | SHA-256 |
|---|---|---|
| audited RBH (§ 1–12 of the validation document) | none — gates declared in the script | — |
| genome-stratified retention | `RBH_STRATIFIED_PLAN.md` | `1fadc6db0ee4310b8004a9c488deb262d9d0815687423200edafa3e7fc77d4ef` |
| forward ambiguity | `FORWARD_AMBIGUITY_PLAN.md` | `e626f6d4386c4e824bd0872ead39ddee0206ef30c60ad6753a644365edb49698` |

Both of those plans declare a 12-target scope in their § 7. **Running them at 92 targets is work those
plans did not authorise**, which is why this plan exists and says so rather than quietly widening
theirs. No estimator, threshold, covariate or tolerance changes.

## 4 Declared primary comparisons

Pilot values are quoted so that agreement or disagreement is visible without reinterpretation.

| quantity | pilot, 12 targets | criterion at 92 targets |
|---|---|---|
| S0, the truth | +0.0268 | reported; it is the 92-target full-universe gap, a different population from the pilot's |
| retention imbalance, discordant vs baseline | 34.5 pp | reported per class, per genome, **per genus** |
| **consequence of retention**, S1 − S0 | **+0.0022** | **must stay within ±0.010 of the pilot's value** |
| **IPW recovery**, S2 − S0 same genome set | **−0.0009** | **\|S2 − S0\| ≤ 0.005**, the tolerance frozen in `RBH_STRATIFIED_PLAN.md` |
| **IPW per genome** | no improvement, 4/9 | **declared as the test this expansion exists for**: improvement over naive in more than half of genomes, sign test reported |
| forward ambiguity, discordant vs baseline | 13.90 % vs 1.94 % | reported per class |
| **consequence of ambiguity**, A1 − A0 | **−0.0030** | **\|A1 − A0\| ≤ 0.005** and no genome changes sign |
| runner-up homology among ambiguous loci | 100.0 % | reported; a value well below 100 % would mean the pilot's paralogy finding does not generalise |
| **genus-level inference** | vacuous, 1 genome per genus | **the new quantity**: genus-clustered interval on ≥ 20 genera with repeated genomes |

## 5 What would count as a failure

Declared now, so that no outcome can be relabelled:

- **The consequence measurements not holding.** If the retention consequence or the ambiguity
  consequence at 92 targets falls outside the bands in § 4, the pilot's reassuring numbers were a
  small-sample artefact and the manuscript may not cite them.
- **IPW still not improving per genome.** Then weighting is not a usable per-genome correction at any
  scale available here, and only the pooled statement survives.
- **Runner-up homology dropping well below 100 %.** Then the paralogy explanation is pilot-specific and
  the asymmetry returns to unexplained.
- **A genus-clustered interval including zero.** The manuscript's headline is a genus-level claim, so an
  RBH-retained genus interval that includes zero would be a substantive negative result and must be
  reported as the headline of this expansion, not buried.

A result that merely confirms the pilot is the expected outcome and is not evidence of anything beyond
the pilot's own scope. The informative outcomes are the genus-level interval and the per-genome IPW
test, because the pilot could not address either.

## 6 Guards

Host C: ≥ 10 GB to start, emergency stop below 8 GB, read from `/mnt/c` and never `df /`. Serial
execution: the audited run completes before the stratified and ambiguity runs start. Outputs isolated
under `results/rbh3_all/`, `results/rbh_stratified_all/` and `results/forward_ambiguity_all/`; ω is
never recomputed; `results/loci.tsv`, `results/locus.json`, `MANUSCRIPT.md` and every frozen plan are
untouched, and the inputs are checksummed before and after.

Nothing in this plan makes RBH biological corroboration of G3. The selection warning in § 7 of
`RBH_ORTHOLOGY_VALIDATION.md` stands at every scale.
