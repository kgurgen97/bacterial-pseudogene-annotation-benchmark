# Is it the caller, or the configuration? A second Pseudofinder setting

**Written 2026-10-27, before `run_pf_configB.sh` or `analysis_pf_config.py` existed and before any
outcome was computed.** Nothing below may be changed because of a result. If a decision here turns out
to be wrong, the fix is a new, separately labelled analysis, not an edit to this file.

## 1. The question this answers, and why it is the first one a reviewer asks

Both frozen plans list the same limitation: **one configuration per caller.** Every result in this
project is reported as caller choice. The S1-versus-S4 sign reversal, the sixtyfold spread in S3
medians, the 3.68× spread in C1 coefficients, all of it is attributed to *which program was used*.

That attribution has never been tested. If changing one documented parameter inside a single caller
moves its numbers as far as switching to a different caller does, then "caller choice" is the wrong
frame and the paper has to say "configuration choice" instead. The claim would not get smaller, but it
would be about something else, and the current wording would be wrong.

## 2. Configuration B, fixed

Config A, already run on all 523 genomes, is:

    pseudofinder.py annotate --genome G --database sprot.dmnd --diamond --threads 4

with every other parameter at its default. Config B changes **one** parameter:

    --length_pseudo 0.50          (default 0.65)

A gene is called a truncation pseudogene when it is shorter than this fraction of the average length of
similar genes. Lowering it from 0.65 to 0.50 makes the caller **more conservative**: a gene must be more
severely truncated to qualify.

**Why this parameter and this value.** It is the one Pseudofinder's own help text uses as its worked
example (`"-l 0.50" will consider genes that are less than 50% of the average length of similar genes`),
so it is the change a user is most likely to make, and the value is the tool's own illustration rather
than one chosen here. Nothing else differs: same genomes, same database, same `--diamond`, same threads,
same counting rule (non-comment lines in `pf_pseudos.gff` against FASTA headers in `pf_intact.faa`).

Scope: all **523** genomes, so every comparison is on the identical corpus.

## 3. Pre-committed predictions, with failure criteria

Let **A/B fold** be the ratio of median pseudogene fraction between the two configurations, and
**between-caller fold** the ratio of the largest to the smallest median fraction across PGAP,
Pseudofinder-A, DFAST and Bakta on the same genomes.

| # | prediction | confirmed if | **failed if** |
|---|---|---|---|
| F1 | sanity: the more conservative setting calls fewer | median fraction B < median fraction A | B ≥ A |
| F2 | **the main test:** configuration moves it less than caller does | A/B fold ≤ half the between-caller fold | A/B fold ≥ the between-caller fold |
| F3 | C3's ordering survives the change | with config B, S3 > S1 and S3 > S4 both at p < 0.05 | either comparison loses significance |
| F4 | **the sharpest:** the S1-versus-S4 sign survives | config B keeps Pseudofinder's negative rank-biserial, endosymbionts lower | the sign flips to positive |
| F5 | C1 keeps its direction | the pooled Spearman stays positive | it turns negative |

**F2 is the one the framing rests on.** If it fails, the paper must be rewritten to attribute its spread
to configuration rather than to caller identity, and the title and abstract of any N2 manuscript change
with it. That is a real possibility and the reason this is being run before the manuscript is written.

**F4 is the one that would hurt most.** The S1-versus-S4 reversal is this project's sharpest result, and
it currently reads "two callers say lower, one says higher". If one documented parameter change flips
Pseudofinder's own sign, the reversal is not between callers at all; it is inside one of them.

## 4. What this cannot do

- One parameter, one alternative value. A grid over `--length_pseudo`, `--shared_hits`, `--evalue` and
  `--hitcap` would map the configuration space; this takes a single documented step through it, so a
  pass on F2 bounds nothing beyond this parameter.
- Only Pseudofinder is re-run. PGAP has no configuration to vary, being NCBI's fixed annotation, but
  Bakta and DFAST do, and their configuration sensitivity is untested.
- The corpus, counting rules and strata are unchanged, so everything the existing limits say still
  applies, including that S1 comparisons involving Bakta run against a floor.
- This tests sensitivity to a setting, not correctness. Neither configuration is validated against a
  ground truth, and nothing here says which is right.

## 5. Deliverables

`run_pf_configB.sh` → `~/research/datasets/n2-pseudofinder-B/` and `results/pf_configB_counts.tsv`;
`analysis_pf_config.py` → `results/pf_config.json`, `results/pf_config.tsv`.
Outcome in `PF_CONFIG_RESULTS.md`, reporting every verdict including failures.
