# Does the reduction gradient survive shared ancestry?

**Written 2026-10-30, before `analysis_phylo.py` existed and before any outcome was computed.**
Nothing below may be changed because of a result. If a decision here turns out to be wrong, the fix is
a new, separately labelled analysis, not an edit to this file.

`DNDS_PLAN.md`, `FRAME_PLAN.md` and `LOCUS_PLAN.md` are untouched.

## 1. The objection

C3 is the one claim that survives every caller and every denominator, so it carries the paper. It is
tested by comparing strata of genomes, and the genomes in a stratum are not independent draws. An
external audit put it plainly: similar median genome size between S3 and S4 does not control for
taxonomy, and a stratum dominated by one lineage reports that lineage's biology as though it were a
general pattern.

The corpus makes this concrete rather than hypothetical. The S1-plus-S4 contrast uses 95 genomes from
50 genera, and **one genus contributes 22 of them**. A test that treats those 22 as 22 independent
observations is overstating its evidence in exactly the way the pooled dN/dS tests did.

## 2. What is compared

**Lineage** comes from `results/corpus_traits.tsv`, which carries genus, family, order, class and
phylum. Genomes with no lineage are dropped and counted. Genus is the primary grouping because it is
the finest level with enough genomes per group; family is reported alongside as the coarser check.

**The contrasts** are the three the paper already makes: S3 against S4, S3 against S1, and S1 against
S4, each on the pseudogene fraction, for each of the four callers.

**Three ways of accounting for ancestry**, all fixed now:

1. **One genome per genus.** Draw one genome at random from each genus, run the contrast, repeat 2 000
   times, and report the proportion of draws reaching p < 0.05 and the median effect. This is the
   strictest reading and costs the most power.
2. **A genus-clustered bootstrap.** Resample whole genera with replacement and report the 95 %
   interval of the rank-biserial.
3. **A permutation restricted within genus.** Shuffle stratum labels only among genomes of the same
   genus, so a genus contributes to the null only through its own members. Genera whose genomes all
   share one stratum contribute nothing and are counted.

## 3. Pre-committed predictions, with failure criteria

| # | prediction | confirmed if | **failed if** |
|---|---|---|---|
| T1 | **the headline survives:** C3's S3 > S4 contrast holds with ancestry accounted for | for at least three of four callers, the genus-clustered interval excludes zero **and** the within-genus permutation gives p < 0.05 | fewer than two callers meet both |
| T2 | **one genome per genus:** the contrast is not carried by a single over-represented lineage | S3 > S4 reaches p < 0.05 in at least half of the 2 000 one-per-genus draws, for at least three callers | it does so in fewer than a quarter of draws for three or more callers |
| T3 | **a prediction of weakness:** the S1-versus-S4 contrast is the most lineage-dependent of the three, because S1 is the most taxonomically narrow stratum | S1 vs S4 loses significance under the within-genus permutation for more callers than S3 vs S4 does | it does not |
| T4 | **how much is there to correct?** the corpus is genuinely clustered | the largest genus contributes at least 10 % of any stratum used | it does not, in which case the correction was unnecessary and should be reported as such |

**T1 is the one that matters.** If it fails, C3 as stated is a statement about the lineages that
happen to be in this corpus and the manuscript must say so.

## 4. What this cannot settle

- Taxonomy is a proxy for phylogeny. A real tree from marker genes would be better and is not built
  here; genus and family are ranks, not branch lengths, and two genera can be far closer than two
  others.
- Dropping to one genome per genus discards most of the data, so T2 is conservative by construction
  and a failure there is weaker evidence against C3 than a failure of T1.
- Strata are lineage-level ecology assigned from the literature, as `CORPUS_EXPANSION_PLAN.md` states.
  Ecology and ancestry are correlated in nature, so a correction for ancestry removes some real
  ecological signal along with the confound. This analysis cannot separate the two and does not claim
  to.
- Nothing here addresses genome size or GC, which § 3.2 handles separately.

## 5. Deliverables

`analysis_phylo.py`, asserting this file's SHA-256 before computing, writing `results/phylo.json`
(per caller, per contrast: the clustered interval, the permutation p, the one-per-genus draw
proportion) and the counts behind T4. Outcome in `PHYLO_RESULTS.md`, reporting every verdict
including failures.
