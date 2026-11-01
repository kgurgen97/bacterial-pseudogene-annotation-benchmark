# Annotation methods shape bacterial pseudogene counts and their interpretation: a comparison of 523 genomes

## 1 Title and authors

**Annotation methods shape bacterial pseudogene counts and their interpretation: a comparison of 523 genomes**

**Gurgen Karapetyan¹˒², Astghik Pepoyan¹**

¹ Armenian National Agrarian University, Yerevan, Armenia  
² Scientific Center of Zoology and Hydroecology, Yerevan, Armenia

**Corresponding author:** Gurgen Karapetyan

Running title: *Annotation-dependent bacterial pseudogene burden*

## 2 Abstract

**Background.** Bacterial pseudogene counts proxy genome decay, yet they are
produced by annotation pipelines whose criteria differ. The size of the resulting between-pipeline
variation has not been measured on a common corpus with a common locus definition.

**Results.** We annotated 523 bacterial genomes under four configurations: PGAP, Pseudofinder, Bakta
and DFAST. Corpus-total pseudogene counts span **19.1-fold** and median per-genome fractions
**77.2-fold**; one reports none in **128 of 523 genomes**. Pseudofinder, run on
PGAP-derived annotation as intended, is **not independent** of PGAP: 267 of 267 PGAP pseudogene CDS
appear in its pseudogene output, none in its intact output. Compared directly, by paired differences of
correlations on the same genomes, whether configurations disagree about the canonical
genome-reduction association **depends on the metric**: on raw counts two of three independent pairs
differ (PGAP − Bakta Δρ = +0.208, 95 % CI [+0.123, +0.297]), while on pseudogene fraction and calls
per Mb none is distinguishable. In 92 genomes with congeneric partners, **4 886 of 103 224**
caller-neutral loci are disputed and carry elevated dN/dS (+0.0302, 25 of 28 genera,
two-sided exact sign-test p = 2.744 × 10⁻⁵). In a small same-gene control (11 cases), frameshifts produced a median mechanical ω increase of **+0.0315**, but its magnitude was not robust under leave-one-out analysis. Independently of dN/dS, alignment evidence finds disruption at **33.3 % of disputed loci
against 3.5 %** of undisputed ones.

**Conclusions.** Reported pseudogene content depends substantially on the annotation configuration. In paired comparisons, no significant between-configuration differences were detected for normalised burden metrics, but these results do not establish equivalence. We make no claim about which configuration is correct or about gene function.

## 3 Keywords

bacterial genome annotation; pseudogenes; genome decay; dN/dS; annotation benchmarking;
reproducibility; PGAP; Pseudofinder; Bakta; DFAST

## 4 Introduction

Pseudogenes — genes disrupted by frameshifts, premature stops, truncation or fragmentation — are used
across bacterial genomics as a measure of genome decay, and their abundance underpins claims about
host restriction, recent niche change and relaxed selection (Kuo and Ochman, 2010; Goodhead and Darby,
2015; McCutcheon and Moran, 2012). In nearly all such work the count is taken from whatever annotation
pipeline produced the genome, and treated as a property of the organism.

Three recent studies bear directly on whether that treatment is safe, and none answers the question
posed here. Jundzill *et al.* (2026) benchmarked prokaryotic annotation tools across thousands of
species, including a large frameshift-perturbation arm, but measured annotation **stability** —
coding density, feature counts, feature length, GO coverage — not pseudogene detection; pseudogenes
are not an endpoint and the word appears once, in a list of feature types. Cooley and Wright (2024)
showed that many purported pseudogenes are artefacts of assembly error, using PGAP alone and counting
per assembly rather than per gene. Syberg-Olsen *et al.* (2022), introducing Pseudofinder, included a
simulation module whose validation is, by the authors' own description, circular: the truth set was
produced with the same parameters used to detect it, on one genome, without quantitative accuracy.

What is missing is a **locus-level** comparison of pseudogene calls across annotation configurations on
one corpus, with a caller-neutral definition of "the same gene", and with explicit controls for the
ways such a comparison can mislead. That is what we provide. We ask three questions:

**Q1.** How large is the between-configuration variation in pseudogene counts on a fixed corpus, and
which measurement of "spread" is being reported?

**Q2.** Does that variation change the conclusion of a downstream comparative analysis?

**Q3.** At loci where configurations disagree, is there independent sequence-level evidence that
distinguishes them — and what are the limits of such evidence?

We state at the outset what this study does **not** do. It does not rank the configurations by
accuracy: that would require per-gene ground truth at corpus scale, which we do not have (§ 13.2).

## 5 Research questions and hypotheses

| id | hypothesis | pre-specified | outcome |
|---|---|---|---|
| H1 | Pseudogene counts differ substantially between configurations on a fixed corpus | yes | supported, § 9.1 |
| H2 | The configuration changes the sign or significance of a downstream comparative result | yes | **metric-dependent**: on raw counts the configurations differ significantly (2 of 3 independent pairs); on **normalised** metrics no significant paired contrasts were detected, although the individual significance verdicts differ; this does not establish equivalence. The sign never changes. § 9.3 |
| H3 | Loci where configurations disagree show elevated dN/dS against loci all configurations call intact | yes | supported with limitations, § 9.4 |
| H4 | Elevated dN/dS at disputed loci is attributable to relaxed selection rather than to the lesion's effect on alignment | yes | **not supported as stated**, § 9.5 |
| H5 | The pseudogene-fraction gradient across ecological strata survives phylogenetic clustering | yes, `PHYLO_PLAN.md` | partly; its compound criterion failed as written, § 9.7 |

## 6 Dataset and genome selection

523 bacterial genomes with complete RefSeq assemblies, stratified by ecology into four strata: S1
extreme obligate endosymbionts, S2 obligate intracellular, S3 recent host restriction, S4 free-living
breadth `[results/corpus523.tsv, results/c3_strata.json]`. 363 genomes carry the trait annotations
required for the downstream re-test (§ 9.3) `[results/literature_retest.tsv]`.

**Figure 1** summarises the design and states its boundaries.

Locus-level dN/dS requires a congeneric partner genome; **92 genomes in 28 genera** have one and form
the dN/dS subset `[results/dnds_pairs.tsv]`.

**Two locus universes are used and must not be mixed** (§ 12.1): the full universe of **103 224**
scored loci, and a **36 535**-locus reference-sensitivity subset in which the concordant-intact
baseline is capped at 400 loci per genome.

## 7 Annotation methods

| configuration | version / mode | input |
|---|---|---|
| NCBI PGAP | as distributed in the RefSeq record | assembly |
| Pseudofinder | default, `annotate` | **PGAP-derived annotation** |
| Bakta | 1.12.2 | assembly |
| DFAST | 1.5.3 | assembly |

**Pseudofinder is a post-processor, not a fourth independent pipeline.** It is run as its authors
intend, on an existing annotation, and we report it as a configuration rather than as an independent
observer. The consequence is measured in § 9.2.

## 8 Definitions, counting rules and statistical methods

### 8.1 What counts as a pseudogene call

Each configuration asserts pseudogenes differently, and we use each one's own assertion rather than a
harmonised definition:

- **PGAP**: a `gene` feature carrying `/pseudo` or `/pseudogene` with no CDS beside it. Looking only at
  CDS features finds none of them.
- **Pseudofinder**: every entry in `pf_pseudos.gff`.
- **Bakta**: two keys. `pseudogene=` matches its own summary line (6 844 corpus-wide); `pseudo=True`
  marks a further 509. **The strict `pseudogene=` definition is primary** throughout; inclusive counts
  are reported separately wherever they change a conclusion.
- **DFAST**: `pseudogene_summary.tsv` rows classified `frameshift`, `internal_stop_codon` or both.
  `partial` is **not** a pseudogene call.

A defect in an earlier version of our own parser matched the substring `pseudo` anywhere in the GFF
attributes column, producing 125 307 spurious DFAST "pseudogenes". It was found, fixed and verified
against each runner's independent counts in all 523 genomes `[benchmark/SUBSTRING_CLASSIFIER_DEFECT.md,
benchmark/results/classifier_fix.json]`. All results here are post-fix.

### 8.2 Denominators

Per-genome pseudogene **fraction** is calls divided by that genome's own CDS count from the same
configuration `[results/denominator.json]`. Corpus **totals** are the sum of calls. These are different
measurements and § 9.1 reports both.

### 8.3 Caller-neutral loci

Intervals from all four configurations are joined into connected components by **50 % reciprocal
overlap**, with the longest interval in a component as its representative. A locus is **discordant**
when at least one configuration calls it a pseudogene and at least one calls it intact, and
**concordant-intact** when all four call it intact `[results/locus.json]`.

### 8.4 dN/dS

Nei–Gojobori with Jukes–Cantor correction, in the orthologue's reading frame. The orthologue is the
best DIAMOND `blastx` hit in a congeneric partner genome (Buchfink *et al.*, 2015). **dN/dS is treated
as one-sided and as an association**: a high value is consistent with a disrupted gene but does not
confirm a call, and § 9.5 shows why.

### 8.5 Inference and units

The gap statistic ("G3") is, per genome, the median ω over discordant loci minus the median over
concordant-intact loci, requiring ≥ 3 loci per side, combined across genomes as the median of
per-genome gaps.

**Two units are used and are always labelled.** Genome-level sign tests treat the genome as
independent; **genus-clustered bootstraps** resample genera and are the headline, because genomes
within a genus are not independent. Where a genus-clustered interval is quoted, a **genus-level**
directional test is quoted with it. The per-genome p-value for the headline gap is 6.2 × 10⁻¹², the
genus-level p-value is 2.7 × 10⁻⁵, and **the genome-level figure overstates the evidence by roughly
seven orders of magnitude** `[results/review_response.json]`.

An earlier version of the analysis drew intact control sets from a single module-level random
generator shared across every genome and configuration, so a change in one configuration's gene count
moved every later draw. This was found and removed; the primary estimator now uses **every** intact
gene with no subsampling and carries no seed, and a 40-seed sweep of the subsampled variant is reported
`[results/seed_sensitivity.json]`.

## 9 Main results

### 9.1 Pseudogene counts span 19- to 77-fold, and the two spreads measure different things

| configuration | corpus total | median per-genome fraction | mean fraction | genomes with zero calls |
|---|---|---|---|---|
| Bakta (strict) | 6 844 | 0.0008225 | 0.0053 | **128** of 523 |
| PGAP | 37 493 | 0.0111 | 0.0228 | 5 |
| DFAST | 45 976 | 0.0123 | 0.0295 | 3 |
| Pseudofinder | 130 748 | 0.0634921 | 0.0721 | 0 |

**19.1-fold** on corpus totals, **77.2-fold** on unrounded median per-genome fractions (0.0634920635 / 0.0008225037)
`[benchmark/results/two_spreads.json]`. **Figure 2.**

The two numbers are not two estimates of one quantity. A corpus total is dominated by pseudogene-rich
genomes; a median fraction is dominated by the typical genome. Almost the whole 4-fold difference
between the two ratios is Bakta's 128 zero-call genomes, which pin its median near a floor while its
mean is 6.5 times higher.

### 9.2 Pseudofinder's agreement with PGAP is structural

Of one genome's 267 PGAP pseudogene CDS, **267 appear in `pf_pseudos.gff` and 0 in `pf_intact.gff`**
`[benchmark/results/pf_input_dependence.json]`. Consequently the 1 742-to-0 asymmetry between
Pseudofinder-pseudogene/PGAP-intact loci and the reverse is **guaranteed by construction and is not
gene-by-gene agreement**. Any analysis treating these two configurations as independent observers —
including inter-annotator agreement statistics — is invalid.

### 9.3 Annotation choice changes the magnitude and the significance of a canonical result, not its sign

We re-tested the canonical genome-reduction claim — pseudogenes accumulate during reduction, so
pseudogene fraction should be **higher in smaller genomes**, i.e. a **negative** correlation with
genome size (Moran and Plague, 2004) — under each configuration and each count metric, on 363 genomes,
with partial correlations controlling for GC content. Design frozen in advance
(`LITERATURE_RETEST_PLAN.md`, sha256 `a576cda9…`).

Pooled, pseudogene **fraction** against genome size `[results/literature_retest.json: C1]`:

| configuration | Spearman ρ | p | partial ρ given GC | p | reading |
|---|---|---|---|---|---|
| PGAP | **+0.126** | 0.0163 | +0.112 | 0.0333 | significant, **opposite** |
| Pseudofinder | **+0.233** | 7.0 × 10⁻⁶ | +0.154 | 0.0033 | significant, **opposite** |
| Bakta | +0.097 | 0.0654 | +0.126 | 0.0164 | raw not significant, partial opposite |
| DFAST | +0.063 | 0.228 | **−0.049** | 0.354 | not significant either way |

**First, a finding about the claim, not about the pipelines: 0 of 4 configurations find a significant
association in the predicted direction, and all four point the same way.** The canonical claim does not
hold in this corpus under any configuration. Across all 48 caller × metric × stratification tests, 30
are significant in the opposite direction and 18 are non-significant. **Figure 3.** That is a result
about this corpus and that claim, and it is **not** evidence that annotation choice changes a
conclusion, because the configurations agree.

**Second, the configurations must be compared directly, and a verdict comparison does not do it.** A
significant result under one configuration and a non-significant one under another is not evidence
that the two effects differ. We therefore computed the **paired difference of correlations on the same
genomes**, Δρ = ρ_A − ρ_B, with a genus cluster bootstrap (B = 5 000) recomputing both correlations in
every replicate, Holm adjustment across the six pairs per metric, and pairs involving Pseudofinder
quarantined because it post-processes PGAP output. 362 genomes, 185 genera. Design frozen before
computation (`CALLER_EFFECT_CONTRAST_PLAN.md`, sha256 `0d1d6458…`)
`[results/caller_effect_contrast.json, DIRECT_CALLER_EFFECT_COMPARISON.md]`.

**On the primary metric, pseudogene fraction, no significant paired magnitude differences were detected among the independent configurations:**

| independent pair | Δρ | 95 % CI | Holm p | |
|---|---|---|---|---|
| PGAP − Bakta | +0.029 | [−0.073, +0.130] | 1.000 | verdicts differ; **effect-size contrast not statistically significant** |
| PGAP − DFAST | +0.063 | [−0.013, +0.134] | 0.312 | verdicts differ, **magnitudes do not** |
| Bakta − DFAST | +0.033 | [−0.080, +0.143] | 1.000 | no difference detected |

So although a reader using PGAP would report a significant relationship (p = 0.016) and one using
DFAST would not (p = 0.228), **those two correlations are not statistically distinguishable from each
other.** We say so rather than resting the claim on the threshold crossing.

**On raw counts, two of three independent pairs do differ:** PGAP − Bakta Δρ = **+0.208**
[+0.123, +0.297], Holm p = 6 × 10⁻⁴; Bakta − DFAST Δρ = **−0.188** [−0.279, −0.098], Holm p = 6 × 10⁻⁴;
PGAP − DFAST +0.021 [−0.032, +0.072], not distinguishable. Here all four configurations are
individually significant, so the **verdicts agree while the magnitudes genuinely differ** — a
difference a verdict comparison would have missed entirely.

**The pattern across metrics is the result.**

| metric | normalisation | independent pairs with distinguishable magnitudes |
|---|---|---|
| **count** | none | **2 of 3** |
| fraction | ÷ the configuration's own CDS count | 0 of 3 |
| per Mb | ÷ genome length | 0 of 3 |

**Whether the configurations appear to disagree about how pseudogene burden scales with genome size
depends on the metric.** On raw counts they demonstrably differ, unsurprisingly given that they emit
between 6 844 and 130 748 calls corpus-wide. Normalising by each configuration's own gene count, or by
genome length, yielded no statistically significant differences in the tested independent-pair contrasts. This does not establish equivalent scaling across tools. This is a narrower claim than "annotation choice changes the conclusion", and it is the one
the data support: it locates where the pipeline-dependence lives and where it does not.

*(The metric pattern itself was not pre-specified as a hypothesis — the three metrics were, the pattern
across them is a post-hoc observation and is labelled as such.)*

### 9.4 Loci where configurations disagree carry elevated dN/dS

103 224 scored loci; **4 886 discordant**, 98 338 concordant-intact. Discordant loci sit above
concordant-intact loci by **+0.0302 at genus level, 25 of 28 genera positive, two-sided exact sign-test p = 2.744 × 10⁻⁵**,
genus-clustered bootstrap CI **[+0.0206, +0.0401]** `[results/review_response.json, results/locus.json]`.
Per genome the gap is +0.0307 over 92 genomes, 78 positive — reported as secondary (§ 8.5).

The association is stable under every structural choice we varied: reciprocal-overlap threshold 0.30
to 0.70 (conditional on the loci already scored), baselines requiring three configurations, removing
every locus with Pseudofinder on the pseudogene side (+0.0384), single-contig genomes only, minimum
loci per genome raised 3 → 50, and leave-one-genus-out `[results/locus_sensitivity.json,
results/overlap_sensitivity.json, results/pf_residual_dependence.json,
results/g3_contig_sensitivity.json]`. **Robustness here is statistical, not biological**: these checks
re-slice the same ω values and cannot detect a bias inside the estimator. § 9.5 looks for one.

### 9.5 Part of the elevation is mechanical; the small synthetic arm limits precision

Using simulated genomes as an exact control — the same gene with and without a known lesion, against
the same orthologue, through the same pipeline, with the design frozen in advance
(`FRAME2_PLAN.md`, sha256 `a926ec3c…`):

| arm | n | median Δω |
|---|---|---|
| untouched genes (determinism check) | 225 | **0.0000, all 225 exactly zero** |
| verified synonymous substitution | 3 | −0.0039 |
| **single frameshift** | 11 | **+0.0315** (sign p = 0.0215, CI [+0.0017, +0.1432]) |
| premature stop, frame intact | 3 | **+0.0001** |

`[results/frame2.json]`. **Figure 4a.** The 11-gene frameshift median is near the prespecified +0.030 decision boundary: leaving out one gene changes it to +0.0246 to +0.0370, with 5/11 leave-one-out estimates below +0.030. The bootstrap confidence interval [+0.0017, +0.1432] is broad. Evidence supports a mechanical effect, not a precise magnitude. The same orthologue was selected in both assemblies for 268 of
268 genes, and the inflation does not track alignment loss (Spearman ρ = −0.018, p = 0.958)
`[benchmark/results/frame2_confound.json]`.

**A synthetic frameshift can induce a positive ω shift of the same order as the observed gap, but the magnitude is uncertain with only 11 cases.** Stratifying disputed
loci by lesion class — possible only for DFAST, the sole configuration publishing a per-call class —
separates the contributions: frameshift-bearing calls +0.0502 (278 loci, 29 genomes, p = 1.5 × 10⁻⁵)
against internal-stop-only calls **+0.0220** (196 loci, 25 genomes, p = 9.1 × 10⁻⁴)
`[results/g3_by_lesion.json]`.

**What this does and does not establish.** Within DFAST's calls, a gap survives in the stratum carrying
no frameshift, so the artefact is not the whole story there. It covers **539 of 4 886 discordant loci
(11.0 %)** and only DFAST's; Pseudofinder's 84.1 % share is unclassified. **It does not establish an
artefact-free effect across configurations**, and no pooled ω gap should be quoted alone.

The frame-artefact analysis originally planned under `FRAME2_PLAN.md` **failed its own pre-registered
gate F4** and the conclusion it would have supported was withdrawn; what is reported above is the
post-hoc stratification run afterwards, labelled as such `[FRAME2_RESULTS.md]`.

### 9.6 Elevated dN/dS does not adjudicate between configurations

All four configurations place their own pseudogene calls above their own intact genes, but only
Pseudofinder's genus-clustered interval excludes zero; the other three lose interval-level significance
once genus structure is respected `[results/seed_sensitivity.json]`. A matched within-genome comparison
of Pseudofinder against DFAST gives **+0.0072, CI [−0.0040, +0.0142], p = 0.296** — **absence of
evidence, not equivalence**: no equivalence margin was pre-specified and power is 0.17–0.60
`[results/caller_comparison.json]`.

### 9.7 A known ecological gradient recovers under all four configurations, with caller-dependent strength

This section is included for what it says about the **configurations**, not as a biological finding:
the gradient it uses is already established in the literature, and we reproduce it as a second
downstream analysis on which the configurations can be compared.

Pseudogene fraction is higher at recent host restriction than in free-living lineages (stratum S3
against S4), and the genus-clustered rank-biserial interval **excludes zero for all four
configurations** — but the **strength differs markedly between them**: PGAP +0.911 [+0.791, +0.981]
and DFAST +0.778 [+0.601, +0.921] against Pseudofinder +0.432 [+0.064, +0.732] and Bakta +0.330
[+0.047, +0.587]. Resampling one genome per genus reaches significance in **100 % of draws for PGAP
and DFAST but 53 % and 42 % for Pseudofinder and Bakta** `[PHYLO_RESULTS.md]`.

**That is the point for this study.** A reader using Bakta would recover the same gradient at less
than half the per-draw reliability of a reader using PGAP. Unlike § 9.3 these configurations are not
formally contrasted against one another, so we claim only that the strengths differ **descriptively**.

**Its pre-registered criterion T1 is recorded as FAILED as written.** T1 was compound and its second
half — a within-genus permutation — is **not testable on this corpus**: no genus contains genomes from
two different strata, so the count of informative genera is zero in all three contrasts and every
permutation reproduces the observed data. The half that could be tested passes four times out of four.
The S1-against-S4 sign reversal does **not** survive: its intervals include zero for three of four
configurations, and it is reported as exploratory only.

## 10 Orthologue and dN/dS sensitivity analyses

Both analyses below use the **36 535-locus subset** (92 genomes, 28 genera), not the 103 224-locus
universe.

### 10.1 Reference selection

Every ω uses an orthologue from a reference that **excludes a mean 9.2 % (maximum 35.0 %) of partner
genes** by Pseudofinder's own judgement — the caller under test also builds the measuring instrument.
Replacing it with each partner's own RefSeq annotation: the control arm reproduces the published
pipeline for **36 535 of 36 535 ω**, orthologue identity agrees at 99.1 %, and the gap moves from
+0.0301 to **+0.0266**, a contribution of **−0.0035 (~12 %), paired sign test p = 5.1 × 10⁻⁴**, with
80 of 92 genomes positive and the genus interval still excluding zero
`[results/reference_sensitivity_all.json]`. **Figure 4b.** The contribution is small, real and
measured; the result is **not** reference-independent.

### 10.2 Orthologue assignment by reciprocal best hits

Replacing the one-way best hit with a reciprocal best hit plus synteny: **93.2 %** of loci are
RBH-confirmed (94.8 % under textbook protein-to-protein RBH where an annotated CDS exists), and the gap
among confirmed loci is **+0.0291, genus-clustered CI [+0.0226, +0.0339], 24 of 28 genera positive,
genus-level sign p = 9.0 × 10⁻⁵** `[results/rbh_audited_all.json, results/review_fixes/review_fixes.json]`.
**Figure 4c.** On a paired within-genome comparison over the same 90 genomes, filtering **does not
detectably change the gap** (median paired difference +0.000000, sign p 0.368; four estimands span
0.000 to −0.004).

**This is an orthologue-assignment sensitivity analysis, not biological validation**, and its
limitations are in § 13.3.

## 11 Independent validation and its limitations

### 11.1 What was attempted

| attempt | outcome |
|---|---|
| independent simulator, Stage A | ran; 65–73 % of truly lesioned genes are annotated intact `[benchmark/results/stage_a_evaluation.json]` |
| independent simulator, Stages B and C | **never ran** — forbidden by a hard disk-capacity constraint on the analysis workstation |
| pre-registered gate **B1** | **FAILED**; honoured without moving the threshold `[benchmark/results/b1_diagnosis.json]` |
| caller-independent lesion classification from annotation-bounded intervals | partial: decides 43.5 % of discordant loci with **zero** false frameshift calls, but recovers only **18 of 413 (4.4 %)** known frameshifts `[results/lesion_from_sequence.json]` |
| caller-independent lesion classification from **extended** intervals | § 11.2 |
| inverse-probability weighting to correct differential RBH retention | **FAILED** both declared tests |
| covariate matching to correct it | **FAILED**; no grid configuration balanced at 92 genomes |

### 11.2 Sequence-based lesion evidence on extended intervals (not independent functional ground truth)

Design frozen before the code existed (`LESION_HARMONISATION_PLAN.md`, sha256 `27552b9d…`), with six
gates and the failure path fixed in advance. 6 360 loci in 12 genomes: 1 871 discordant, 4 489
concordant-intact `[results/lesion_extended.json]`.

**The problem addressed.** Our earlier caller-independent classifier, run on annotation-bounded
intervals, recovered only **18 of 413 (4.4 %)** known frameshifts. The mechanism is that a caller
detecting a frameshift reports the *truncated* feature, so the post-lesion fragment lies outside the
interval and the two-frame signature cannot be seen. Each span was therefore extended by
**max(300 bp, half the orthologue length)** before evidence was collected.

**Controls first, named for what they are.** None of these is independent biological ground truth, and
we do not present them as such (`INDEPENDENT_VALIDATION_AUDIT.md`).

| control | n | result | what it establishes |
|---|---|---|---|
| **artificial-lesion detection control**: a 1 bp deletion seeded at 40 % of the gene | 96 | **96.9 %** classified `frameshift` | that the classifier detects a frameshift **it is known to contain**. **Not** real-world biological sensitivity |
| the same genes unmodified | 96 | 7.3 % show any lesion | background rate of the same pipeline with no lesion present |
| **annotation-consensus control**: loci all four configurations call intact | **4 489** | **4.5 % show any lesion**, Wilson [3.9, 5.1] | that lesions are rarely reported where four pipelines agree there is none. These are **not independently verified true negatives** |

The annotation-consensus control is the one the earlier analysis lacked entirely, so its specificity
was unmeasured: a classifier that labels intact genes as lesioned is useless regardless of its recall.
Its limitation is that its negatives are a consensus of the same four pipelines under study.

**Against DFAST's published per-call class** (806 labelled loci):

| metric | annotation-bounded (earlier) | **extended** |
|---|---|---|
| frameshift recall | 18/413 = **4.4 %** | **307/466 = 65.9 %** |
| false frameshift on internal-stop loci | 0 % | **0.9 %** (3/340) |
| internal-stop recall | 280/340 = **82.4 %** | **224/340 = 65.9 %** |
| discordant loci decided | 43.5 % | **84.6 %** |

**Extension buys frameshift sensitivity and costs internal-stop sensitivity.** Frameshift recall rises
fifteen-fold; internal-stop recall falls, because for 90 of 340 internal-stop loci the extended span
recovers a full-length alignment past the stop and the locus is classified `intact`. This is a real
trade-off, not a tuning artefact, and it is why the plan's gate **G5 (internal-stop recall ≥ 70 %,
no regression) FAILS at 65.9 %**. Under § 12 of the plan that blocks expansion to 92 genomes; under
§ 11 the method remains usable, because the gates that govern usability (G1 frameshift recall ≥ 40 %,
G3 false-lesion rate ≤ 10 %) both pass.

**The result the classification was for.** Lesion evidence is far more common at disputed loci than at
loci all four configurations call intact, per genome and with a genus-clustered interval:

| | median per-genome lesion rate |
|---|---|
| discordant loci | **33.3 %** |
| concordant-intact loci | **3.5 %** |
| **difference** | **+31.4 percentage points**, 12 of 12 genomes, sign p = 4.9 × 10⁻⁴, genus-clustered CI **[+6.5, +46.2] pp** |

`[results/lesion_extended.json: lesion_rate_contrast, results/lesion_extended/lesion_rate_by_genome.tsv]`.
**Figure 5.**

**This is sequence-based evidence of lesions that does not use dN/dS, not independent biological validation of gene function.** It is caller-
independent — the evidence is alignment geometry against an orthologue, and DFAST's class is used only
afterwards, for agreement — and it is not vulnerable to the frameshift artefact of § 9.5, because it
does not estimate a substitution rate.

**It is sequence evidence, not evidence of lost function.** A frameshift or premature stop detected by
alignment geometry is a statement about the sequence, not a demonstration that the gene is
non-functional: proteogenomic work in other organisms shows that loci called pseudogenes can still be
expressed. We make no loss-of-function claim anywhere in this study, and § 11.4 records why
independent functional ground truth was not obtainable here.

**What it does not show.** Among discordant loci, **41.2 % still classify as `intact`**: on extended
evidence, four in ten disputed loci show a full-length, single-frame, stop-free alignment to their
orthologue. The comparable figure for undisputed intact genes is 72.7 %, so disputed loci are far less
often intact — but a substantial minority carry no detectable lesion, consistent with Cooley and
Wright (2024). We do not claim these are false positives of the calling configuration: absence of
detectable evidence under this method is not evidence of absence, and 22.8 % of undisputed intact genes
are themselves `ambiguous` here.

Per the frozen plan, **caller disagreement is not decomposed into definitional, detection and boundary
components**: that decomposition would require the internal-stop arm to hold up, and it did not.

### 11.3 Why independent biological ground truth was not obtainable

We audited the literature and public resources for pseudogenes supported independently of the four
configurations, requiring traceable provenance, verified strain identity and lesion coordinates
(`INDEPENDENT_VALIDATION_AUDIT.md`). **No benchmark-grade set exists for this corpus**, for three
specific reasons:

1. **No curated database of experimentally verified bacterial pseudogenes was found.** PHI-base and
   PATRIC curate virulence genes, not lesion-level pseudogene truth.
2. **The curated sets that do exist are for organisms we cannot evaluate.** Manually curated
   *Salmonella* pseudogene lists, and the proteogenomic evidence on pseudogene expression in
   *Salmonella* and *E. coli*, concern genomes with **no congeneric partner in this corpus** — the
   corpus holds one *Salmonella* and no *Escherichia*, so neither is in the 92-genome subset where
   caller-neutral loci exist. Evaluating them would require annotating further genomes.
3. **The best experimentally grounded cases are deletions, which have no locus to classify.** The
   *Shigella* antivirulence genes — *cadA*, *nadA*, *nadB*, *speG*, *ompT* — carry genuine
   loss-of-function evidence from complementation: reintroducing *cadA* attenuates *S. flexneri*
   virulence. But these are **large genomic deletions**. A deleted gene is absent from the assembly:
   there is no interval, no alignment geometry, and nothing for any method here to classify. The
   strongest available experimental evidence is therefore outside the evaluable class, by the nature
   of the biology rather than of the method.

**Published pseudogene lists for these genomes are excluded as circular**: *M. leprae*'s ~1 100,
*Y. pestis* CO92's 149 and *Shigella*'s few hundred each come from that genome's own annotation, so
using them would evaluate one annotation against another.

We therefore make **no loss-of-function claim**, and we label every control in § 11.2 for what it is:
an artificial-lesion detection control and an annotation-consensus control.

### 11.4 Why dN/dS cannot serve as the validator

The lesion that defines a pseudogene is also a perturbation of the alignment from which dN/dS is
computed (§ 9.5). Using ω to decide which configuration is right is therefore circular for exactly the
class of call — frameshifts — that dominates disagreement.

## 12 Discussion

### 12.1 The headline is a measurement problem, not a tool ranking

A pseudogene count, as reported in a genome record, is substantially a property of the pipeline.
Studies using pseudogene counts across genomes annotated by different pipelines, or compared against
literature values obtained under a different pipeline, are comparing partly incommensurable quantities.

The practical consequence (§ 9.3) is narrower than "conclusions reverse", and we state the narrow
version because it is the one the direct comparison supports. All four configurations agree on
direction, so the configuration does not flip the sign. **On the normalised metrics the configurations
are not even distinguishable from one another**, despite their significance verdicts differing — a
verdict comparison would have suggested otherwise, and we report the paired contrast instead. **On raw
counts they do differ**, substantially and after adjustment.

The practical guidance follows directly: **report pseudogene burden as a fraction of that pipeline's
own gene count rather than as a raw count**, and the analysed correlations do not show statistically detectable differences between independent configurations on normalised metrics. What remains — a 77-fold spread in the fraction itself, and one pipeline reporting nothing
in 128 of 523 genomes — is a difference in level, not in how burden scales with genome size.

We deliberately do **not** conclude which configuration is best. The configurations encode different
definitions — Bakta's strict `pseudogene=` key against Pseudofinder's inclusive post-processing — and a
definitional difference is not an accuracy difference.

### 12.2 What the dN/dS result means

Disputed loci do carry elevated ω, at genus level, robustly to our analytical choices, to
reference selection and to verified orthology. But ω is an **association and plausibility measure with
substantial methodological limitations**, not confirmation of a call: a small synthetic frameshift arm produced an effect of similar order, but its magnitude was uncertain (11 genes, with 5 of 11 leave-one-out medians below the prespecified threshold). The useful statement is narrow and worth having — **disagreement between
pipelines is not random with respect to evolutionary signal** — and the unsupported statement is the
obvious one, that elevated ω shows the pseudogene-calling configuration was right.

### 12.3 Relation to prior work

Our results are consistent with Cooley and Wright (2024) in direction — many pseudogene calls are
fragile — while differing in cause: they isolate assembly error with one annotator, we isolate the
annotator with fixed assemblies. They are orthogonal to Jundzill *et al.* (2026), who measure
annotation stability under perturbation rather than pseudogene detection, and whose finding that one
pipeline is most **resilient** to introduced frameshifts cuts in a direction worth noting: for a
pseudogene caller, resilience and detection pull in opposite directions, and an annotation whose
feature counts barely move when genes are deliberately broken is one that is not reporting the
breakage. We extend Syberg-Olsen *et al.* (2022) by supplying the non-circular cross-tool comparison
their simulation module could not.

## 13 Limitations

**13.1 No independent biological ground truth, and none obtainable here.** Stage A of the simulator
ran; Stages B and C did not. A systematic audit (§ 11.3, `INDEPENDENT_VALIDATION_AUDIT.md`) found no
curated set of experimentally verified bacterial pseudogenes; the curated sets that exist are for
organisms outside the evaluable subset; and the best experimentally grounded cases — the *Shigella*
antivirulence genes — are **deletions**, which have no locus to classify. **No ranking of the four
configurations by accuracy is supported, and no claim is made that any called pseudogene has lost
function.** Our controls are an artificial-lesion detection control and an annotation-consensus
control, not independently verified truth.

**13.2 The dN/dS subset is not the corpus.** 92 of 523 genomes have a congeneric partner, and they are
not a random sample of the corpus. All dN/dS results are conditional on that subset.

**13.3 Orthologue assignment is verified only for the loci that reciprocate.** RBH retains 97.2 % of
concordant-intact and **67.5 %** of discordant loci, a 29.7-point imbalance; excluded discordant loci
have higher ω (genome-clustered +0.0479 over 34 eligible genomes, p = 3.9 × 10⁻⁵). The gap among
excluded loci is unstable and imprecise, so **the direction of the residual bias cannot be signed**.
Both statistical corrections failed. **17.21 %** of discordant loci have a near-tied forward runner-up
against 1.03 % of intact loci, and in all 1 166 such cases the runner-up is a strong sequence match to
the winner — consistent with enrichment in multi-copy families, though **association is not verified
paralogy** without independent orthology evidence, which this study lacks.

**13.4 Baseline contamination exists and its direction is untested.** 65–73 % of truly lesioned genes
are annotated intact, so real pseudogenes sit in the concordant-intact baseline. We do **not** claim
this attenuates the gap: the missed genes are not a random sample of pseudogenes but those whose lesion
leaves a long surviving ORF.

**13.5 Lesion stratification reaches 11.0 % of discordant loci**, DFAST only.

**13.6 The phylogenetic criterion T1 failed as written** and its permutation half is untestable on this
corpus (§ 9.7).

**13.7 Bakta's result depends on which of its two keys is counted.** Under its own strict definition the
caller-level ω result is not significant (p = 0.076); including `pseudo=True` gives p = 0.036. Both are
reported and neither was chosen on its p-value.

## 14 Conclusions

1. Bacterial pseudogene counts vary **19- to 77-fold** between annotation configurations on a fixed
   corpus, with one configuration reporting none at all in a quarter of genomes.
2. One widely used configuration is **not independent** of another, and their agreement is structural.
3. On the canonical genome-reduction claim, **no configuration finds the predicted association**.
   Compared directly on the same genomes, the configurations **differ significantly on raw counts**
   (2 of 3 independent pairs) but are **not statistically distinguishable on normalised metrics in these tests**, even where their
   significance verdicts differ. This absence of significant paired contrasts on normalised metrics does not demonstrate equivalence or eliminate between-pipeline differences in measured pseudogene burden.
4. Loci where configurations disagree carry **elevated dN/dS at genus level**, robust to analytical
   choice and reciprocal-search sensitivity — but **synthetic frameshifts can mechanically inflate ω, with the magnitude uncertain in this small control arm**. Thus dN/dS cannot adjudicate the disagreement. Independently of dN/dS, those loci
   carry alignment-geometry evidence of disruption far more often than loci all four configurations
   call intact (**33.3 % against 3.5 %**). This is **sequence** evidence; whether the genes have lost
   function is not something these data can establish.
5. Pseudogene counts should be reported with the pipeline, version and assertion key that produced
   them, and should not be compared across pipelines without a locus-level reconciliation.

## 15 Data availability

All 523 assemblies are public RefSeq records; accessions are in `results/corpus523.tsv`. Derived
per-locus tables, per-genome summaries and all result JSONs are in the repository under `results/` and
`benchmark/results/`, indexed by `results/result_registry.tsv`.
**Archive DOI:** pending. A persistent archive of the derived tables is required before journal submission.

## 16 Code availability

All analysis code is in the repository, including the frozen pre-registered plans
(`benchmark/SIMULATION_PLAN.md`, `FRAME2_PLAN.md`, `PHYLO_PLAN.md`, `RBH_STRATIFIED_PLAN.md`,
`FORWARD_AMBIGUITY_PLAN.md`, `RBH_EXPANSION_PLAN.md`, `LESION_HARMONISATION_PLAN.md`), each asserting
its own SHA-256 at runtime, and the verification suite (§ 17).
**Archive DOI:** pending. The public repository link and licensing should be verified at submission.

## 17 Reproducibility and verification

Original result files are supplied under `results/`, `benchmark/results/`, and `supplementary/`, with the five source figures under `figures/`. The release preserves frozen plans with their original hashes and records failed quality gates without recategorisation. Unit tests and figure-generation scripts are included, but full genome-wide reruns require public RefSeq assemblies, external executables, and local intermediate data that are not all bundled. See `REPRODUCIBILITY.md`.

## 18 Author declarations (items requiring author confirmation)

**Authors and affiliations.** Gurgen Karapetyan¹˒² and Astghik Pepoyan¹. ¹ Armenian National Agrarian University, Yerevan, Armenia; ² Scientific Center of Zoology and Hydroecology, Yerevan, Armenia.

**Corresponding author.** Gurgen Karapetyan. Email and ORCID to be confirmed before submission.

**Funding.** Grant 25AA-4A009 is the only grant attributed to this study; the authors should confirm the exact funding-agency wording and the grant's relevance to this work.

**Author contributions, competing interests and acknowledgements.** Require explicit confirmation by both authors.

**Ethics and consent.** The analysis uses public bacterial genome assemblies and includes no newly recruited human or animal participants; the authors should confirm the applicable journal declaration.

**Data and code.** A persistent release DOI and final archive URL have not yet been assigned.

## 19 References

Benjamini, Y. and Hochberg, Y. (1995) Controlling the false discovery rate. *Journal of the Royal
Statistical Society B* **57**, 289–300. doi:10.1111/j.2517-6161.1995.tb02031.x

Buchfink, B., Xie, C. and Huson, D.H. (2015) Fast and sensitive protein alignment using DIAMOND.
*Nature Methods* **12**, 59–60. doi:10.1038/nmeth.3176

Cole, S.T., Eiglmeier, K., Parkhill, J., James, K.D. *et al.* (2001) Massive gene decay in the leprosy
bacillus. *Nature* **409**, 1007–1011. doi:10.1038/35059006

Cooley, N.P. and Wright, E.S. (2024) Many purported pseudogenes in bacterial genomes are bona fide
genes. *BMC Genomics* **25**, 365. doi:10.1186/s12864-024-10137-0 *(PMID 38622536, PMC11017572)*

Goodhead, I. and Darby, A.C. (2015) Taking the pseudo out of pseudogenes. *Current Opinion in
Microbiology* **23**, 102–109. doi:10.1016/j.mib.2014.11.012

Jundzill, M., Hölzer, M., Mangul, S., Marquet, M., Ehricht, R., Lohde, M., Spott, R., Makarewicz, O.,
Pletz, M.W. and Brandt, C. (2026) Large-scale benchmarking of prokaryotic annotation tools across
thousands of species. *Genome Biology* **27**, 284. doi:10.1186/s13059-026-04262-0
*(PMID 42706541, PMC13548610)*

Kuo, C.-H. and Ochman, H. (2010) The extinction dynamics of bacterial pseudogenes. *PLoS Genetics*
**6**, e1001050. doi:10.1371/journal.pgen.1001050

Lerat, E. and Ochman, H. (2004) Ψ-Φ: exploring the outer limits of bacterial pseudogenes. *Genome
Research* **14**, 2273–2278. doi:10.1101/gr.2925604

Madin, J.S., Nielsen, D.A., Brbic, M., Corkrey, R. *et al.* (2020) A synthesis of bacterial and archaeal
phenotypic trait data. *Scientific Data* **7**, 170. doi:10.1038/s41597-020-0497-4

McCutcheon, J.P. and Moran, N.A. (2012) Extreme genome reduction in symbiotic bacteria. *Nature Reviews
Microbiology* **10**, 13–26. doi:10.1038/nrmicro2670

Moran, N.A. and Plague, G.R. (2004) Genomic changes following host restriction in bacteria. *Current
Opinion in Genetics & Development* **14**, 627–633. doi:10.1016/j.gde.2004.09.003

Nei, M. and Gojobori, T. (1986) Simple methods for estimating the numbers of synonymous and
nonsynonymous nucleotide substitutions. *Molecular Biology and Evolution* **3**, 418–426.
doi:10.1093/oxfordjournals.molbev.a040410

Schwengers, O., Jelonek, L., Dieckmann, M.A., Beyvers, S. *et al.* (2021) Bakta: rapid and standardized
annotation of bacterial genomes. *Microbial Genomics* **7**, 000685. doi:10.1099/mgen.0.000685

Syberg-Olsen, M.J., Garber, A.I., Keeling, P.J., McCutcheon, J.P. and Husnik, F. (2022) Pseudofinder:
detection of pseudogenes in prokaryotic genomes. *Molecular Biology and Evolution* **39**, msac153.
doi:10.1093/molbev/msac153

Tanizawa, Y., Fujisawa, T. and Nakamura, Y. (2018) DFAST: a flexible prokaryotic genome annotation
pipeline. *Bioinformatics* **34**, 1037–1039. doi:10.1093/bioinformatics/btx713

Tatusova, T., DiCuccio, M., Badretdin, A., Chetvernin, V. *et al.* (2016) NCBI prokaryotic genome
annotation pipeline. *Nucleic Acids Research* **44**, 6614–6624. doi:10.1093/nar/gkw569

*(Full author lists for works cited with "et al." are in `MANUSCRIPT.md`; they are preserved and will be
expanded to the target journal's style at submission.)*

## 20 Supplementary materials

| item | content |
|---|---|
| S1 | Per-genome pseudogene counts and fractions, four configurations, 523 genomes |
| S2 | Per-genus RBH outcome, 28 genera (`results/rbh3_all/by_genus.tsv`) |
| S3 | The 48 literature re-tests with ρ, partial ρ, p and verdict |
| S4 | Seed sensitivity, 40 seeds, four configurations |
| S5 | Overlap-threshold sweep, with its conditional nature stated |
| S6 | FRAME2 arms and the F4 gate failure |
| S7 | Lesion classification and controls |
| S8 | Full gate register, including every failed gate |
| S9 | Provenance: result registry, plan hashes, verification output |

---

### Figures

| figure | content | status |
|---|---|---|
| **1** | study design and annotation workflow | built, `figures/figure1_study_design.*` |
| **2** | caller-dependent abundance and disagreement | built, `figures/figure2_caller_disagreement.*` |
| **3** | Metric-dependent genome-size association and direct caller comparison | built, `figures/figure3_downstream_consequence.*` |
| **4** | dN/dS: frameshift artefact, reference and orthologue sensitivity | built, `figures/figure4_dnds_sensitivity.*` |
| **5** | Sequence-based lesion evidence on extended intervals | built, `figures/figure5_lesion_evidence.*` |

