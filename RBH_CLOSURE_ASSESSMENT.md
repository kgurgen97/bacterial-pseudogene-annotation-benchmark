# RBH closure assessment

Written at commit `fa49b94`. Sources: `MANUSCRIPT_CANDIDATE_RBH_V2.md`,
`RBH_ORTHOLOGY_VALIDATION.md`, `results/rbh_audited_all.json`,
`results/rbh_stratified_all.json`, `results/forward_ambiguity_all.json`,
`results/review_fixes/review_fixes.json`, `results/manuscript_verify/`.

---

## 1 Verdict

> **Ready to incorporate into the full manuscript — as an orthologue-assignment sensitivity analysis,
> with its failed gates intact. Not ready, and never will be from these data, as biological
> validation of G3.**

Three conditions make it incorporable, and all three now hold:

1. **Every number traces to committed evidence** and is machine-re-derived: 19/19 numerical claims,
   12/12 semantic checks, 13/13 negative controls caught, plus the three pre-existing report
   verifiers. Verification code itself carries controls because **three defects were found in it**
   during review (a missing `hashlib` import, unfolded unicode superscripts, and whole-word quote
   markers) — each one a case of the checker being wrong rather than the analysis.
2. **The claim is scoped to what the design supports.** The candidate states in its masthead, § A.4 and
   § D.1 that this is a sensitivity analysis; it says the paralogue-pairing concern is **reduced, not
   eliminated**; and the four older places that said "tested and rejected" are corrected with the
   original wording recorded.
3. **No failed gate is hidden.** Scientific gate FAIL overall; retention 29.7 pp FAIL; matching
   balance FAIL under the frozen criterion; IPW same-genome-set recovery −0.0087 FAIL; IPW improving
   11 of 30 genomes FAIL. The misleading script-level IPW PASS is named and excluded.

## 2 What the analysis established

| finding | evidence | unit |
|---|---|---|
| Orthologue assignment reciprocates for most loci | 93.2 % modified RBH; 94.8 % standard where applicable | locus |
| The gap persists among RBH-confirmed loci | **+0.0291**, genus CI **[+0.0226, +0.0339]** | genus-clustered |
| Direction is consistent across genera | **24 of 28 positive**, genus sign p **9.0 × 10⁻⁵**, Wilcoxon 2.1 × 10⁻⁷ | genus |
| Filtering does not detectably change the per-genome gap | **median paired difference +0.000000** over the same 90 genomes (44 lower / 35 higher / 11 unchanged, sign p 0.368, Wilcoxon p 0.13); genus-level paired −0.000318 (17/28, p 0.108); mean paired −0.003779; difference of marginal medians −0.002551 | genome (paired), genus |
| Ambiguous loci do not drive the result | dropping all of them: **+0.00146**, no genome flips | genome |
| Exclusion is non-random and its ω direction is now established | excluded discordant ω higher: genome **+0.0479** (34/92 units, p 3.9 × 10⁻⁵); genus **+0.0283** (19/28, p 0.0044) | genome and genus |
| Forward ambiguity is concentrated and family-associated | 17.21 % vs 1.03 %; all 1 166 near-ties have a homologous runner-up | locus (descriptive) |
| CDS absence does **not** explain modified-RBH attrition | CDS-less loci confirmed **more**: 73.8 % vs 61.8 %, χ² p 7.5 × 10⁻¹⁹ | locus |

The strongest single sentence the data support: **replacing a one-way best hit with a reciprocal best
hit plus synteny, across 92 genomes in 28 genera, leaves the dN/dS gap at caller-disputed loci
unchanged at the median on a paired within-genome comparison (0.0000, sign p 0.37; four estimands span
0.000 to −0.004), positive in 24 of 28 genera, with a genus-clustered interval excluding zero.**

## 3 What failed, and why that is part of the result

| failure | what it means |
|---|---|
| **Retention imbalance 29.7 pp** | RBH keeps 97.2 % of intact loci and 67.5 % of discordant ones. The filter is informative about the thing being tested |
| **IPW cannot correct it per genome** | fittable in 30 of 92 genomes; same-set bias −0.0087, outside the frozen ±0.005; median error *worse* than uncorrected (0.0072 vs 0.0062). Only the pooled statement survives |
| **Matching never balances at scale** | all ten grid cells imbalanced, max \|SMD\| 9.0. The criterion was declared as a **maximum over units**, which is not comparable between 12 and 92 genomes — median \|SMD\| is 0.070/0.032/0.024. Both reported; the matched arms are cited as support for nothing |
| **E3 unstable** | +0.0423 per genome (p 0.215) but +0.0053 per genus with 11/18 positive; 32/92 genomes eligible. The residual bias from filtering therefore **cannot be signed** |
| **Clustered ω test is eligibility-selected** | only genomes that exclude a lot can be tested, which is why the clustered difference (+0.0479) is three times the pooled one (+0.0161) |

## 4 Questions outside this analysis's scope

**Not addressed by any RBH result, and not addressable by more RBH:**

1. **Whether elevated ω at these loci is biological at all.** ω is an association and plausibility
   measure. The FRAME2 control shows a frameshift **alone** raises measured ω by +0.0315 — the size of
   the whole published effect — so an elevated value can come from the lesion's effect on the
   alignment. RBH tests orthologue *assignment*, not alignment behaviour, and leaves this untouched.
2. **Caller accuracy.** Stage A of the independent simulator ran; **Stage B and Stage C never did**,
   being forbidden under the host disk guard. No claim about which caller is right is supported at
   corpus scale, and the pre-registered gates B1 and F4 **failed and remain failed**.
3. **True orthology.** Both search directions are DIAMOND best hits over annotation-derived proteomes,
   and the partner reference is RefSeq-derived. There is **no independent orthology ground truth** in
   this study, so "RBH-confirmed" means "reciprocates under these searches", not "is an orthologue".
4. **The lesion-class composition of the gap beyond DFAST.** The internal-stop stratification covers
   **539 of 4 886 discordant loci (11.0 %)** and only DFAST's calls; Pseudofinder's 84.1 % share is
   unclassified. Whether the gap survives the frameshift artefact **across callers** is unknown.
5. **Baseline contamination's direction.** 65–73 % of truly lesioned genes are annotated intact, so
   real pseudogenes sit in the baseline. The naive attenuation argument is withdrawn as unsound, and
   the direction is untested.
6. **What the excluded loci would have contributed.** E3 is too imprecise to say, at either scale.
7. **Whether genuine paralogy can ever be resolved here.** 17.21 % of discordant loci have a
   near-tied, homologous runner-up. That is association consistent with multi-copy families; assigning
   those loci would need evidence this study does not have, and **no reciprocity rule, threshold or
   aligner removes real paralogy.**

## 5 Recommended handling in the full manuscript

- Place it as a **Methods subsection plus a short Results paragraph**, with Supplementary Table S-RBH
  (28 genera). Not a headline result.
- Quote **genus-level** figures beside the genus-clustered interval: +0.0291, [+0.0226, +0.0339],
  24/28 genera, sign p 9.0 × 10⁻⁵. Genome-level p-values (7.8 × 10⁻¹⁴) are **not** genus-adjusted and
  should be labelled if used at all.
- State the retention imbalance **in the same breath** as the confirmation rate. 93.2 % alone is
  misleading.
- Keep **E1/E2/E3** labelled and never present E2 as a corrected E1.
- Carry § D.1–D.8 of the candidate into Limitations essentially verbatim; they are the honest cost of
  the result.
- Do **not** write that paralogue misassignment is excluded.

## 6 One-line closure

**The orthologue-assignment concern is as closed as these data can close it; the ω-artefact,
caller-accuracy and true-orthology questions are untouched and need the simulator stages that the disk
guard forbids.**
