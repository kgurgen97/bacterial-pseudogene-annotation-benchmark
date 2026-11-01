# A ground truth: does selection agree with the caller?

**Written 2026-10-28, before `analysis_dnds.py` existed and before any outcome was computed.** Nothing
below may be changed because of a result. If a decision here turns out to be wrong, the fix is a new,
separately labelled analysis, not an edit to this file.

## 1. The question, and why the project cannot answer it yet

Everything in this project measures disagreement. The four callers span 77-fold in median pseudogene
fraction, they split on the sign of a published claim, and the manuscript says plainly that none of this
identifies which caller is right. A reviewer will ask exactly that, and the honest answer today is that
we do not know.

A pseudogene is a sequence no longer maintained by purifying selection. That leaves a signature
independent of any annotation pipeline: the ratio of non-synonymous to synonymous substitution, dN/dS,
which sits well below 1 for a gene under purifying selection and drifts toward 1 once selection relaxes.
Comparing each caller's calls against that signature asks whether what it labels a pseudogene behaves
like one.

**Why this is legitimate here, which is not obvious.** Pseudofinder can use dN/dS as one of its own
criteria, which would make this circular. It does so only when given `--reference`, invoking PAML. Our
runs passed `annotate --genome --outprefix --database --diamond --threads` and nothing else, and the
run logs contain no mention of dN/dS. Our Pseudofinder calls were therefore made without it, so dN/dS is
independent of all four callers. This was checked before the plan was written.

## 2. What is compared

**Targets.** Genomes with a congeneric partner of a different species inside the corpus. Genus is the
first word of the organism name, except that `Candidatus X` yields genus `X`, because grouping every
`Candidatus` together would pair unrelated endosymbionts. Where several partners qualify, the one with
the **lowest accession** is used: an arbitrary rule, fixed now so it cannot be tuned.

Only the 160 expansion genomes carry an organism name in `corpus523.tsv`, so the original 363 are out of
scope. A survey before writing this found roughly one hundred usable targets spread across all four
strata.

**Sequences.** For each caller and each target genome, two sets:

- its **pseudogene calls**, as genomic intervals from PGAP's `/pseudo` qualifiers, Bakta's GFF3,
  DFAST's GFF and Pseudofinder's `pf_pseudos.gff`;
- an equal-sized random sample of that caller's **intact** coding sequences, as the control.

**Orthologues.** Each sequence is matched against the partner genome's intact CDS by best DIAMOND hit,
requiring at least 50 % identity and 50 % coverage of the partner gene. Unmatched sequences are dropped
and counted.

**How dN/dS is computed on a sequence that may be frameshifted.** This is the methodological crux. A
pseudogene can be truncated or shifted, so it has no reading frame of its own, and translating it is
meaningless. The candidate is therefore aligned to its functional orthologue **in the orthologue's
frame**, and only gap-free aligned codons contribute. dN/dS is the Nei-Gojobori counting estimate over
those codons, implemented directly rather than through a library, and validated against a worked example
before use.

This has a consequence that must be stated in advance: a heavily degraded pseudogene leaves few gap-free
codons, so the estimate is noisy or impossible exactly where degradation is greatest. The measure is
therefore biased toward mildly degraded sequences, and a sequence with fewer than **30** usable codons is
dropped.

## 3. Pre-committed predictions, with failure criteria

| # | prediction | confirmed if | **failed if** |
|---|---|---|---|
| G1 | **sanity:** dN/dS carries signal here at all | for every caller, median dN/dS of its pseudogene calls exceeds that of its intact controls | it does not, for two or more callers |
| G2 | the callers differ in how well their calls match selection | the gap between pseudogene and intact median dN/dS differs by ≥ 0.15 between the best and worst caller | all four lie within 0.05 |
| G3 | **the question:** where callers disagree, selection takes a side | among sequences one caller calls pseudogene and another calls intact, median dN/dS differs from the all-intact baseline at p < 0.05 | no difference |
| G4 | **a prediction of weakness:** the signal is weakest where it is most needed | median dN/dS separation is smaller in S3, recent restriction, than in S1 | S3 separation is the largest of the four strata |

**G1 is a gate.** If it fails, dN/dS carries no usable signal on this corpus and nothing else here may be
interpreted. That outcome would be reported, and the manuscript would keep saying the project has no
ground truth.

**G4 is deliberately a prediction of failure of the method, not of the callers.** Recent pseudogenes have
had little time to accumulate substitutions, so S3 is where dN/dS should be least able to adjudicate.
Confirming it bounds what the whole analysis can claim.

## 4. What this cannot do

- dN/dS is a ground truth for *relaxed selection*, not for *annotation correctness*. A gene can be
  intact and weakly constrained, or recently broken and still look constrained. The measure is strong
  evidence in one direction only: high dN/dS supports a pseudogene call, low dN/dS does not refute it.
- Only the 160 expansion genomes are in scope; the original 363 have no organism name in the corpus
  table and are excluded.
- One partner per target, chosen by an arbitrary rule. A different partner would give different
  divergence and so different power.
- Nei-Gojobori is a counting estimator without correction for transition/transversion bias or codon
  usage. It is adequate for ranking callers and is not a substitute for a maximum-likelihood estimate.
- Sequences with fewer than 30 gap-free codons are dropped, which removes the most degraded cases.
- Bakta's floor applies here as everywhere: it makes few calls, so its sample will be small.

## 5. Deliverables

`analysis_dnds.py` → `results/dnds_pairs.tsv` (one row per sequence: caller, call, dN, dS, dN/dS,
codons used, stratum), `results/dnds.json` (per-caller summaries and each prediction's verdict).
Outcome in `DNDS_RESULTS.md`, reporting every verdict including failures.
