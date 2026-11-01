# A caller-neutral locus set, and G3 at last

**Written 2026-10-30, before `analysis_locus.py` existed and before any outcome was computed.**
Nothing below may be changed because of a result. If a decision here turns out to be wrong, the fix is
a new, separately labelled analysis, not an edit to this file.

`DNDS_PLAN.md` and `FRAME_PLAN.md` are untouched.

## 1. Why this exists

`DNDS_PLAN.md` registered G3: among the sequences **one caller calls a pseudogene and another calls
intact**, does dN/dS differ from the all-intact baseline? It is the only registered prediction that
adjudicates between callers gene by gene, and it is the one the project could not run. The reason
given, correctly, was that the four annotations share no gene identity: each emits its own intervals
under its own contig names.

That reason was sufficient to stop the original analysis and it is not sufficient forever. The
callers annotate the same assemblies, and § 3.7 already reconciles contig names by length. Once names
are reconciled the intervals live in one coordinate system, and overlap defines a locus. This builds
that locus set and runs G3 as `DNDS_PLAN.md` wrote it.

## 2. How a locus is defined

Fixed now, before any locus is built.

**Coordinates.** Contig names are reconciled exactly as `analysis_dnds.py` does: by length where
lengths are unique, falling back to order, with the fall-back counted and reported.

**A locus** is a connected component of intervals on the same reconciled contig and the same strand,
where two intervals are linked when they overlap by at least **50 % of each** (reciprocal overlap).
Connected components are used rather than pairwise clusters because a locus may be called with
slightly different boundaries by each tool; the cost is that a chain of partial overlaps can merge two
genes, and the number of loci containing more than one interval from the same caller is reported as
the measure of that.

**A caller's verdict at a locus** is `pseudo` if one of its pseudogene intervals is in the component,
`intact` if one of its intact intervals is, `absent` if neither. A locus where one caller contributes
both a pseudogene and an intact interval is **ambiguous** and is excluded, with the count reported.

**A discordant locus** has at least one caller saying `pseudo` and at least one saying `intact`.
**A concordant-intact locus** has at least two callers saying `intact` and none saying `pseudo`; these
form the baseline.

**dN/dS at a locus** is computed on the longest interval in the component, by the same route as
§ 2.5: aligned to the partner genome's best DIAMOND hit at 50 % identity and 50 % coverage, in the
partner's reading frame, Nei-Gojobori with the Jukes-Cantor correction, 30 gap-free codons minimum.

**Scope.** The same 99 target/partner pairs as § 3.7, so the locus analysis and the sequence analysis
speak about the same genomes.

## 3. Pre-committed predictions, with failure criteria

| # | prediction | confirmed if | **failed if** |
|---|---|---|---|
| L1 | **gate:** the locus set is usable | at least 500 discordant loci across the corpus, and fewer than 20 % of intervals join no locus | fewer than 100 discordant loci, or more than 40 % of intervals unjoined |
| G3 | **the registered prediction, from `DNDS_PLAN.md` § 3:** where callers disagree, selection takes a side | median dN/dS at discordant loci differs from the concordant-intact baseline at p < 0.05, by a genome-clustered test | it does not |
| L2 | **the adjudication:** where Pseudofinder and PGAP disagree, one of them is closer to selection | at loci where Pseudofinder says pseudo and PGAP says intact, median dN/dS differs from loci where PGAP says pseudo and Pseudofinder says intact, p < 0.05 | the two sets do not differ |
| L3 | **a prediction of weakness:** overlap is a crude gene identity | more than 5 % of loci contain two or more intervals from the same caller, which would mean the components are merging genes | fewer than 5 %, in which case the merging worry was unfounded |

**L1 is a gate.** If it fails, the locus set is not fit to carry G3 and nothing below it may be read.

**All tests are clustered by genome**, one value per genome where the design allows, because the
pooled per-sequence tests in § 3.7 were shown to overstate evidence by orders of magnitude.

## 4. What this cannot settle

- Overlap is not orthology. Two tools calling overlapping intervals may still mean different genes,
  and L3 measures how often the components merge rather than fixing it.
- dN/dS stays one-sided: high supports a pseudogene call, low does not refute one.
- A caller that never annotates a region contributes `absent`, not evidence. Bakta will be `absent`
  often, so the discordant set will under-represent it.
- The partner choice and the estimator are inherited from `DNDS_PLAN.md` with the limits stated there.
- L2 compares two sets of loci, not two sets of genes matched to each other, so it says which caller's
  disagreements look more like relaxed selection, not which caller is right about any one gene.

## 5. Deliverables

`analysis_locus.py`, asserting this file's SHA-256 before computing, writing `results/loci.tsv` (one
row per locus: genome, contig, span, each caller's verdict, dN/dS, codons) and `results/locus.json`
(counts, each prediction's verdict). Outcome in `LOCUS_RESULTS.md`, reporting every verdict including
failures.
