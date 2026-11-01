# Is the selection signal in § 3.7 an artefact of the reading frame?

**Written 2026-10-29, before `analysis_frame.py` existed and before any outcome was computed.**
Nothing below may be changed because of a result. If a decision here turns out to be wrong, the fix is
a new, separately labelled analysis, not an edit to this file.

This does not amend `DNDS_PLAN.md`, which stays frozen as it is. It is a new, post hoc analysis of a
question that plan did not ask.

## 1. The objection

A second review of the manuscript raised it and it is the strongest objection made so far.

`analysis_dnds.py` aligns each candidate to its functional orthologue **in the orthologue's frame**,
because a truncated or frameshifted pseudogene has no frame of its own. That is the right choice and it
has a consequence the plan did not consider. Past a frameshift, the candidate read in the partner's
frame yields codons that are close to random with respect to the partner's. Non-synonymous differences
accumulate mechanically, dN rises, and dN/dS rises with it, **with no selective history involved at
all**.

This bears directly on the published result. Frameshift detection is a core criterion for PGAP and for
Pseudofinder, the two callers with the largest gaps, +0.030 and +0.037. If the artefact dominates, the
claim that three callers label sequences that have genuinely stopped being maintained does not hold,
and § 3.7 would have to be withdrawn.

## 2. The test

The existing pipeline aligns the translated candidate against the translated partner and keeps every
column where neither side is a gap. That alignment already locates the damage: an indel in the
candidate shows up as a gap column, and everything 3' of it is suspect.

**The split.** For each sequence, record the index of the **first gap column in the protein alignment**.
Codon pairs collected before it are `pre`; those collected after it are `post`. Recompute dN/dS on `pre`
alone, with the same Nei-Gojobori estimator, the same Jukes-Cantor correction and the same 30-codon
floor, which is applied to the `pre` set rather than to the whole alignment.

Both classes are split the same way. Intact controls carry real indels too, so the comparison stays
like for like.

**What the split cannot do.** A frameshift does not always open a gap. The candidate is translated in
one frame across its whole length, so after a shift the protein diverges and the aligner may render
that as a run of mismatches instead. `pre` is therefore a proxy for "before the damage", not a
guarantee, and F2 below tests the same question without relying on gap detection.

## 3. Pre-committed predictions, with failure criteria

Full-alignment gaps, from `results/dnds.json`, are the reference: PGAP **+0.030**, Pseudofinder
**+0.037**, Bakta **+0.023**, DFAST **+0.002**.

| # | prediction | the artefact is NOT dominating if | **the artefact IS dominating if** |
|---|---|---|---|
| F1 | **the main test:** the gap survives on pre-damage codons | for PGAP and Pseudofinder, the `pre` gap keeps at least half the full gap, so ≥ +0.015 and ≥ +0.019 | either `pre` gap falls below a quarter of its full gap, so < +0.0075 or < +0.0094, or reverses sign |
| F2 | **no gap detection needed:** ω rises along the alignment faster in calls than in controls | the rise, median ω of the last 15 codon pairs minus the first 15, exceeds the controls' rise by < 0.02 | it exceeds it by ≥ 0.02 |
| F3 | **sanity:** the split captures damage | median codons before the first gap is higher in intact controls than in pseudogene calls | it is not, in which case F1 is uninterpretable and must be reported as such |
| F4 | **self-invalidation:** the split must not manufacture signal | DFAST's `pre` gap stays below +0.010, as its full gap is +0.002 at p = 0.52 | DFAST's `pre` gap becomes large and significant, in which case the split is creating differences and **this whole analysis is void** |

A result between a quarter and a half of the full gap in F1 is **INCONCLUSIVE**, not a pass, and is to
be reported that way.

**F4 is a gate.** If it fires, nothing else here may be interpreted.

## 4. What this cannot settle

- It tests whether the § 3.7 gap is inflated by frame damage. It does not make dN/dS two-sided, and
  nothing here measures recall.
- `pre` sets are shorter than full alignments, so the 30-codon floor removes more sequences and the
  surviving set is less degraded still. That pushes against finding a gap, which makes a surviving gap
  the more convincing outcome and a vanished one harder to attribute.
- Sequences whose first gap falls within the first 30 codons drop out entirely. Their number is
  reported per caller and per class, because losing them unevenly between classes would itself bias the
  comparison.
- One partner per target and the Nei-Gojobori estimator are inherited from `DNDS_PLAN.md`, with the
  limitations stated there.

## 5. Deliverables

`analysis_frame.py`, asserting this file's SHA-256 before computing, writing
`results/frame_pairs.tsv` (one row per sequence: caller, call, pre and post dN/dS, codon counts, first
gap index) and `results/frame.json` (per-caller `pre` gaps and each prediction's verdict). Outcome in
`FRAME_RESULTS.md`, reporting every verdict including failures.
