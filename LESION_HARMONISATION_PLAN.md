# Frozen plan: caller-independent lesion evidence on extended intervals

**Frozen before any code for it exists.** `analysis_lesion_extended.py` will assert this file's
SHA-256 at runtime. A frozen plan is never edited because of a result; if a result demands a different
design, the remedy is a new, separately labelled plan.

---

## 1 Hypothesis

Caller disagreement at discordant loci can be partly decomposed into **definitional** differences
(which lesions a pipeline counts as a pseudogene), **detection** differences (whether it finds the
lesion at all) and **boundary** differences (where it places the feature) — and that decomposition can
be driven by sequence evidence collected **independently of any caller**, provided the evidence is
gathered on **extended** genomic intervals rather than annotation-bounded ones.

## 2 The exact scientific question

Does extending each locus's interval beyond its annotated boundaries recover the lesion evidence that
annotation-bounded analysis misses — specifically, does it raise frameshift recall from the **4.4 %**
(18 of 413) achieved by `analysis_lesion_from_sequence.py` to a level that supports classifying
discordant loci by lesion type?

**Why it failed before, stated as the mechanism this design targets.** A caller that detects a
frameshift typically reports the *truncated* feature, so the post-lesion fragment lies outside the
interval, and a frameshift signature — two alignment segments in different reading frames — cannot be
observed. The interval, not the detector, was the limiting factor.

## 3 Expected contribution

If it works: the study gains a **caller-independent lesion classification** with usable coverage, which
converts "the pipelines disagree" into "they disagree because they define X differently and detect Y
differently". That is the strongest contribution available from existing data.

If it fails: the study gains a documented, quantified statement that **lesion truth cannot be recovered
from annotation-derived intervals even when extended**, which bounds every sequence-level claim in the
paper and in the published literature that relies on annotation intervals.

Either outcome is reportable. The design is not contingent on the sign.

## 4 Inputs

| input | source | role |
|---|---|---|
| caller-neutral loci, 12-genome pilot | `results/reference_sensitivity/per_locus.tsv` (sha256 `9ec66048bcd30fac…`) | the loci to classify |
| partner proteomes | `~/research/datasets/n2-refseq-reference/*.faa` (read-only, sentinel-guarded) | orthologue reference |
| contig sequences | `~/research/datasets/n2-pseudofinder/*/pf_contigs.fasta` | for interval extension |
| DFAST per-call lesion class | `~/research/datasets/n2-dfast/*/pseudogene_summary.tsv` | **external validation only**, never an input to classification |

No new downloads. No re-annotation. ω is not used and not recomputed.

## 5 Method

**Interval extension.** For a locus on contig *c* spanning [s, e] with orthologue nucleotide length
*L*, the extended span is [s − m, e + m] clipped to the contig, with
**m = max(300, ceil(0.5 × L))**. The margin is tied to the orthologue so that a post-lesion fragment
of plausible size is captured, and floored at 300 bp so that short genes are not left unextended.

**Evidence collection.** `diamond blastx --max-target-seqs 1 --max-hsps 0 --evalue 1e-5` of the
extended span against the partner proteome. All HSPs are retained.

**Classification**, applied in this order, using HSPs of ≥ 20 aa against the single best subject:

| class | rule |
|---|---|
| `frameshift` | ≥ 2 distinct reading frames `(qstart−1) mod 3`, and the HSPs in different frames cover **non-overlapping** segments of the subject (≥ 20 aa each) |
| `internal_stop` | all strong HSPs in one frame, and a stop codon occurs in that frame before 90 % of the subject's length is covered |
| `truncation` | one frame, no premature stop, and total subject coverage < 50 % |
| `fragmentation` | ≥ 2 non-overlapping HSPs in the **same** frame separated by ≥ 30 subject residues |
| `intact` | one frame, no premature stop, subject coverage ≥ 90 % |
| `ambiguous` | anything else |

## 6 Primary outcome

**Frameshift recall** against DFAST's published `frameshift` class, on loci where DFAST assigns one.
Secondary: internal-stop recall, decided fraction, and the false-frameshift rate.

## 7 Controls

| control | construction | what it bounds |
|---|---|---|
| **positive, synthetic** | single 1-bp deletion introduced at 40 % of the length of 100 randomly chosen concordant-intact loci, classified through the identical path | does the pipeline detect a lesion it is known to contain |
| **negative, biological** | 400 concordant-intact loci per genome, unmodified — all four callers call these intact | the **false-lesion rate**. The previous analysis had no such control, so its specificity was unmeasured |
| **negative, synthetic** | the same 100 loci unmodified | frame and boundary noise in the absence of any lesion |
| **external** | DFAST's own class, used only after classification | agreement, not truth |

The negative biological control is the one that decides whether the method is usable: a classifier that
labels intact genes as lesioned has no value regardless of its recall.

## 8 Statistical estimand and units

Classification is per **locus**. Recall, precision and the false-lesion rate are proportions over loci,
reported with Wilson 95 % intervals. Any statement comparing callers is summarised per **genome** and
per **genus**, never pooled across loci, because loci within a genome are not independent.

## 9 Biases and confounders, named in advance

1. **Extension imports neighbouring genes.** A longer span may align to a different gene and manufacture
   an apparent multi-frame signature. Mitigation: all HSPs must be against the **same best subject**,
   and the frameshift rule requires non-overlapping subject segments.
2. **Overlapping loci.** Extended spans of adjacent loci will overlap; each locus is classified on its
   own span and no locus is informed by another's.
3. **DFAST's class is itself a caller output.** It is used for external agreement only; it cannot
   establish truth, and disagreement with it is not automatically an error.
4. **The reference is RefSeq-derived**, so independence is from Pseudofinder, not from annotation.
5. **Short genes** get a proportionally larger extension, which may make them more likely to be
   classified `fragmentation`. Reported as a stratified diagnostic by locus length.

## 10 Success criteria — declared now

| gate | criterion |
|---|---|
| **G1, necessary** | frameshift recall **≥ 40 %** against DFAST frameshift calls |
| **G2** | false-frameshift rate on DFAST internal-stop loci **≤ 10 %** |
| **G3, necessary** | false-lesion rate on the biological negative control **≤ 10 %** |
| **G4** | decided fraction (not `ambiguous`) **≥ 60 %** of discordant loci |
| **G5** | internal-stop recall **≥ 70 %** — no regression from the 82.4 % (280/340) already achieved |
| **G6** | synthetic positive control: ≥ 80 % of the 100 seeded 1-bp deletions classified `frameshift` |

## 11 Failure criteria — also declared now

- **G1 or G3 fails** → the method is **not usable**. Report the measured recall and false-lesion rate,
  state that lesion truth cannot be recovered from annotation-derived intervals even when extended, and
  **do not decompose caller disagreement by lesion class**. Do not widen the margin and retry more than
  once, and if a second margin is tried it is reported as a second, separately labelled arm.
- **G6 fails** → the implementation is suspect before any interpretation; debug the classifier, not the
  biology.
- Any gate may fail without invalidating the rest of the manuscript: this is an addition, not a
  dependency.

## 12 Resource estimate and stop conditions

Pilot: the same 12 genomes, ~1 871 discordant plus 400 baseline loci each. One DIAMOND `blastx` per
genome over extended spans; measured footprint of the comparable earlier run was a few MB per genome
and under a minute of compute. Expected total **< 100 MB transient, < 30 minutes**.

Stop conditions: host C: below **10 GB** free at start or **8 GB** at any checkpoint, read from
`/mnt/c` and never `df /`; any protected file changing digest; any gate in § 10 failing, which stops
expansion rather than the run.

Expansion to 92 genomes only if **G1, G3, G5 and G6 all pass**, and the expansion is a separate
labelled run writing to its own directory.

## 13 Scope

This plan authorises the 12-genome pilot and nothing else. It does not authorise re-annotation, new
downloads, Stage B/C simulation, or any change to ω.
