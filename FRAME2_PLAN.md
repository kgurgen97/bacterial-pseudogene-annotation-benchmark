# FRAME2_PLAN v1 — frozen 2026-10-31, before any ω is computed

**This file is frozen.** Its SHA-256 is asserted at runtime by `analysis_frame2.py`, which refuses to
run if it changed. It is never edited because of a result. `FRAME_PLAN.md` is untouched and the
original frame analysis stays void.

## 1. Why a new plan at all

`FRAME_PLAN.md`'s F4 was a gate: *"DFAST's `pre` gap stays below +0.010, as its full gap is +0.002 at
p = 0.52"*, failing if *"DFAST's `pre` gap becomes large and significant, in which case the split is
creating differences and this whole analysis is void."* After the classifier correction DFAST's full
gap is +0.023 at p = 1.9 × 10⁻⁹ and its pre gap is +0.018. **The gate fires.**

The gate used DFAST as a negative control *because DFAST had no signal*, and that absence was an
artefact of a substring parser. A gate resting on a defective number cannot survive the number being
fixed. The threshold is not moved and F4 is not reinterpreted. What is needed is a negative control
that does not depend on any caller's behaviour.

## 2. The question, unchanged

Does the procedure used to align a frameshifted gene to an orthologue **mechanically inflate**
estimated dN/dS, independently of any real relaxation of selection?

This matters because § 3.7 and G3 rest on elevated ω at pseudogene calls, and a frameshift makes the
downstream codons translate out of register. If the aligner renders that as amino-acid mismatches,
dN rises for a purely mechanical reason.

## 3. The design, and why its control is exact

The simulated assemblies already built for Stage A supply something the natural corpus cannot: **the
same gene, with and without a known frameshift, against the same orthologue.**

For a gene *g* with a simulated lesion, compute

    Δω(g) = ω(g in the modified assembly) − ω(g in the parental assembly)

through the identical pipeline — DIAMOND `blastx` against the partner's proteins, the project's
`codon_pairs`, the project's NG86. Everything is held constant except the lesion. Δω is therefore the
inflation attributable to the lesion and the alignment of it, with no confounding by gene identity,
orthologue choice, divergence, or caller behaviour.

**Experimental unit: the gene.** Tests are paired.

## 4. Arms, fixed now

| arm | lesions | role |
|---|---|---|
| **F** frameshift | `M2_fs_ins1`, `M3_fs_del1`, `M4_fs_indel2`, high tier only | the primary test |
| **S** premature stop | `M1_stop_sub`, high tier only | comparator: a stop does **not** shift the frame |
| **N1** untouched | eligible genes carrying no lesion in that replicate | exactness control |
| **N2** synonymous | `C1_synonymous` | biological-null control: the protein is byte-identical |

## 5. Predictions and the bars that decide them, fixed now

| # | prediction | bar |
|---|---|---|
| **N1** | **gate.** An untouched gene's ω is identical in both assemblies | every untouched gene must give Δω exactly 0. A single non-zero Δω means the pipeline is not deterministic across assemblies and **the whole analysis is void** |
| **N2** | a synonymous change does not inflate ω | median \|Δω\| < 0.010 |
| **P1** | **primary.** A frameshift inflates ω mechanically | median Δω over arm F, sign test over genes, bootstrap 95 % CI |
| **P2** | the inflation is specific to frame change | median Δω in arm F exceeds arm S |

**P1's interpretation is fixed before the data are seen**, so the answer cannot be chosen afterwards:

- median Δω ≥ **+0.030** — as large as the gaps § 3.7 and G3 report. The artefact could then account
  for those gaps and both are **seriously threatened**.
- median Δω < **+0.010** — too small to account for them. The artefact exists but does not explain
  the findings.
- between +0.010 and +0.030 — **INCONCLUSIVE**, reported as such, with the CI.

If N1 fails, nothing else is reported. If arm F yields fewer than 5 non-tied genes, P1 is **NOT
TESTABLE** and that is the result.

## 6. Known limitations, stated before running

- **One genome.** Only `GCF_000006925.2` among the two simulated genomes has a dN/dS partner
  (`GCF_000012005.1`). The other has none, so this is a single-genome, single-partner analysis.
- **Small n.** Roughly 15–25 high-tier frameshift genes. The exact pairing buys power — the only
  difference between the two measurements is one or two bases — but this can detect a large
  inflation, not resolve a small one. The bars in § 5 are set accordingly.
- **A simulated frameshift is recent and single.** Natural pseudogenes are old and multiply damaged,
  so the mechanical inflation measured here is a *lower* bound on what an old, heavily damaged gene
  might suffer.
- This measures the **pipeline's** behaviour, not biology. It cannot say whether any natural
  pseudogene call is right.
- Genes whose ω is undefined in either assembly (no orthologue hit, dS ≤ 0, below the codon floor)
  drop out, and the number dropped is reported per arm.

## 7. Outputs

`FRAME2_RESULTS.md` · `results/frame2.json` · `results/frame2_genes.tsv` with one row per gene and
both ω values · `benchmark/tests/test_frame2.py`.

Post hoc additions, if any, are labelled post hoc and kept out of § 5's verdicts.
