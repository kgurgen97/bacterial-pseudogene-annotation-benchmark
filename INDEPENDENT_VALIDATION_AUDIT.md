# Independent pseudogene ground truth: feasibility audit

**Question.** Can we independently demonstrate that our sequence-level lesion evidence corresponds to
real gene disruptions — using evidence that is *not* derived from PGAP, Pseudofinder, Bakta or DFAST?

**Verdict.** **No, not at benchmark grade, within this corpus and without new annotation.** The reasons
are specific and are recorded below rather than waved at. The manuscript is revised accordingly, and
its existing controls are relabelled to say exactly what they are.

---

## 1 Evidence categories, as the audit required

| | category | what would qualify | available to us? |
|---|---|---|---|
| **A** | experimentally demonstrated loss of function **with a verified causal genomic lesion** | complementation or knock-in restoring function, plus the lesion's coordinates | **partly — but see § 3** |
| **B** | independently sequence-validated frameshift or nonsense mutation with matched strain and assembly coordinates | resequencing or curation independent of the annotation pipelines | **no** |
| **C** | high-confidence intact gene controls | independently verified functional genes in the same assemblies | **no** |
| **D** | annotation-derived pseudogene labels | PGAP/Pseudofinder/Bakta/DFAST output, or a genome paper's own annotation | **abundant — and explicitly excluded as circular** |

## 2 What was searched

| source | result |
|---|---|
| curated databases of experimentally verified bacterial pseudogenes | **none found.** PHI-base curates virulence genes, not pseudogenes; PATRIC curates virulence factors; neither catalogues lesion-level pseudogene truth |
| Nuccio & Bäumler manually curated *Salmonella* pseudogene sets | exists and is used elsewhere as a truth set — **unusable here**, see § 3.1 |
| proteogenomic evidence on pseudogene loci (*Salmonella*, NAR 2022; *E. coli*) | real experimental evidence, but about **expression**, and for organisms we cannot evaluate (§ 3.1). It cuts *against* pseudogene calls rather than confirming them |
| *Shigella* antivirulence-gene literature (Maurelli "black holes"; Bliven & Maurelli review) | five documented genes — **cadA, nadA, nadB, speG, ompT** — with functional complementation evidence. **Mostly deletions**, see § 3.2 |
| our own corpus composition | contains the canonical degradation organisms: *Shigella* × 4, *M. leprae* Kyoto-2, *S.* Typhi Ty2, *Y. pestis* CO92, *Buchnera*, *Sodalis* × 4, *Coxiella* × 2 |

## 3 Why each promising route fails, specifically

### 3.1 The curated *Salmonella* set cannot be evaluated without new annotation

Our lesion classifier operates on **caller-neutral loci**, which exist only for genomes annotated by
all four configurations *and* possessing a congeneric partner for the orthologue alignment — the
**92-genome subset**. The corpus contains exactly **one** *Salmonella* (Ty2) and **no** *Escherichia*,
so Ty2 has no congeneric partner and is not in that subset. Evaluating it would require annotating
additional *Salmonella* genomes with four pipelines, which the execution rules forbid and which would
also reopen a completed workstream.

The same applies to the *E. coli* proteogenomics evidence: no *Escherichia* is in the corpus.

### 3.2 The best Category A evidence is of **deletions**, which have no locus to classify

This is the decisive point and it is a property of the biology, not of our method.

The *Shigella* antivirulence genes are the strongest experimentally grounded cases available:
reintroducing *cadA* into *S. flexneri* 2a **attenuates virulence**, and restoring functional *nadA*
and *nadB* likewise attenuates *S. flexneri* 5a. That is genuine Category A evidence of loss of
function.

But Maurelli's "black holes" are **large genomic deletions**. A deleted gene is *absent from the
assembly*: there is no interval, no locus, and no alignment geometry to classify. Our method detects
**lesions in genes that are still present** — frameshifts, premature stops, truncations. **It is
structurally incapable of being tested against a deletion**, and a deletion is not what any of the four
configurations reports as a pseudogene either.

So the best-documented experimental loss-of-function cases in the most relevant organism in our corpus
are, by their nature, outside the evaluable class.

### 3.3 Published pseudogene lists for these genomes are Category D

*M. leprae*'s ~1 100 pseudogenes, *Y. pestis* CO92's 149, *Shigella*'s few hundred — each comes from
that genome's own annotation effort. Using them as truth would evaluate one annotation against
another, which is the circularity this audit exists to avoid. They are excluded.

## 4 What our existing controls actually are — relabelled

The audit's central warning applies directly to us, and the manuscript has been corrected.

| control | what it was called | what it **is** | what it supports |
|---|---|---|---|
| 96 genes with a seeded 1 bp deletion | "synthetic positive control" | an **artificial-lesion detection control** | that the classifier detects a frameshift **it is known to contain** (96.9 %). **Not** real-world biological sensitivity |
| 4 489 loci all four configurations call intact | "biological negative control" | an **annotation-consensus control** | that the classifier rarely reports a lesion where four pipelines agree there is none (4.5 %). **Not** independently verified true negatives |
| DFAST's per-call lesion class | external validation | **another caller's output** | agreement, never accuracy |

None of these is independent biological ground truth, and the manuscript now says so in those words.

## 5 What this means for the manuscript's claims

**Unaffected — these are directly observed facts about annotation output, and need no external truth:**

- the 19.1-fold and 77.2-fold count spreads;
- Pseudofinder's structural dependence on PGAP (267 of 267);
- the 4 886 discordant loci and their distribution;
- that the configurations differ on raw counts and not on normalised metrics;
- that a frameshift alone inflates ω by +0.0315 — this is an **exact same-gene control with a lesion
  we introduced**, so its truth is construction, not annotation.

**Qualified — sequence evidence, not demonstrated biology:**

- the lesion-rate contrast (33.3 % vs 3.5 %) shows **alignment-geometry evidence of disruption** is far
  more common at disputed loci. It does **not** show those genes have lost function.

**Removed or never claimed:**

- any statement that a called pseudogene is non-functional;
- any ranking of the configurations by accuracy.

## 6 Stop rule applied

The audit's stop rule says: if genuinely independent ground truth cannot be found at sufficient
quality, document that clearly and revise the manuscript to focus on reproducible annotation
disagreement, separating verified sequence evidence from biological function.

**That is what has been done.** No bounded pilot is pre-registered, because there is no qualifying
evaluation set to pre-register against. Manufacturing one from annotation output would be precisely
the circularity the task forbids.

## 7 What a future study would need

Stated concretely so the gap is actionable rather than rhetorical:

1. **A strain with resequenced, coordinate-resolved lesions** — long-read assemblies of isolates whose
   specific frameshifts and nonsense mutations have been verified independently of the annotation.
2. **Functional assays at those loci**, not genome-wide claims: complementation, proteomics or
   transcriptomics tied to the individual gene.
3. **Congeneric partners in the same corpus**, so the orthologue-based evidence can be computed at all.
4. Realistically **tens, not thousands**, of such loci — and the audit is right that a modest number of
   high-quality examples beats thousands of circular ones.

The most tractable concrete target is a proteogenomics dataset for an organism with a congeneric
partner already annotated here, which would give Category A-adjacent evidence (expression, hence
*not* complete loss of function) at pseudogene-called loci. That is a separate study.

## 8 Honest summary for the editor

> We do not claim to have validated pseudogene predictions against experimental ground truth. We claim
> to have measured, with controls, how much bacterial pseudogene annotation depends on the pipeline
> that produced it — and we show that the loci where pipelines disagree carry sequence-level evidence
> of disruption far more often than loci where they agree. Whether those genes have lost function is a
> question our data cannot answer, and we say so.
