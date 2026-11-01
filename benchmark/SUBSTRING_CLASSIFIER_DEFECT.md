# Every DFAST "pseudogene" in this study is a product name containing the letters "pseudo"

2026-10-31. Found while writing the benchmark's output parser, which required reading what Bakta and
DFAST actually emit. **This invalidates DFAST's entire pseudogene call set and roughly a third of
Bakta's.** It is the most serious defect found in this project. No frozen plan is touched and
`MANUSCRIPT.md` is not edited, but the manuscript cannot be submitted as it stands.

## The defect

`analysis_dnds.py`, in `gff_intervals`, the function that turns a caller's GFF into pseudogene and
intact interval lists:

```python
is_pseudo = ("pseudo" in f[8].lower()) or f[2] == "pseudogene"
```

`f[8]` is the GFF attributes column. It holds `ID`, `Name`, `locus_tag`, `product`, `gene`, `Dbxref`,
`inference` and free-text `note`. A substring search for `pseudo` across all of that matches any gene
whose **product name or organism-derived annotation happens to contain those six letters**.

It should have tested for a pseudogene *assertion* — a `pseudogene=` or `pseudo=` attribute, or
feature type `pseudogene` — which is what the GenBank path in the same file does correctly:

```python
is_pseudo = ("pseudo" in ft.qualifiers) or ("pseudogene" in ft.qualifiers)   # dict KEY lookup
```

One is a key lookup on a qualifier dictionary. The other is a substring search on a concatenated
text blob. They sit forty lines apart.

## Which callers it touches

| caller | parsing path | sound? |
|---|---|---|
| PGAP | `gbff_intervals`, qualifier **key** lookup on GenBank features | **yes** |
| Pseudofinder | separate `pf_pseudos.gff` / `pf_intact.gff` files, no classification needed | **yes** (but see `PF_INPUT_DEPENDENCE.md`) |
| **Bakta** | `gff_intervals`, substring | **no** |
| **DFAST** | `gff_intervals`, substring | **no** |

## Magnitude, measured on the study's own output

`check_substring_classifier.py` replays the rule over the parent study's caller output and splits its
hits into real assertions and name matches. **All 523 genomes of the corpus:**

| caller | features called pseudogene | carrying a real `pseudo=` assertion | matched only by a name |
|---|---|---|---|
| **DFAST** | 125 307 | **0** | **125 307 — 100 %** |
| **Bakta** | 12 068 | 7 353 | 4 715 — 39.1 % |

**DFAST's GFF carries no pseudogene assertion** — no `pseudogene` feature type, no `pseudo=`
attribute — because DFAST asserts its pseudogenes in a **separate file**,
`pseudogene_summary.tsv`, with a `classification` column and explicit coordinates. Across the corpus
that file holds **45 976** disrupted calls. The interval path looked only at the GFF, so it missed
all 45 976 and substituted 125 307 name matches. The two errors do not partly cancel: they are
disjoint sets, in opposite directions, and the wrong one is 2.7 times larger.

### Correction, same day

An earlier revision of this file asserted that DFAST "emits no pseudogene assertion at all in this
configuration" and that its `PseudoGeneDetection` output was absent from the runs. **Both were
wrong.** `pseudogene_summary.tsv` is present for all 523 genomes, and the transient
`PseudoGeneDetection` working directory is simply cleaned up at the end of a run. The defect is real
and confined to the interval path; the claim that DFAST has no pseudogene output was not.

What triggered the spurious hits, pooled over all 523 genomes. These are **word occurrences, not
features** — one product name can carry two, as `pseudouridine(32) … pseudouridine(746)` does — so
they sum to more than the feature counts above:

| word | occurrences, DFAST | occurrences, Bakta |
|---|---|---|
| `pseudomonas` | 91 676 | — |
| `pseudotuberculosis` | 15 843 | — |
| `pseudouridine` | 11 342 | 8 010 |
| `pseudomallei` | 10 463 | — |
| `pseudopilin` | 2 024 | 1 198 |
| `pseudomonadota` | 1 465 | — |
| `pseudouridylate` | 943 | 104 |
| `pseudoazurin` | — | 100 |

DFAST's output names homologues by source organism, which is why *Pseudomonas*, *Yersinia
pseudotuberculosis* and *Burkholderia pseudomallei* dominate it and Bakta's contamination is almost
entirely the tRNA pseudouridine synthases.

A representative "DFAST pseudogene":

```
sequence1  CDS  57645-58274  product=bifunctional tRNA pseudouridine(32) synthase/
                             23S rRNA pseudouridine(746) synthase RluA
```

That is an intact, conserved, essential housekeeping enzyme.

## It explains every DFAST result in the project, and that is the point

`DFAST_RESULTS.md` reported four things about DFAST's 2 639 unique calls and treated them as a
puzzle. Each is exactly what this defect predicts:

| observed | explanation |
|---|---|
| DFAST's unique set is the **longest** (median 292 codons vs 281 baseline) | synthases and large named enzymes are long |
| it is **depleted of ω > 1 six-fold** (1 locus of 2 639, ~7 expected), and length-matching did not explain it | they are ordinary genes under normal constraint, not decaying ones |
| selection **does not corroborate** them (+0.006, p = 0.06; 15/26 genomes in the reduced strata) | there is nothing to corroborate; they were never broken |
| DFAST is **sole** pseudogene caller at 2 639 loci and opposed to all three others | no other caller agrees because there is nothing to agree with |

And the stratum concentration — DFAST-only loci at 91.5 per 1 000 in S1 and 44.3 in S3 against 12.4
in S4 — is an artefact of **organism names**. `pseudotuberculosis` alone contributes 3 204 hits, and
*Yersinia pseudotuberculosis* sits in S3_recent_restriction. The "concentration in reduced genomes"
is a concentration of genus names containing the string.

The honest reading of `DFAST_RESULTS.md` is therefore that it successfully characterised a parsing
bug and mistook it for a property of an annotation tool. Its own § 5 said the ω > 1 depletion was
"the one result that is harder to explain" by a caller tuned to recent disruption. It was harder to
explain because the explanation was not biological.

## The project already knew this trap, in the other half of the code

This is the part that locates the error precisely. `run_dfast_all.sh` carries this in its header:

> do not grep the GFF for "pseudo". It matches Pseudomonas, pseudomallei, pseudouridine and
> pseudopilin. The summary TSV is the only authoritative source.

and `run_bakta_all.sh`:

> Do NOT grep the TSV for "pseudo": that matches products such as "23S rRNA
> (pseudouridine(1915)-N(3))-methyltransferase", which gave 21 instead of 7

Both runners take their **per-genome counts** from each caller's own summary file, and both are
right. The warning never reached the **interval** path in `analysis_dnds.py`. So the project contains
a correct implementation and an incorrect one, and which results are affected follows exactly from
which of the two they were built on.

## What is invalidated, and what is not

**Sound — count-based, from the callers' own summary files:** `results/dfast_counts.tsv`,
`results/bakta_counts.tsv` and everything downstream of them. That means
`LITERATURE_RETEST_RESULTS.md` including the 19-fold spread in raw counts,
`analysis_caller_dispersion.py`, `analysis_c3_strata.py`'s C3 stratum test, `analysis_denominator.py`,
`PF_CONFIG_RESULTS.md`, and `PHYLO_RESULTS.md`, which reads `corpus523.tsv`. The § 3.5
PGAP-versus-DFAST sign reversal is a count result and is **not** affected.

**Invalidated — interval-based, through `gff_intervals`:**

- **`DFAST_RESULTS.md` in full.** The pair matrix, the sole-caller comparison, the shape table and
  the stratum table all come from `loci.tsv`.
- **`LOCUS_RESULTS.md`'s discordant set.** DFAST appears on the pseudogene side of 3 169 of 7 565
  discordant loci and Bakta of 686, so **G3** (+0.022 per genome, p = 4 × 10⁻¹²), **L1** and **L3**
  must be recomputed before they can be quoted.
- **`DNDS_RESULTS.md`'s** G1, G2 and G3, and **`FRAME_RESULTS.md`'s** F-tests, both of which call
  `calls_for`.
- **`CLUSTERING_RESULTS.md`**, which reads `dnds_pairs.tsv`.

## What survives

- **PGAP's calls**, parsed by the sound GenBank path.
- **Pseudofinder's calls**, read from its own two files.
- Therefore the strongest result in the project survives: the loci where Pseudofinder alone calls a
  pseudogene and PGAP calls the gene an intact CDS sit well above the concordant-intact baseline.
  Both sides of that comparison come from sound code paths.

  **Updated after the recomputation, 2026-10-31.** This file originally quoted the pre-fix figures —
  1 955 loci, +0.035, 60 of 67 genomes, p = 1 × 10⁻¹¹ — and called them "unchanged". An adversarial
  review caught that: the set itself moved by 11 %. Corrected: **1 746 loci, +0.0379, 62 of 68
  genomes, p = 8.2 × 10⁻¹³** (`results/caller_comparison.json`); the L2 count against PGAP-intact is
  **1 742** (`results/locus.json`). The direction and strength survive — they strengthen slightly —
  but "unchanged" was wrong.
- `PF_CONFIG_RESULTS.md`, which varied Pseudofinder's own parameters on a fixed input.
- The **simulation benchmark** is unaffected and becomes more useful, not less: it has its own
  parser, built from the vocabulary actually present in the output, and DFAST's real disruption
  signals (`note=...partial/frameshifted/internal stop`) are in its plan as endpoint B rather than
  being read as pseudogene assertions.

## The fix, and how it was verified

`gff_asserts_pseudo` now tests for an attribute whose **key** is `pseudo` or `pseudogene`, and
`dfast_intervals` reads DFAST's pseudogenes from `pseudogene_summary.tsv` — taking the disrupted
classes `frameshift`, `internal_stop_codon` and `frameshift,internal_stop_codon`, and **excluding
`partial`**, which marks low coverage against a reference and is a length judgement rather than a
disruption. That is the same decision `run_dfast_all.sh` made for the counts.

The fix has an independent oracle: the runners' awk counts, written months earlier from different
files. `check_classifier_fix.py` compares them across all 523 genomes:

| | corrected | defective |
|---|---|---|
| DFAST pseudogene intervals | **45 976** | 125 307 |
| Bakta pseudogene intervals | **7 353** | 12 068 |

**Zero mismatches in 523 genomes.** DFAST's corrected count equals `dfast_counts.tsv`'s
`dfast_disrupted` exactly, genome by genome. Bakta decomposes into 6 844 rows carrying
`pseudogene=`, which matches its own `pseudogenes:` summary line exactly, plus **509** rows carrying
`pseudo=True`, a second key Bakta uses and its summary line omits. Both are genuine assertions so
both are kept as intervals, and the decomposition is reported rather than hidden inside a total.

Bakta also marks a handful of **tRNA** features `pseudogene=unknown`; those are excluded by the
feature-type filter, which is what makes the strict subset match the summary.

## What must still happen

1. Recompute the interval analyses above and compare against what the write-ups published.
2. Only then revisit the manuscript.

This is not a defect the benchmark was built to find. It was found because building the benchmark
required reading the callers' real output rather than trusting a rule written months earlier, which
is an argument for the benchmark that the novelty audit did not anticipate.

## Files

`check_substring_classifier.py` → `results/substring_classifier.tsv`,
`results/substring_classifier.json` · the defect is at `analysis_dnds.py` in `gff_intervals`, and the
correct pattern is in `gbff_intervals` in the same file.
