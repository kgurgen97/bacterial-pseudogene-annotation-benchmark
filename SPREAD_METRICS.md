# The 19-fold and the 77-fold are different measurements, and both are right

2026-10-31. Written because the two numbers appear in the same sentence of § 3.1 and then separately
elsewhere, where "77-fold" alone can read as *the* spread. Recomputed from the count files by
`benchmark/check_two_spreads.py`; **`benchmark/results/two_spreads.json`** holds the output — this line
previously gave the path as `results/two_spreads.json`, where no such file exists, while § Files below
always had it right. Nothing here was
touched by the classifier defect — these come from each caller's own summary output, which was always
parsed correctly.

## The two numbers

| caller | corpus total | median per-genome fraction | genomes with **zero** calls |
|---|---|---|---|
| Bakta | 6 844 | **0.0008** | **128** of 523 |
| PGAP | 37 493 | 0.0111 | 5 |
| DFAST | 45 976 | 0.0123 | 3 |
| Pseudofinder | 130 748 | **0.0635** | 0 |

- **Corpus totals span 19.1-fold** — 130 748 against 6 844. Every pseudogene call in the corpus
  summed, then the largest caller over the smallest.
- **Median per-genome fractions span 77.2-fold** — 0.0635 against 0.0008. Per genome, calls divided
  by that genome's gene count; the median of that across 523 genomes; then largest over smallest.

## Why they differ by 4-fold, stated rather than hand-waved

They are not two estimates of one quantity. A corpus total is dominated by the pseudogene-rich
genomes, because a genome contributing 2 000 calls contributes all 2 000. A median fraction is
dominated by the *typical* genome and is indifferent to how extreme the rich ones are.

The gap is almost entirely one fact: **Bakta reports no pseudogene at all in 128 of 523 genomes.**
Its median fraction therefore sits essentially on a floor (0.0008, against a mean of 0.0053 — a
seven-fold gap between its own mean and median), while its corpus total of 6 844 is not small. A
ratio taken against a floor is larger than a ratio of sums. Pseudofinder, at the other end, reports
zero in **no** genome.

So the honest reading is that the 77.2-fold figure measures disagreement *about the typical genome*
and is inflated by one caller's silence, while the 19.1-fold figure measures disagreement *about the
corpus* and is not. Neither is wrong; quoting one without naming which it is, is.

## How each should be written

- Always name the metric: "77.2-fold in **median per-genome pseudogene fraction**" and "19.1-fold in
  **corpus totals**", never a bare "77-fold spread between callers".
- Where only one number is given — the abstract, the cover letter, § 4 — it must carry its metric
  name, and the floor should be mentioned at least once, because 128 genomes of silence is the reason
  the larger number is larger.
- The 1.5-fold Pseudofinder configuration shift (`PF_CONFIG_RESULTS.md`) is a shift in **median
  fraction**, so it is comparable to 77.2-fold and not to 19.1-fold. The existing sentence that sets
  1.5-fold against 77-fold is therefore comparing like with like, and should say so.

## One definitional choice that moves the 19.1-fold figure

`bakta_counts.tsv` takes Bakta's total from Bakta's own `pseudogenes:` summary line, which counts the
**6 844** CDS rows carrying `pseudogene=`. Bakta's GFF also marks **509** CDS rows `pseudo=True`, a
second key its summary line omits (`benchmark/results/parser_audit.json`). Both are genuine
assertions, so the corrected *interval* analyses use all **7 353**.

That leaves one caller with two defensible totals, and the spread depends on which is used:

| Bakta total used | corpus-total spread |
|---|---|
| 6 844, Bakta's own summary definition | **19.1-fold** |
| 7 353, every assertion in its GFF | **17.8-fold** |

The difference is 6.9 % of Bakta's calls and it does not change any sign, ordering or verdict. But the
count-based tables and the interval-based analyses currently use different definitions for the same
caller, and a reviewer will find that. It needs one decision, applied consistently, and recorded —
which is a decision about the manuscript and is deliberately not taken here.

## Files

`benchmark/check_two_spreads.py` → `benchmark/results/two_spreads.json` · inputs
`results/bakta_counts.tsv`, `results/dfast_counts.tsv`, `results/corpus523.tsv`, all from the
callers' own summary output.
