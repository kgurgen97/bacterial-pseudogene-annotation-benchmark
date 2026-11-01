#!/usr/bin/env python3
"""The 19-fold and the 77-fold are different metrics. Both can be right.

MANUSCRIPT.md says "The median pseudogene fraction spans 77.2-fold across the
four, and the corpus totals span nineteenfold." Those are not in tension and
they are not alternative estimates of one quantity:

  corpus totals     sum over all 523 genomes of raw pseudogene calls, then
                    max/min across callers. Dominated by the pseudogene-rich
                    genomes, since a genome with 2,000 calls contributes 2,000.

  median fraction   per genome, pseudogene calls divided by that genome's gene
                    count; then the MEDIAN of that across genomes, then max/min
                    across callers. Dominated by the typical genome, and
                    insensitive to how extreme the rich ones are.

They diverge here for a specific, statable reason: Bakta reports zero
pseudogenes in many genomes, so its MEDIAN fraction sits on a floor while its
TOTAL is still non-trivial. A ratio of medians against a floor is larger than a
ratio of sums.

This recomputes both from the count files, which come from each caller's own
summary output and were never touched by the classifier defect, and prints the
arithmetic so the two numbers can be cited without being conflated.

    python3 benchmark/check_two_spreads.py
"""
import json
import pathlib

import numpy as np
import pandas as pd

P = pathlib.Path(__file__).resolve().parent.parent
R = P / "results"


def load():
    """Per-genome pseudogene counts and gene denominators for all four callers."""
    frames = {}
    bk = pd.read_csv(R / "bakta_counts.tsv", sep="\t")
    frames["Bakta"] = pd.DataFrame({
        "accession": bk.accession, "n": bk.bakta_pseudo,
        "genes": bk.bakta_cds + bk.bakta_pseudo})
    df = pd.read_csv(R / "dfast_counts.tsv", sep="\t")
    frames["DFAST"] = pd.DataFrame({
        "accession": df.accession, "n": df.dfast_disrupted,
        "genes": df.dfast_cds})
    # PGAP and Pseudofinder come from the corpus table the other analyses use
    c = pd.read_csv(R / "corpus523.tsv", sep="\t")
    cols = {x.lower(): x for x in c.columns}

    def pick(*cands):
        for x in cands:
            if x in cols:
                return cols[x]
        return None

    for caller, ncands, gcands in (
            ("PGAP", ("pgap_pseudo", "pseudo", "n_pseudo"),
             ("pgap_cds", "cds", "n_cds")),
            ("Pseudofinder", ("pf_pseudo", "pseudofinder_pseudo"),
             ("pf_intact", "pseudofinder_intact", "pf_cds"))):
        nc, gc = pick(*ncands), pick(*gcands)
        if nc and gc:
            frames[caller] = pd.DataFrame({
                "accession": c.accession, "n": c[nc],
                "genes": c[nc] + c[gc] if "intact" in gc else c[gc]})
    return frames


def main() -> None:
    frames = load()
    print(f"callers with per-genome counts: {sorted(frames)}\n")
    rows = []
    for caller, t in sorted(frames.items()):
        t = t.dropna()
        t = t[t.genes > 0]
        frac = t.n / t.genes
        rows.append({"caller": caller, "genomes": len(t),
                     "corpus_total": int(t.n.sum()),
                     "median_fraction": float(np.median(frac)),
                     "mean_fraction": float(np.mean(frac)),
                     "genomes_with_zero": int((t.n == 0).sum())})
    d = pd.DataFrame(rows).sort_values("median_fraction")
    print(d.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    tot = d.corpus_total
    med = d.median_fraction
    r_tot = tot.max() / max(tot.min(), 1)
    nonzero_med = med[med > 0]
    r_med = med.max() / nonzero_med.min() if len(nonzero_med) else float("nan")

    print(f"\ncorpus totals      : {int(tot.min())} .. {int(tot.max())}"
          f"   ratio {r_tot:.1f}-fold")
    print(f"median fractions   : {med.min():.4f} .. {med.max():.4f}"
          f"   ratio {r_med:.1f}-fold")
    print(f"\n  the ratio of medians is {r_med / r_tot:.1f}x the ratio of totals")
    floor = d[d.median_fraction == med.min()]
    if len(floor):
        f = floor.iloc[0]
        print(f"  the smaller median belongs to {f.caller}, which reports zero "
              f"pseudogenes in\n  {f.genomes_with_zero} of {f.genomes} genomes "
              f"-- a median against a floor, which is why a ratio of\n  medians "
              f"exceeds a ratio of sums. Both numbers are correct; they answer "
              f"different\n  questions and must be labelled in the text.")

    out = P / "benchmark/results/two_spreads.json"
    out.write_text(json.dumps(
        {"per_caller": rows,
         "corpus_total_ratio": float(r_tot),
         "median_fraction_ratio": float(r_med)}, indent=2, sort_keys=True) + "\n")
    print(f"\nwrote {out.relative_to(P)}")


if __name__ == "__main__":
    main()
