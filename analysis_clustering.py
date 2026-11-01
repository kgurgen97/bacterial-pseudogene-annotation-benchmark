#!/usr/bin/env python3
"""POST HOC, 2026-10-30. How much of the dN/dS significance is pseudoreplication?

Not covered by any frozen plan, and labelled post hoc wherever reported. No plan
was frozen because this computes no new estimand: it re-expresses the section 3.7
contrast with the clustering respected, which an audit pointed out the original
test ignores.

analysis_dnds.py pools 33,850 sequences from 99 genomes and 29 genera and runs
one Mann-Whitney, which treats sequences from the same genome as independent.
The median genome contributes 166 of them and one contributes 3,589, so the
p-values are inflated by an unknown but certainly large factor. The direction of
the effect is a separate question from the size of its p-value, and this asks
whether the direction survives when the genome, not the sequence, is the unit.

Three ways of respecting the cluster structure:

  1. per-genome gaps, then a sign test and a Wilcoxon over genomes
  2. a cluster bootstrap resampling whole genomes
  3. the same resampling whole genera, which is the more conservative unit
"""
import json
import pathlib

import numpy as np
import pandas as pd
from scipy import stats

P = pathlib.Path.home() / "research/projects/n2-pseudogene-caller-bias"
R = P / "results"
RNG = np.random.default_rng(20261030)
CALLERS = ["Pseudofinder", "PGAP", "Bakta", "DFAST"]
MIN_PER_SIDE = 5          # a genome needs this many calls and controls to count

t = pd.read_csv(R / "dnds_pairs.tsv", sep="\t")
print(f"{len(t)} sequences, {t.accession.nunique()} genomes, "
      f"{t.genus.nunique()} genera")
print(f"sequences per genome: median {int(t.groupby('accession').size().median())}, "
      f"max {int(t.groupby('accession').size().max())}\n")

out = {"note": "post hoc; re-expresses section 3.7 with clustering respected",
       "min_per_side": MIN_PER_SIDE, "per_caller": {}}


def genome_gaps(sub):
    """One gap per genome, from genomes carrying enough of both classes."""
    g = []
    for acc, s in sub.groupby("accession"):
        a = s[s.call == "pseudo"].omega
        b = s[s.call == "intact"].omega
        if len(a) >= MIN_PER_SIDE and len(b) >= MIN_PER_SIDE:
            g.append((acc, float(np.median(a) - np.median(b))))
    return pd.DataFrame(g, columns=["accession", "gap"])


def cluster_boot(sub, unit, B=4000):
    """Resample whole clusters, recomputing the pooled gap each time."""
    keys = sub[unit].unique()
    byk = {k: s for k, s in sub.groupby(unit)}
    vals = np.empty(B)
    for i in range(B):
        draw = RNG.choice(keys, len(keys), replace=True)
        c = pd.concat([byk[k] for k in draw])
        a = c[c.call == "pseudo"].omega
        b = c[c.call == "intact"].omega
        vals[i] = np.median(a) - np.median(b) if len(a) and len(b) else np.nan
    v = vals[~np.isnan(vals)]
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


print(f"{'caller':13s} {'n gen':>6s} {'median':>8s} {'pos':>7s} "
      f"{'sign p':>9s} {'wilcox p':>9s} {'genome CI':>20s} {'genus CI':>20s}")
for c in CALLERS:
    sub = t[t.caller == c]
    gg = genome_gaps(sub)
    if len(gg) < 5:
        print(f"{c:13s} too few genomes with both classes ({len(gg)})")
        continue
    pos = int((gg.gap > 0).sum())
    sign_p = stats.binomtest(pos, len(gg), 0.5, alternative="greater").pvalue
    wil_p = stats.wilcoxon(gg.gap, alternative="greater").pvalue
    glo, ghi = cluster_boot(sub, "accession")
    nlo, nhi = cluster_boot(sub, "genus")
    print(f"{c:13s} {len(gg):6d} {np.median(gg.gap):+8.3f} "
          f"{pos:3d}/{len(gg):<3d} {sign_p:9.3g} {wil_p:9.3g} "
          f"[{glo:+.3f},{ghi:+.3f}]".ljust(0) + f"  [{nlo:+.3f},{nhi:+.3f}]")
    out["per_caller"][c] = dict(
        n_genomes=len(gg), median_genome_gap=float(np.median(gg.gap)),
        n_positive=pos, sign_test_p=float(sign_p), wilcoxon_p=float(wil_p),
        genome_cluster_ci=[glo, ghi], genus_cluster_ci=[nlo, nhi],
        pooled_p_for_comparison=None)

full = json.loads((R / "dnds.json").read_text())["per_caller"]
for c in out["per_caller"]:
    out["per_caller"][c]["pooled_p_for_comparison"] = full[c]["p"]
    out["per_caller"][c]["pooled_gap_for_comparison"] = full[c]["gap"]

print("\nwhat changes when the genome rather than the sequence is the unit")
for c in out["per_caller"]:
    o = out["per_caller"][c]
    lo, hi = o["genus_cluster_ci"]
    excl = "excludes 0" if lo > 0 or hi < 0 else "INCLUDES 0"
    print(f"  {c:13s} pooled p {o['pooled_p_for_comparison']:.3g} -> "
          f"sign test p {o['sign_test_p']:.3g}; genus-cluster CI {excl}")

(R / "clustering.json").write_text(json.dumps(out, indent=1))
print("\nwrote results/clustering.json")
