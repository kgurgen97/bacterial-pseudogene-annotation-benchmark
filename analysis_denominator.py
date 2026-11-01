#!/usr/bin/env python3
"""POST HOC, 2026-10-29. Does the sign reversal depend on the denominator?

Not covered by any frozen plan, and labelled post hoc wherever it is reported.
No plan was frozen for it because it computes no new quantity: the counting rule
in LITERATURE_RETEST_PLAN.md already specified that raw counts and counts per Mb
be computed alongside the fraction, and this reports the contrasts on all three
instead of on the fraction alone. There is nothing here to pre-commit; a metric
cannot be chosen after the fact when all three are shown.

A second review of the manuscript observed that the fraction is compositional and
that the S1-versus-S4 contrast spans a tenfold difference in gene count, so the
denominator could be carrying the result. It is.

Also counts what DFAST's excluded `partial` class would add, which bears on the
section 3.6 claim that the spread is caller identity rather than configuration.
"""
import json
import pathlib

import numpy as np
import pandas as pd
from scipy import stats

P = pathlib.Path.home() / "research/projects/n2-pseudogene-caller-bias"
R = P / "results"
RNG = np.random.default_rng(20261029)

COUNT = {"PGAP": "pgap_pseudo", "Pseudofinder": "pf_pseudo",
         "DFAST": "dfast_disrupted", "Bakta": "bakta_pseudo"}
DENOM = {"PGAP": "pgap_cds", "Pseudofinder": "pf_intact",
         "DFAST": "dfast_cds", "Bakta": "bakta_cds"}
CALLERS = ["Pseudofinder", "PGAP", "DFAST", "Bakta"]

d = pd.read_csv(R / "corpus523.tsv", sep="\t")
d["mb"] = pd.to_numeric(d.size_bp, errors="coerce") / 1e6


def rb(x, y):
    """Rank-biserial and its p, the same effect size section 3.5 uses."""
    u = stats.mannwhitneyu(x, y, alternative="two-sided")
    return 2 * u.statistic / (len(x) * len(y)) - 1, u.pvalue


def boot(x, y, B=4000):
    out = np.empty(B)
    for b in range(B):
        out[b] = rb(RNG.choice(x, len(x), True), RNG.choice(y, len(y), True))[0]
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def series(sub, c, kind):
    if kind == "fraction":
        s = sub[COUNT[c]] / (sub[COUNT[c]] + sub[DENOM[c]])
    elif kind == "count":
        s = sub[COUNT[c]]
    else:
        s = sub[COUNT[c]] / sub["mb"]
    return s.replace([np.inf, -np.inf], np.nan).dropna().values


out = {"contrasts": {}, "note": "post hoc; no plan frozen, see the docstring"}
for g1, g2, label in (
        ("S1_extreme_obligate", "S4_free_living_breadth", "S1_vs_S4"),
        ("S3_recent_restriction", "S4_free_living_breadth", "S3_vs_S4"),
        ("S3_recent_restriction", "S1_extreme_obligate", "S3_vs_S1")):
    a, b = d[d.stratum == g1], d[d.stratum == g2]
    print(f"\n=== {label}   negative means the first group is lower")
    out["contrasts"][label] = {}
    for c in CALLERS:
        out["contrasts"][label][c] = {}
        line = f"  {c:13s}"
        for kind in ("fraction", "count", "per_Mb"):
            x, y = series(a, c, kind), series(b, c, kind)
            r, p = rb(x, y)
            lo, hi = boot(x, y)
            out["contrasts"][label][c][kind] = dict(
                rb=float(r), p=float(p), ci_lo=lo, ci_hi=hi,
                n1=len(x), n2=len(y),
                med1=float(np.median(x)), med2=float(np.median(y)))
            line += f"  {kind} {r:+.3f} [{lo:+.3f},{hi:+.3f}] p={p:.3g}"
        print(line)

# do the callers agree on sign, per metric?
print("\n=== sign agreement across the four callers")
for label in out["contrasts"]:
    for kind in ("fraction", "count", "per_Mb"):
        signs = {c: np.sign(out["contrasts"][label][c][kind]["rb"])
                 for c in CALLERS}
        sig = [c for c in CALLERS
               if out["contrasts"][label][c][kind]["p"] < 0.05]
        agree = len(set(signs.values())) == 1
        out["contrasts"][label].setdefault("_agreement", {})[kind] = dict(
            all_same_sign=bool(agree), n_significant=len(sig),
            significant=sig)
        print(f"  {label:10s} {kind:9s} same sign: {str(agree):5s}  "
              f"significant in {len(sig)} of 4")

print("\n=== DFAST's counting rule, against the configuration effect")
t = pd.read_csv(R / "dfast_counts.tsv", sep="\t")
counted = int(t.dfast_disrupted.sum())
partial = int(t.dfast_partial.sum())
cfg = json.loads((R / "pf_config.json").read_text())
out["dfast_rule"] = dict(
    counted=counted, excluded_partial=partial,
    ratio=partial / counted,
    shift_if_included=(counted + partial) / counted,
    pf_configuration_fold=cfg["within_caller_fold"],
    between_caller_fold=cfg["between_caller_fold"])
print(f"  counted as pseudogenes        {counted:>9,}")
print(f"  excluded as partial           {partial:>9,}")
print(f"  including them would give     {counted + partial:>9,}"
      f"  a {(counted + partial) / counted:.2f}x shift")
print(f"  the Pseudofinder parameter moves its caller {cfg['within_caller_fold']:.1f}x")
print(f"  the spread between callers is {cfg['between_caller_fold']:.1f}x")

(R / "denominator.json").write_text(json.dumps(out, indent=1))
print("\nwrote results/denominator.json")
