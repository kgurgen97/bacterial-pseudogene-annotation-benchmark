#!/usr/bin/env python3
"""POST HOC, added 2026-10-21 after the staged C1/C2 run. Not in LITERATURE_RETEST_PLAN.md.

Why: on the reliable 74-genome set, DFAST detects the canonical host-restriction
signature (p 0.030), PGAP is borderline (0.058) and Pseudofinder misses it entirely
(0.495) - despite Pseudofinder calling about six times more pseudogenes. The printed
IQRs suggest the reason: Pseudofinder's fraction sits in a narrow band around 6 % for
almost every genome, so it has little between-genome variance left to carry a signal.

This quantifies that. A caller that assigns nearly the same fraction to every genome
cannot discriminate between groups however many pseudogenes it finds, so dispersion -
not count - is what decides whether an association is detectable.
"""
import json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

P = Path.home() / "research/projects/n2-pseudogene-caller-bias"
d = pd.read_csv(P / "results/literature_retest.tsv", sep="\t")
cols = {c.replace("frac_", ""): c for c in d.columns if c.startswith("frac_")}
NAME = {"pgap": "PGAP", "pfinder": "Pseudofinder", "bakta": "Bakta", "dfast": "DFAST"}

print(f"{'caller':14s} {'n':>5s} {'median':>9s} {'IQR':>19s} {'CV':>7s} "
      f"{'IQR/median':>11s} {'Gini':>7s} {'range/med':>10s}")
rows = []
for key, col in cols.items():
    v = d[col].dropna()
    if len(v) < 50:
        continue
    q1, q3 = v.quantile(.25), v.quantile(.75)
    cv = v.std() / v.mean()
    s = np.sort(v.to_numpy())
    n = len(s)
    gini = (2 * np.arange(1, n + 1) - n - 1).dot(s) / (n * s.sum())
    rows.append(dict(caller=NAME.get(key, key), n=len(v), median=float(v.median()),
                     q1=float(q1), q3=float(q3), cv=float(cv),
                     iqr_over_median=float((q3 - q1) / v.median()),
                     gini=float(gini),
                     range_over_median=float((v.max() - v.min()) / v.median())))
    r = rows[-1]
    print(f"{r['caller']:14s} {r['n']:5d} {r['median']:9.4f} "
          f"{q1:8.4f}-{q3:<10.4f} {cv:7.3f} {r['iqr_over_median']:11.3f} "
          f"{gini:7.3f} {r['range_over_median']:10.1f}")

t = pd.DataFrame(rows)
print("\nRelative dispersion is what matters for detecting a group difference, not the")
print("absolute number of calls. Ranked by IQR/median (higher = more discriminating):")
for _, r in t.sort_values("iqr_over_median", ascending=False).iterrows():
    print(f"  {r.caller:14s} IQR/median {r.iqr_over_median:.3f}   "
          f"CV {r.cv:.3f}   Gini {r.gini:.3f}")

# does dispersion track whether the caller found C2?
c2 = json.load((P / "results/literature_retest.json").open()).get("C2_detail", {})
prim = {k.split("|")[1]: v for k, v in c2.items() if k.startswith("primary")}
if prim:
    print("\nAgainst the C2 result on the reliable 74-genome set:")
    print(f"  {'caller':14s} {'IQR/median':>11s} {'rank-biserial':>14s} {'MWU p':>10s}")
    for _, r in t.sort_values("iqr_over_median", ascending=False).iterrows():
        v = prim.get(r.caller)
        if v:
            print(f"  {r.caller:14s} {r.iqr_over_median:11.3f} "
                  f"{v['rbc']:+14.3f} {v['p']:10.3g}")
    xs = [r.iqr_over_median for _, r in t.iterrows() if r.caller in prim]
    ys = [abs(prim[r.caller]["rbc"]) for _, r in t.iterrows() if r.caller in prim]
    if len(xs) >= 3:
        rho, p_ = stats.spearmanr(xs, ys)
        print(f"  Spearman(dispersion, |effect found|) = {rho:+.3f}  p = {p_:.3g}  "
              f"(n={len(xs)} callers - far too few to test, shown as description only)")

print("\n=== how much do the callers agree per genome?")
keys = [k for k in cols if d[cols[k]].notna().sum() > 300]
for i, a in enumerate(keys):
    for b in keys[i + 1:]:
        ok = d[cols[a]].notna() & d[cols[b]].notna()
        r, p_ = stats.spearmanr(d.loc[ok, cols[a]], d.loc[ok, cols[b]])
        ratio = (d.loc[ok, cols[b]] / d.loc[ok, cols[a]]).median()
        print(f"  {NAME.get(a,a):14s} vs {NAME.get(b,b):14s} "
              f"Spearman {r:+.3f} (p {p_:.2g})   median ratio {ratio:6.2f}x")

t.to_csv(P / "results/caller_dispersion.tsv", sep="\t", index=False)
json.dump(rows, open(P / "results/caller_dispersion.json", "w"), indent=1)
print("\nwrote results/caller_dispersion.{tsv,json}")
