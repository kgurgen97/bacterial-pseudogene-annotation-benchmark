#!/usr/bin/env python3
"""POST HOC, 2026-10-30. Does the reduction gradient survive shared ancestry?

Fixed in PHYLO_PLAN.md before this file existed; the plan's SHA-256 is asserted
below. No other plan is touched.

C3 carries the paper and is tested across strata of genomes that are not
independent draws: the S1-plus-S4 contrast uses 95 genomes from 50 genera and
one genus supplies 22 of them.
"""
import hashlib
import json
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import stats

P = pathlib.Path(__file__).resolve().parent
R = P / "results"
PLAN = P / "PHYLO_PLAN.md"
FROZEN = "9b8e4b2ec030300ec844ad8dd15eeb65cccd13d32f5dc23af26436d1ac6a095b"

h = hashlib.sha256(PLAN.read_bytes()).hexdigest()
if h != FROZEN:
    sys.exit(f"PHYLO_PLAN.md changed since it was frozen\n  expected {FROZEN}\n"
             f"  found    {h}")
print(f"plan sha256 {h[:16]}… matches the frozen plan")

RNG = np.random.default_rng(20261030)
DRAWS = 2000           # plan section 2
PERMS = 5000
CALLERS = ["Pseudofinder", "PGAP", "DFAST", "Bakta"]
COUNT = {"PGAP": "pgap_pseudo", "Pseudofinder": "pf_pseudo",
         "DFAST": "dfast_disrupted", "Bakta": "bakta_pseudo"}
DENOM = {"PGAP": "pgap_cds", "Pseudofinder": "pf_intact",
         "DFAST": "dfast_cds", "Bakta": "bakta_cds"}
CONTRASTS = [("S3_recent_restriction", "S4_free_living_breadth", "S3_vs_S4"),
             ("S3_recent_restriction", "S1_extreme_obligate", "S3_vs_S1"),
             ("S1_extreme_obligate", "S4_free_living_breadth", "S1_vs_S4")]

# DEVIATION from PHYLO_PLAN.md section 2, reported in PHYLO_RESULTS.md.
# The plan names results/corpus_traits.tsv as the lineage source. That file
# covers the 363 original genomes and the four strata are the 160 expansion
# genomes, so the two sets do not intersect and the planned merge leaves every
# contrast empty. Genus is therefore parsed from the organism name by the same
# rule analysis_dnds.py uses, which resolves all 160. Family is lost with the
# named source, so the coarser check the plan promised cannot be run.
import re


def genus_of(o):
    if not isinstance(o, str):
        return None
    w = re.sub(r"[\[\]]", "", o).split()
    if not w:
        return None
    return w[1] if w[0] == "Candidatus" and len(w) > 1 else w[0]


d = pd.read_csv(R / "corpus523.tsv", sep="\t")
before = len(d)
d["genus"] = d.organism.map(genus_of)
d = d[d.genus.notna()]
print(f"{len(d)} of {before} genomes carry a genus; {d.genus.nunique()} genera")
print("  DEVIATION: genus parsed from the organism name, not corpus_traits.tsv,"
      " which does not cover the stratified genomes; family unavailable\n")

for c in CALLERS:
    d[f"f_{c}"] = d[COUNT[c]] / (d[COUNT[c]] + d[DENOM[c]])


def rb(x, y):
    u = stats.mannwhitneyu(x, y, alternative="two-sided")
    return 2 * u.statistic / (len(x) * len(y)) - 1, u.pvalue


out = {"plan_sha256": h, "n_genomes": len(d), "n_genera": int(d.genus.nunique()),
       "contrasts": {}, "verdicts": {}}

# ---------------------------------------------------------------- T4 first
print("=== T4  how clustered is the corpus?")
worst = 0.0
for g1, g2, label in CONTRASTS:
    sub = d[d.stratum.isin([g1, g2])]
    for s in (g1, g2):
        ss = sub[sub.stratum == s]
        if not len(ss):
            continue
        top = ss.genus.value_counts()
        frac = top.iloc[0] / len(ss)
        worst = max(worst, frac)
        print(f"  {s:26s} n {len(ss):3d}  largest genus {top.index[0]} "
              f"{top.iloc[0]} ({100 * frac:.0f}%)")
out["max_genus_share"] = float(worst)

for g1, g2, label in CONTRASTS:
    sub = d[d.stratum.isin([g1, g2])].copy()
    sub["grp"] = (sub.stratum == g1).astype(int)
    out["contrasts"][label] = {}
    print(f"\n=== {label}")
    for c in CALLERS:
        col = f"f_{c}"
        s = sub[[col, "grp", "genus"]].dropna()
        x = s[s.grp == 1][col].values
        y = s[s.grp == 0][col].values
        if len(x) < 5 or len(y) < 5:
            continue
        obs, p_raw = rb(x, y)

        # genus-clustered bootstrap
        genera = s.genus.unique()
        bykey = {k: v for k, v in s.groupby("genus")}
        vals = []
        for _ in range(DRAWS):
            draw = RNG.choice(genera, len(genera), replace=True)
            cc = pd.concat([bykey[k] for k in draw])
            a, b = cc[cc.grp == 1][col].values, cc[cc.grp == 0][col].values
            if len(a) >= 3 and len(b) >= 3:
                vals.append(rb(a, b)[0])
        lo, hi = (float(np.percentile(vals, 2.5)),
                  float(np.percentile(vals, 97.5))) if len(vals) > 50 else (np.nan, np.nan)

        # permutation restricted within genus
        informative = sum(1 for k, v in bykey.items() if v.grp.nunique() > 1)
        null = []
        for _ in range(PERMS):
            lab = s.groupby("genus").grp.transform(
                lambda v: RNG.permutation(v.values))
            a, b = s[lab == 1][col].values, s[lab == 0][col].values
            if len(a) >= 3 and len(b) >= 3:
                null.append(rb(a, b)[0])
        p_perm = (1 + sum(1 for v in null if abs(v) >= abs(obs))) / (1 + len(null)) \
            if null else np.nan

        # one genome per genus
        hits = 0
        for _ in range(DRAWS):
            pick = s.groupby("genus").sample(1, random_state=int(RNG.integers(1e9)))
            a, b = pick[pick.grp == 1][col].values, pick[pick.grp == 0][col].values
            if len(a) >= 3 and len(b) >= 3 and rb(a, b)[1] < 0.05:
                hits += 1
        frac_sig = hits / DRAWS

        out["contrasts"][label][c] = dict(
            rb=float(obs), p_pooled=float(p_raw), ci_lo=lo, ci_hi=hi,
            p_within_genus_perm=float(p_perm),
            frac_draws_significant=frac_sig,
            n_informative_genera=informative, n1=len(x), n2=len(y))
        print(f"  {c:13s} rb {obs:+.3f} (pooled p {p_raw:.3g})  "
              f"genus CI [{lo:+.3f},{hi:+.3f}]  perm p {p_perm:.4f}  "
              f"one-per-genus sig in {100 * frac_sig:.0f}% of draws  "
              f"({informative} informative genera)")


def verdict(k, pred, ok, bad, val):
    v = "CONFIRMED" if ok else ("FAILED" if bad else "INCONCLUSIVE")
    out["verdicts"][k] = dict(prediction=pred, value=val, verdict=v)
    print(f"  {k}  {v:13s} {val}")


print("\nverdicts, as fixed in PHYLO_PLAN.md")
c34 = out["contrasts"].get("S3_vs_S4", {})
meets = [c for c, j in c34.items()
         if not np.isnan(j["ci_lo"]) and j["ci_lo"] * j["ci_hi"] > 0
         and j["p_within_genus_perm"] < 0.05]
verdict("T1", "C3's S3 > S4 holds with ancestry accounted for",
        len(meets) >= 3, len(meets) < 2,
        f"{len(meets)} of {len(c34)} callers have a clustered interval "
        f"excluding zero and a within-genus permutation p < 0.05: "
        f"{', '.join(meets) if meets else 'none'}")

half = [c for c, j in c34.items() if j["frac_draws_significant"] >= 0.5]
quarter = [c for c, j in c34.items() if j["frac_draws_significant"] < 0.25]
verdict("T2", "the contrast is not carried by one over-represented lineage",
        len(half) >= 3, len(quarter) >= 3,
        "one-per-genus draws reaching p < 0.05: "
        + ", ".join(f"{c} {100 * j['frac_draws_significant']:.0f}%"
                    for c, j in c34.items()))

c14 = out["contrasts"].get("S1_vs_S4", {})
lost14 = sum(1 for j in c14.values() if j["p_within_genus_perm"] >= 0.05)
lost34 = sum(1 for j in c34.values() if j["p_within_genus_perm"] >= 0.05)
verdict("T3", "S1 vs S4 is the most lineage-dependent of the three",
        lost14 > lost34, lost14 <= lost34,
        f"callers losing significance under the within-genus permutation: "
        f"S1 vs S4 {lost14}, S3 vs S4 {lost34}")

verdict("T4", "the corpus is genuinely clustered",
        worst >= 0.10, worst < 0.10,
        f"the largest genus share of any stratum is {100 * worst:.0f}%")

(R / "phylo.json").write_text(json.dumps(out, indent=1))
nbad = sum(1 for v in out["verdicts"].values() if v["verdict"] != "CONFIRMED")
print(f"\nwrote results/phylo.json - {nbad} of {len(out['verdicts'])} not confirmed")
