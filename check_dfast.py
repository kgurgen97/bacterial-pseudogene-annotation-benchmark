#!/usr/bin/env python3
"""Two ways the DFAST result could be an artefact, and one arithmetic check.

C1  DFAST's unique loci are depleted of omega > 1 (0.04% against a baseline
    0.25%). They are also the longest set (median 292 codons). Longer genes
    give omega less room to take an extreme value, so the depletion may be
    nothing but length. Compare against a length-matched slice of the
    baseline before claiming anything about it.

C2  DFAST's unique calls concentrate in S1 and S3, the reduced and recently
    restricted genomes that section 3.5 argues about. Is the weak selection
    signal weak everywhere, or only there? This is an exploratory split on a
    grouping the project already uses (S1+S3 against S2+S4), reported with no
    verdict attached.

C3  The pair matrix says 1 990 loci have Pseudofinder on the pseudogene side
    and PGAP on the intact side; LOCUS_RESULTS.md says 1 955 for the
    sole-Pseudofinder version, and the sole-Pseudofinder set is 1 960. Those
    three numbers must reconcile exactly or one of them is wrong.
"""
import pathlib

import numpy as np
import pandas as pd
from scipy import stats

P = pathlib.Path(__file__).resolve().parent
t = pd.read_csv(P / "results/loci.tsv", sep="\t")
dis = t[t.kind == "discordant"].copy()
con = t[t.kind == "concordant_intact"].copy()
dfo = dis[dis.pseudo_callers == "DFAST"]

print("=== C1  is the omega>1 depletion in DFAST's unique set just length?")
lo, hi = dfo.codons.quantile(0.25), dfo.codons.quantile(0.75)
m = con[(con.codons >= lo) & (con.codons <= hi)]
print(f"  DFAST-only          n {len(dfo):6d}  median codons "
      f"{dfo.codons.median():5.0f}  frac omega>1 {(dfo.omega > 1).mean():.4f}")
print(f"  baseline, all       n {len(con):6d}  median codons "
      f"{con.codons.median():5.0f}  frac omega>1 {(con.omega > 1).mean():.4f}")
print(f"  baseline, {lo:.0f}-{hi:.0f} codons  n {len(m):6d}  median codons "
      f"{m.codons.median():5.0f}  frac omega>1 {(m.omega > 1).mean():.4f}")
print("  -> the length-matched baseline is the only fair comparison; if it is")
print("     close to DFAST's figure, the depletion is length and not biology.")

print("\n=== C2  exploratory split, no verdict")


def sign_test(sub, label):
    g = []
    for acc in set(sub.accession) & set(con.accession):
        x = sub.loc[sub.accession == acc, "omega"]
        y = con.loc[con.accession == acc, "omega"]
        if len(x) >= 3 and len(y) >= 3:
            g.append(float(np.median(x) - np.median(y)))
    nz = [v for v in g if v != 0.0]
    if len(nz) < 5:
        print(f"  {label:34s} {len(sub):5d} loci, only {len(nz)} genomes, not tested")
        return
    pos = sum(1 for v in nz if v > 0)
    p = stats.binomtest(pos, len(nz), 0.5).pvalue
    print(f"  {label:34s} {len(sub):5d} loci  median {np.median(g):+.4f}  "
          f"{pos}/{len(nz)} positive  p {p:.3g}")


RED = ["S1_extreme_obligate", "S3_recent_restriction"]
sign_test(dfo, "DFAST-only, all strata")
sign_test(dfo[dfo.stratum.isin(RED)], "DFAST-only, S1+S3 reduced")
sign_test(dfo[~dfo.stratum.isin(RED)], "DFAST-only, S2+S4")
pfo = dis[dis.pseudo_callers == "Pseudofinder"]
sign_test(pfo[pfo.stratum.isin(RED)], "Pseudofinder-only, S1+S3 (contrast)")
sign_test(pfo[~pfo.stratum.isin(RED)], "Pseudofinder-only, S2+S4 (contrast)")

print("\n=== C3  do 1990, 1960 and 1955 reconcile?")
sets = dis.intact_callers.fillna("").apply(lambda s: set(s.split("+")) - {""})
any_pf = dis.pseudo_callers.str.contains("Pseudofinder") & sets.apply(
    lambda s: "PGAP" in s)
sole_pf = dis.pseudo_callers == "Pseudofinder"
sole_pf_vs_pgap = sole_pf & sets.apply(lambda s: "PGAP" in s)
print(f"  Pseudofinder anywhere on pseudo side, PGAP intact : {int(any_pf.sum())}")
print(f"  Pseudofinder SOLE pseudogene caller               : {int(sole_pf.sum())}")
print(f"  Pseudofinder sole, PGAP on the intact side        : {int(sole_pf_vs_pgap.sum())}")
rest = dis[sole_pf & ~sole_pf_vs_pgap]
print(f"  the remainder, sole PF but PGAP NOT intact        : {len(rest)}"
      f"  intact sides {rest.intact_callers.value_counts().to_dict()}")
ok = (int(sole_pf.sum()) == int(sole_pf_vs_pgap.sum()) + len(rest)
      and int(sole_pf_vs_pgap.sum()) == 1955)
print(f"  reconciles and matches the write-up: {ok}")
if not ok:
    raise SystemExit("the three numbers do not reconcile")
