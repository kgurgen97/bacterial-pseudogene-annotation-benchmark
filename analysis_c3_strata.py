#!/usr/bin/env python3
"""C3, the non-monotonicity test, exactly as CORPUS_EXPANSION_PLAN.md specifies.

Plan sha256: 1b1ee35c19bd66388ba7b952de41abab840e08ebbd4656bf08e47aef83e8c5fa

C3: pseudogene fraction is non-monotonic in genome reduction - low in free-living,
HIGHEST in recently host-restricted lineages, low again in extremely reduced obligate
endosymbionts. Predicted ordering S3 > S2 > S4 ~ S1.

Declared decision rule, from the plan: "Non-monotonicity is declared confirmed only if
S3 > S1 significantly AND S3 > S4 significantly for at least three of four callers."

C3 is tested on the 160 additions, because the strata are defined for them. The original
363 are shown alongside as context only, never pooled into a stratum.
"""
import json, re
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

P = Path.home() / "research/projects/n2-pseudogene-caller-bias"
R = P / "results"
PF = Path.home() / "research/datasets/n2-pseudofinder"
ALPHA = 0.05

# ---------------------------------------------------------------- Pseudofinder counts
# Rule verified against the recorded values for GCF_000008325.1: non-comment lines in
# pf_pseudos.gff = 194 = pf_pseudo, and FASTA headers in pf_intact.faa = 2849 = pf_intact.
def pf_counts(acc):
    g, f = PF / acc / "pf_pseudos.gff", PF / acc / "pf_intact.faa"
    if not (g.exists() and f.exists()):
        return None, None
    ps = sum(1 for l in g.open(encoding="utf-8", errors="replace")
             if l.strip() and not l.startswith("#"))
    it = sum(1 for l in f.open(encoding="utf-8", errors="replace") if l.startswith(">"))
    return ps, it


old = pd.read_csv(R / "pseudogene_counts.tsv", sep="\t")
exp = pd.read_csv(P / "manifests/expansion_pgap.tsv", sep="\t")

rows = []
for acc in sorted(set(old.accession) | set(exp.accession)):
    ps, it = pf_counts(acc)
    rows.append(dict(accession=acc, pf_pseudo=ps, pf_intact=it))
pf = pd.DataFrame(rows)

# integrity check: recomputed values must match what was recorded for the original 363
chk = old[["accession", "pf_pseudo", "pf_intact"]].merge(pf, on="accession",
                                                         suffixes=("_rec", "_new"))
bad = chk[(chk.pf_pseudo_rec != chk.pf_pseudo_new) | (chk.pf_intact_rec != chk.pf_intact_new)]
print(f"Pseudofinder recount vs recorded, original corpus: {len(chk)} compared, "
      f"{len(bad)} mismatched")
if len(bad):
    print(bad.head(8).to_string(index=False))

# ---------------------------------------------------------------- assemble 523
d = pd.DataFrame({"accession": sorted(set(old.accession) | set(exp.accession))})
d = d.merge(pf, on="accession", how="left")
d = d.merge(old[["accession", "pgap_pseudo", "pgap_cds"]], on="accession", how="left")
e = exp[["accession", "pgap_pseudo", "pgap_cds", "stratum", "species", "size_bp",
         "gc_percent", "organism"]].copy()
d = d.merge(e, on="accession", how="left", suffixes=("", "_exp"))
d["pgap_pseudo"] = d.pgap_pseudo.fillna(d.pgap_pseudo_exp)
d["pgap_cds"] = d.pgap_cds.fillna(d.pgap_cds_exp)
d["stratum"] = d.stratum.fillna("original_corpus")

df_ = pd.read_csv(R / "dfast_counts.tsv", sep="\t")
d = d.merge(df_[["accession", "dfast_disrupted", "dfast_cds"]], on="accession", how="left")
bk = pd.read_csv(R / "bakta_counts.tsv", sep="\t")
d = d.merge(bk[["accession", "bakta_pseudo", "bakta_cds"]], on="accession", how="left")

d["frac_PGAP"] = d.pgap_pseudo / (d.pgap_pseudo + d.pgap_cds)
d["frac_Pseudofinder"] = d.pf_pseudo / (d.pf_pseudo + d.pf_intact)
d["frac_DFAST"] = d.dfast_disrupted / (d.dfast_disrupted + d.dfast_cds)
d["frac_Bakta"] = d.bakta_pseudo / (d.bakta_pseudo + d.bakta_cds)

S = ["S4_free_living_breadth", "S1_extreme_obligate",
     "S2_obligate_intracellular", "S3_recent_restriction"]
exp_only = d[d.stratum.isin(S)]

# C3 is tested on the additions alone, so what matters is a caller's coverage of
# THOSE genomes, not of the 523. The earlier gate was ">= 400 of 523", which Bakta
# would have crossed while holding only part of the additions - and unevenly, since
# it works through the corpus in its own order. A caller is admitted here only if it
# covers every addition, so all callers are compared on one genome set as the plan
# requires. Coverage is printed either way, so an exclusion is visible rather than
# silent.
FRACS = [c for c in d.columns if c.startswith("frac_")]
print(f"\ngenomes: {len(d)}   additions C3 is tested on: {len(exp_only)}")
print(f"  {'caller':22s} {'of 523':>7s} {'of additions':>13s}  admitted")
CALLERS = []
for c in FRACS:
    nm = c.replace("frac_", "")
    n_all = int(d[c].notna().sum())
    n_exp = int(exp_only[c].notna().sum())
    ok = n_exp == len(exp_only)
    if ok:
        CALLERS.append(nm)
    print(f"  {nm:22s} {n_all:7d} {n_exp:13d}  {'yes' if ok else 'NO - incomplete'}")
assert CALLERS, "no caller covers all the additions; nothing can be tested yet"
print(f"\nC3 is tested on the {len(exp_only)} additions; "
      f"the {int((d.stratum=='original_corpus').sum())} original genomes are context only")

# ---------------------------------------------------------------- C3
print("\n" + "=" * 92)
print("C3  pseudogene fraction by ecological stratum.  PREDICTED S3 > S2 > S4 ~ S1")
print("=" * 92)
out = {"n_additions": int(len(exp_only)), "callers": CALLERS, "per_caller": {}}
for nm in CALLERS:
    col = f"frac_{nm}"
    print(f"\n--- {nm}")
    print(f"  {'stratum':28s} {'n':>4s} {'median':>9s} {'IQR':>19s} {'median Mb':>10s}")
    med = {}
    for s in S:
        g = exp_only[exp_only.stratum == s]
        v = g[col].dropna()
        if not len(v):
            print(f"  {s:28s}    0         -"); continue
        med[s] = v.median()
        print(f"  {s:28s} {len(v):4d} {v.median():9.4f} "
              f"{v.quantile(.25):8.4f}-{v.quantile(.75):<10.4f} "
              f"{g.loc[v.index,'size_bp'].median()/1e6:10.2f}")
    groups = [exp_only.loc[exp_only.stratum == s, col].dropna() for s in S]
    groups = [g for g in groups if len(g) >= 5]
    kw = stats.kruskal(*groups) if len(groups) >= 3 else None
    if kw:
        print(f"  Kruskal-Wallis across strata: H = {kw.statistic:.2f}, p = {kw.pvalue:.3g}")
    res = {"kruskal_p": float(kw.pvalue) if kw else None, "medians": {k: float(v) for k, v in med.items()}}
    print(f"\n  pairwise Mann-Whitney (two-sided)")
    pairs = [("S3_recent_restriction", "S1_extreme_obligate"),
             ("S3_recent_restriction", "S4_free_living_breadth"),
             ("S3_recent_restriction", "S2_obligate_intracellular"),
             ("S2_obligate_intracellular", "S4_free_living_breadth"),
             ("S1_extreme_obligate", "S4_free_living_breadth")]
    for a, b in pairs:
        x = exp_only.loc[exp_only.stratum == a, col].dropna()
        y = exp_only.loc[exp_only.stratum == b, col].dropna()
        if len(x) < 5 or len(y) < 5:
            continue
        u = stats.mannwhitneyu(x, y, alternative="two-sided")
        rbc = 2 * u.statistic / (len(x) * len(y)) - 1
        sign = ">" if x.median() > y.median() else "<"
        star = "*" if u.pvalue < ALPHA else " "
        print(f"    {a.split('_')[0]} {sign} {b.split('_')[0]}   "
              f"{x.median():.4f} vs {y.median():.4f}   p = {u.pvalue:9.3g}{star}  rbc {rbc:+.3f}")
        res[f"{a.split('_')[0]}_vs_{b.split('_')[0]}"] = dict(
            p=float(u.pvalue), rbc=float(rbc), higher=a.split('_')[0] if sign == ">" else b.split('_')[0])
    # the declared rule
    k31 = res.get("S3_vs_S1", {}); k34 = res.get("S3_vs_S4", {})
    ok = (k31.get("higher") == "S3" and k31.get("p", 1) < ALPHA
          and k34.get("higher") == "S3" and k34.get("p", 1) < ALPHA)
    res["non_monotonic_confirmed"] = bool(ok)
    print(f"  -> declared rule (S3>S1 and S3>S4, both p<0.05): "
          f"{'MET' if ok else 'NOT MET'}")
    out["per_caller"][nm] = res

n_ok = sum(1 for v in out["per_caller"].values() if v["non_monotonic_confirmed"])
print("\n" + "=" * 92)
print(f"C3 VERDICT: the rule is met by {n_ok} of {len(CALLERS)} callers available "
      f"({', '.join(CALLERS)})")
print(f"  the plan requires at least three of four. "
      f"{'Bakta still incomplete.' if 'Bakta' not in CALLERS else ''}")
out["n_callers_confirming"] = n_ok

print("\n=== context: the original 363, which carry no stratum")
for nm in CALLERS:
    v = d.loc[d.stratum == "original_corpus", f"frac_{nm}"].dropna()
    print(f"  {nm:14s} n={len(v):4d}  median {v.median():.4f}")

d.to_csv(R / "corpus523.tsv", sep="\t", index=False)
json.dump(out, open(R / "c3_strata.json", "w"), indent=1)
print(f"\nwrote results/corpus523.tsv and results/c3_strata.json")
