#!/usr/bin/env python3
"""POST HOC. Not registered in any plan. DFAST is the caller nobody has explained.

The locus run (LOCUS_RESULTS.md) found Pseudofinder and PGAP are nested rather
than opposed, and it found something it did not pursue: DFAST is the sole
pseudogene caller at 2 639 discordant loci, more than any other caller, and it
appears on the intact side of PGAP's 113 sole calls. Section 3.5's clustered
interval for DFAST includes zero, so DFAST is also the caller whose answer about
reduced genomes is the least decided.

Three questions, in order of how much they could change the paper:

  D1  Is DFAST nested with the others, like PGAP inside Pseudofinder, or is it
      calling a different set of genes? The full ordered pair matrix answers
      this; L2's 1955-against-0 was one cell of it.
  D2  Are DFAST's unique calls backed by selection? Same test as the
      Pseudofinder-only post hoc: sole-caller set against the
      concordant-intact baseline, one difference of medians per genome, sign
      test over genomes. Run for all four callers so they are comparable.
  D3  Does DFAST's unique set look different in kind - shorter, longer, nearer
      the neutral line - and does it concentrate in a stratum? That is where a
      caller-specific artefact would show.

Nothing here is a registered prediction and nothing here gets a verdict. The
output is an observation to be reviewed.
"""
import json
import pathlib

import numpy as np
import pandas as pd
from scipy import stats

P = pathlib.Path(__file__).resolve().parent
CALLERS = ["PGAP", "Pseudofinder", "Bakta", "DFAST"]

t = pd.read_csv(P / "results/loci.tsv", sep="\t")
dis = t[t.kind == "discordant"].copy()
con = t[t.kind == "concordant_intact"].copy()
for col in ("pseudo_callers", "intact_callers"):
    dis[col + "_set"] = dis[col].fillna("").apply(
        lambda s: set(x for x in s.split("+") if x))

print(f"{len(t)} loci: {len(dis)} discordant, {len(con)} concordant-intact")
print(f"{dis.accession.nunique()} genomes carry a discordant locus\n")

out = {"n_loci": int(len(t)), "n_discordant": int(len(dis))}


def sign_test_vs_baseline(sub, min_loci=3):
    """One difference of medians per genome, then a sign test over genomes.

    The same procedure as G3 and the Pseudofinder-only post hoc: exact zeros
    are dropped because a sign test requires it, and at least five non-tied
    genomes are needed before a p-value is reported at all.
    """
    g = []
    for acc in set(sub.accession) & set(con.accession):
        x = sub.loc[sub.accession == acc, "omega"]
        y = con.loc[con.accession == acc, "omega"]
        if len(x) >= min_loci and len(y) >= min_loci:
            g.append(float(np.median(x) - np.median(y)))
    nz = [v for v in g if v != 0.0]
    if len(nz) < 5:
        return {"n_genomes": len(g), "p": None,
                "note": f"only {len(nz)} non-tied genomes, too few to test"}
    pos = sum(1 for v in nz if v > 0)
    return {"n_genomes": len(g), "n_tied_dropped": len(g) - len(nz),
            "median_diff": float(np.median(g)), "n_positive": pos,
            "n_tested": len(nz),
            "p": float(stats.binomtest(pos, len(nz), 0.5).pvalue)}


# ---------------------------------------------------------------- D1
print("=== D1  ordered pairs: rows call pseudogene, columns call intact")
print("    (a caller that does not annotate the region is absent from both,")
print("     so every cell is an ACTIVE disagreement, not a silence)\n")
M = pd.DataFrame(0, index=CALLERS, columns=CALLERS, dtype=int)
for _, r in dis.iterrows():
    for a in r.pseudo_callers_set:
        for b in r.intact_callers_set:
            M.loc[a, b] += 1
print(M.to_string())

print("\n    the same cells read as asymmetry, A-pseudo/B-intact vs the reverse")
pairs = {}
for i, a in enumerate(CALLERS):
    for b in CALLERS[i + 1:]:
        fwd, rev = int(M.loc[a, b]), int(M.loc[b, a])
        tot = fwd + rev
        lab = "nested" if min(fwd, rev) == 0 else (
            "near-nested" if min(fwd, rev) / tot < 0.05 else "opposed")
        pairs[f"{a}|{b}"] = {"a_pseudo_b_intact": fwd,
                             "b_pseudo_a_intact": rev, "shape": lab}
        print(f"      {a:13s} vs {b:13s}  {fwd:6d} / {rev:6d}   {lab}")
out["pair_matrix"] = {a: {b: int(M.loc[a, b]) for b in CALLERS} for a in CALLERS}
out["pairs"] = pairs

# ---------------------------------------------------------------- D2
print("\n=== D2  each caller's SOLE pseudogene calls against the")
print("        concordant-intact baseline, clustered by genome")
sole = {}
for c in CALLERS:
    sub = dis[dis.pseudo_callers == c]
    r = sign_test_vs_baseline(sub)
    r["n_loci"] = int(len(sub))
    r["n_genomes_with_loci"] = int(sub.accession.nunique())
    sole[c] = r
    if r.get("p") is None:
        print(f"  {c:13s} {len(sub):6d} loci  -  {r['note']}")
    else:
        print(f"  {c:13s} {len(sub):6d} loci  median {r['median_diff']:+.4f}"
              f"  {r['n_positive']}/{r['n_tested']} positive"
              f"  p {r['p']:.3g}  ({r['n_tied_dropped']} tied dropped)")
out["sole_caller_vs_baseline"] = sole

# ---------------------------------------------------------------- D3
print("\n=== D3  does DFAST's unique set look different in kind?")
rows = []
for c in CALLERS:
    sub = dis[dis.pseudo_callers == c]
    if not len(sub):
        continue
    rows.append({"caller": c, "n": len(sub),
                 "median_omega": float(np.median(sub.omega)),
                 "median_codons": float(np.median(sub.codons)),
                 "frac_omega_above_1": float((sub.omega > 1).mean())})
rows.append({"caller": "(baseline)", "n": len(con),
             "median_omega": float(np.median(con.omega)),
             "median_codons": float(np.median(con.codons)),
             "frac_omega_above_1": float((con.omega > 1).mean())})
d3 = pd.DataFrame(rows)
print(d3.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
out["sole_caller_shape"] = rows

print("\n    DFAST-only loci by stratum, against that stratum's share of loci")
dfo = dis[dis.pseudo_callers == "DFAST"]
by = pd.DataFrame({
    "dfast_only": dfo.groupby("stratum").size(),
    "all_loci": t.groupby("stratum").size(),
}).fillna(0)
by["per_1000_loci"] = 1000 * by.dfast_only / by.all_loci
print(by.astype({"dfast_only": int, "all_loci": int}).to_string(
    float_format=lambda v: f"{v:.1f}"))
out["dfast_only_by_stratum"] = {
    k: {"dfast_only": int(v.dfast_only), "all_loci": int(v.all_loci),
        "per_1000_loci": float(v.per_1000_loci)} for k, v in by.iterrows()}

print("\n    how often each caller is on the INTACT side of someone else's call")
for c in CALLERS:
    n_i = int(dis.intact_callers_set.apply(lambda s: c in s).sum())
    n_p = int(dis.pseudo_callers_set.apply(lambda s: c in s).sum())
    print(f"      {c:13s} pseudogene side {n_p:6d}   intact side {n_i:6d}")
    out.setdefault("sides", {})[c] = {"pseudo_side": n_p, "intact_side": n_i}

pathlib.Path(P / "results/dfast.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("\nwrote results/dfast.json - POST HOC, no verdicts")
