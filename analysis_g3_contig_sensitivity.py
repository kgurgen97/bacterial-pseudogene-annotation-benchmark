#!/usr/bin/env python3
"""Does the corrected G3 depend on the 26% of lookups resolved by contig order?

LOCUS_RESULTS.md records that 230,627 of 886,139 interval lookups took the
order-based fall-back in contig_key:

    if length is not None and len(by_len.get(length, [])) == 1:
        return by_len[length][0], False      # unique length: safe
    if index is not None:
        return order[index], True            # ambiguous: resolved by ORDER

Only Bakta and DFAST reach that branch, because only they are given a region
length; PGAP matches after its version suffix is stripped and Pseudofinder's
names are the reference names. So the exposure is real but bounded.

The decisive observation is that the fall-back cannot pick the wrong contig in
a genome that has only one. `order[index]` is then the only contig there is.
So G3 recomputed on single-replicon genomes alone is free of this hazard
entirely, and if it holds there the fall-back cannot be what produces it.

Reported as a sensitivity analysis. LOCUS_PLAN.md is not touched and G3's
registered verdict is not restated from this; it is a robustness check on the
corrected result.

    python3 analysis_g3_contig_sensitivity.py
"""
import ast
import collections
import json
import pathlib
import types

import numpy as np
import pandas as pd
from scipy import stats

P = pathlib.Path(__file__).resolve().parent
R = P / "results"


def borrow(path):
    tree = ast.parse(path.read_text())
    keep = [n for n in tree.body
            if isinstance(n, (ast.Import, ast.ImportFrom, ast.FunctionDef,
                              ast.ClassDef, ast.Assign, ast.AnnAssign))]
    m = types.ModuleType("b")
    m.__dict__["__file__"] = str(path)
    exec(compile(ast.Module(body=keep, type_ignores=[]), str(path), "exec"),
         m.__dict__)
    return m


X = borrow(P / "analysis_dnds.py")


def g3(sub):
    """G3 exactly as LOCUS_PLAN.md fixes it: one value per genome, sign test,
    exact ties dropped, at least five non-tied genomes."""
    gaps = []
    for acc, g in sub.groupby("accession"):
        d = g.loc[g.kind == "discordant", "omega"]
        c = g.loc[g.kind == "concordant_intact", "omega"]
        if len(d) >= 3 and len(c) >= 3:
            gaps.append(float(np.median(d) - np.median(c)))
    nz = [v for v in gaps if v != 0.0]
    if len(nz) < 5:
        return {"n_genomes": len(gaps), "p": None,
                "note": f"only {len(nz)} non-tied genomes"}
    pos = sum(1 for v in nz if v > 0)
    return {"n_genomes": len(gaps), "n_tied_dropped": len(gaps) - len(nz),
            "median_gap": float(np.median(gaps)), "n_positive": pos,
            "n_tested": len(nz),
            "p": float(stats.binomtest(pos, len(nz), 0.5).pvalue)}


def main() -> None:
    t = pd.read_csv(R / "loci.tsv", sep="\t")
    accs = sorted(t.accession.unique())
    print(f"{len(t)} loci over {len(accs)} genomes\n")

    # how many contigs each genome has, and whether their lengths are unique
    shape = {}
    for acc in accs:
        contigs = X.genome_seqs(acc)
        lens = [len(s) for s in contigs.values()]
        dup = len(lens) - len(set(lens))
        shape[acc] = {"n_contigs": len(contigs), "n_duplicate_lengths": dup}
    n1 = [a for a in accs if shape[a]["n_contigs"] == 1]
    nm = [a for a in accs if shape[a]["n_contigs"] > 1]
    uniq = [a for a in nm if shape[a]["n_duplicate_lengths"] == 0]
    print(f"  single-contig genomes : {len(n1)}  "
          f"(order fall-back cannot choose wrongly)")
    print(f"  multi-contig genomes  : {len(nm)}")
    print(f"    of those, all contig lengths distinct: {len(uniq)}")
    print(f"    with at least two contigs of equal length: "
          f"{len(nm) - len(uniq)}\n")

    arms = {
        "all genomes": t,
        "single-contig only": t[t.accession.isin(n1)],
        "multi-contig only": t[t.accession.isin(nm)],
        "multi-contig, all lengths distinct": t[t.accession.isin(uniq)],
    }
    out = {"shape": shape, "arms": {}}
    for name, sub in arms.items():
        r = g3(sub)
        out["arms"][name] = r
        if r.get("p") is None:
            print(f"  {name:36s} {r['n_genomes']:3d} genomes  "
                  f"NOT TESTABLE: {r['note']}")
        else:
            print(f"  {name:36s} {r['n_genomes']:3d} genomes  "
                  f"{r['median_gap']:+.4f}  {r['n_positive']}/{r['n_tested']} "
                  f"positive  p {r['p']:.3g}")

    a, s = out["arms"]["all genomes"], out["arms"]["single-contig only"]
    print("\nverdict")
    if s.get("p") is None:
        print("  NOT TESTABLE on single-contig genomes alone; the hazard")
        print("  cannot be excluded this way and needs another design.")
    elif s["p"] < 0.05 and s["median_gap"] > 0:
        print(f"  G3 HOLDS where the order fall-back cannot err: "
              f"{s['median_gap']:+.4f} at p {s['p']:.3g} on {s['n_genomes']} "
              f"single-contig genomes,")
        print(f"  against {a['median_gap']:+.4f} at p {a['p']:.3g} on all "
              f"{a['n_genomes']}. The 26% of order-resolved lookups are")
        print("  therefore not what produces G3.")
    else:
        print("  G3 does NOT hold once the hazard is removed. The order")
        print("  fall-back may be carrying the result and this must be")
        print("  reported as a threat to it.")

    (R / "g3_contig_sensitivity.json").write_text(
        json.dumps(out, indent=2, sort_keys=True) + "\n")
    print("\nwrote results/g3_contig_sensitivity.json")


if __name__ == "__main__":
    main()
