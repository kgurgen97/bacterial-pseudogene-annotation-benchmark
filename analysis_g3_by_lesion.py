#!/usr/bin/env python3
"""Does G3 survive where the frameshift artefact cannot operate?

POST HOC. Not registered. Not a verdict. This is the test FRAME2_RESULTS.md
named as the obvious next step, and it is the crux of whether the paper's dN/dS
corroboration survives.

FRAME2 established two things on simulated lesions:
  a frameshift inflates omega by a median +0.0315, purely mechanically
  a premature stop inflates it by +0.0001

So if G3's gap is produced by the artefact, it should live in the
frameshift-bearing calls and vanish in the internal-stop-only calls. If instead
the internal-stop calls show the gap too, the artefact cannot be the
explanation, because for them there is almost nothing to explain it with.

DFAST is the only caller that publishes a per-call lesion class, in
pseudogene_summary.tsv: frameshift, internal_stop_codon, or both. Its 540
sole-caller loci can therefore be split. PGAP, Pseudofinder and Bakta do not
expose a comparable class, so this test is DFAST-only and that is a real limit
on its reach.

    python3 analysis_g3_by_lesion.py
"""
import ast
import collections
import json
import pathlib
import re
import types

import numpy as np
import pandas as pd
from scipy import stats

P = pathlib.Path(__file__).resolve().parent
R = P / "results"
D = pathlib.Path.home() / "research/datasets"
LOC = re.compile(r"^([^:]+):<?(\d+)\.\.>?(\d+)\(([+-])\)")
DISRUPTED = {"frameshift", "internal_stop_codon",
             "frameshift,internal_stop_codon"}


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


def dfast_classes(acc):
    """-> [(real_contig, start0, end0, classification)] in reference coordinates.

    Resolved with the same contig_key the locus analysis uses, and with the
    same per-caller `seen` bookkeeping, so the coordinates line up with
    loci.tsv rather than merely being close.
    """
    contigs = X.genome_seqs(acc)
    if not contigs:
        return []
    resolve = X.contig_key(contigs)
    reglen = X.gff_region_lengths(D / "n2-dfast" / acc / "genome.gff")
    p = D / "n2-dfast" / acc / "pseudogene_summary.tsv"
    if not p.exists():
        return []
    out, seen = [], {}
    for line in p.read_text(errors="replace").splitlines()[1:]:
        f = line.split("\t")
        if len(f) < 4 or f[3] not in DISRUPTED:
            continue
        m = LOC.match(f[2])
        if not m:
            continue
        cg = m.group(1)
        if cg not in seen:
            seen[cg] = len(seen)
        real, _ = resolve(cg, reglen.get(cg), seen[cg])
        if real is None:
            continue
        out.append((real, int(m.group(2)) - 1, int(m.group(3)), f[3]))
    return out


def g3(sub, baseline, min_loci=3):
    gaps = {}
    for acc in sorted(set(sub.accession)):
        d = sub.loc[sub.accession == acc, "omega"]
        c = baseline.loc[baseline.accession == acc, "omega"]
        if len(d) >= min_loci and len(c) >= min_loci:
            gaps[acc] = float(np.median(d) - np.median(c))
    nz = [v for v in gaps.values() if v != 0.0]
    if len(nz) < 5:
        return {"n_genomes": len(gaps), "p": None, "median_gap": None,
                "note": f"only {len(nz)} non-tied genomes"}
    pos = sum(1 for v in nz if v > 0)
    return {"n_genomes": len(gaps), "n_tested": len(nz),
            "median_gap": float(np.median(list(gaps.values()))),
            "n_positive": pos,
            "p": float(stats.binomtest(pos, len(nz), 0.5).pvalue)}


def main() -> None:
    t = pd.read_csv(R / "loci.tsv", sep="\t")
    base = t[t.kind == "concordant_intact"]
    dfo = t[(t.kind == "discordant") & (t.pseudo_callers == "DFAST")].copy()
    print(f"{len(dfo)} DFAST-only discordant loci to classify")

    # attach DFAST's lesion class by >=50% reciprocal overlap
    byacc = collections.defaultdict(list)
    for acc in sorted(set(dfo.accession)):
        byacc[acc] = dfast_classes(acc)
    lab = []
    for r in dfo.itertuples():
        best, bov = None, 0
        qlen = r.end - r.start
        for cg, s, e, cl in byacc.get(r.accession, []):
            if cg != r.contig:
                continue
            ov = min(r.end, e) - max(r.start, s)
            if ov <= 0:
                continue
            if ov >= 0.5 * qlen and ov >= 0.5 * (e - s) and ov > bov:
                best, bov = cl, ov
        lab.append(best)
    dfo["lesion_class"] = lab
    n_lab = int(dfo.lesion_class.notna().sum())
    print(f"  classified: {n_lab} ({100 * n_lab / len(dfo):.1f}%)")
    print(dfo.lesion_class.value_counts(dropna=False).to_string())

    groups = {
        "frameshift only": dfo[dfo.lesion_class == "frameshift"],
        "internal stop only": dfo[dfo.lesion_class == "internal_stop_codon"],
        "both": dfo[dfo.lesion_class == "frameshift,internal_stop_codon"],
        "any frameshift": dfo[dfo.lesion_class.isin(
            ["frameshift", "frameshift,internal_stop_codon"])],
        "no frameshift": dfo[dfo.lesion_class == "internal_stop_codon"],
        "all DFAST-only": dfo,
    }
    out = {"note": "POST HOC; DFAST-only, the only caller publishing a class",
           "n_dfast_only": int(len(dfo)), "n_classified": n_lab,
           "class_counts": {str(k): int(v) for k, v in
                            dfo.lesion_class.value_counts(dropna=False).items()},
           "groups": {}}
    print("\nG3-style gap against the concordant-intact baseline")
    for nm, sub in groups.items():
        r = g3(sub, base)
        out["groups"][nm] = {**r, "n_loci": int(len(sub))}
        if r.get("p") is None:
            print(f"  {nm:22s} {len(sub):5d} loci  {r['n_genomes']:3d} genomes"
                  f"  NOT TESTABLE ({r['note']})")
        else:
            print(f"  {nm:22s} {len(sub):5d} loci  {r['n_genomes']:3d} genomes"
                  f"  {r['median_gap']:+.4f}  {r['n_positive']}/{r['n_tested']}"
                  f"  p {r['p']:.3g}")

    fs = out["groups"]["any frameshift"]
    ns = out["groups"]["no frameshift"]
    print("\nverdict")
    if ns.get("p") is None:
        print("  NOT TESTABLE. The internal-stop-only set is too small to")
        print("  carry a test, so the frameshift artefact can be neither")
        print("  confirmed nor excluded as G3's explanation by this route.")
        out["verdict"] = "NOT_TESTABLE"
    elif ns["p"] < 0.05 and ns["median_gap"] > 0:
        print("  THE ARTEFACT DOES NOT EXPLAIN IT. Calls with no frameshift,")
        print("  where FRAME2 measured a mechanical inflation of +0.0001, show")
        print(f"  a gap of {ns['median_gap']:+.4f} at p {ns['p']:.3g}. Something")
        print("  other than frameshift misalignment is elevating omega.")
        out["verdict"] = "ARTEFACT_INSUFFICIENT"
    else:
        print("  CONSISTENT WITH THE ARTEFACT. The gap is present in")
        print("  frameshift-bearing calls and absent without a frameshift,")
        print("  which is what a mechanical explanation predicts.")
        out["verdict"] = "CONSISTENT_WITH_ARTEFACT"

    (R / "g3_by_lesion.json").write_text(
        json.dumps(out, indent=2, sort_keys=True) + "\n")
    print("\nwrote results/g3_by_lesion.json")


if __name__ == "__main__":
    main()
