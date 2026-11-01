#!/usr/bin/env python3
"""Does the 50 % reciprocal-overlap threshold decide G3? Isolated rebuild.

POST HOC. Not registered. LOCUS_PLAN.md is untouched, and nothing here writes
results/loci.tsv, results/locus.json or any other validated output: everything
goes under results/overlap_sensitivity/.

This is the one sensitivity check FINAL_LOCUS_SENSITIVITY.md recorded as
outstanding. The obstacle is cost: the full locus analysis took 4.7 hours
because it computes ω for 103 224 loci, and five thresholds would be a day.

The way round it is that the threshold changes the PARTITION, not the ω of a
given span. So this rebuilds only the partition — reading the 886 139 caller
intervals, resolving contigs, and running analysis_locus.py's own union-find
with RECIPROCAL overridden — then matches each resulting locus against the
validated loci.tsv by (accession, contig, start, end) and reuses the ω already
computed for it.

What that buys and what it cannot: G3 is recomputed exactly on every locus that
survives the change of threshold, and the number of loci that appear or vanish
is reported. Loci that are NEW at a threshold have no ω and are excluded, so the
G3 figures below are for the stable subset. That is a real sensitivity result,
and its limitation is stated rather than hidden.

    python3 analysis_overlap_sensitivity.py [--probe N] [--thresholds 0.3,0.5,0.7]
"""
from __future__ import annotations

import argparse
import ast
import collections
import csv
import json
import pathlib
import time
import types

import numpy as np
import pandas as pd
from scipy import stats

P = pathlib.Path(__file__).resolve().parent
R = P / "results"
OUT = R / "overlap_sensitivity"
D = pathlib.Path.home() / "research/datasets"
MIN_LOCI = 3


def borrow(path):
    tree = ast.parse(path.read_text())
    keep = [n for n in tree.body
            if isinstance(n, (ast.Import, ast.ImportFrom, ast.FunctionDef,
                              ast.ClassDef, ast.Assign, ast.AnnAssign))]
    m = types.ModuleType(path.stem + "_borrowed")
    m.__dict__["__file__"] = str(path)
    exec(compile(ast.Module(body=keep, type_ignores=[]), str(path), "exec"),
         m.__dict__)
    return m


L = borrow(P / "analysis_locus.py")     # brings its own X = analysis_dnds
X = L.X
CALLERS = L.CALLERS


def build_partition(pairs, recip, probe=None):
    """Locus spans and caller verdicts at one threshold. No omega computed."""
    L.RECIPROCAL = recip                      # analysis_locus.components reads this
    rows = []
    stats_ = collections.Counter()
    for k, (acc, partner, stratum, gen) in enumerate(pairs, 1):
        contigs = X.genome_seqs(acc)
        if not contigs:
            continue
        resolve = X.contig_key(contigs)
        reglen = {
            "Bakta": X.gff_region_lengths(D / "n2-bakta" / acc / f"{acc}.gff3"),
            "DFAST": X.gff_region_lengths(D / "n2-dfast" / acc / "genome.gff"),
        }
        ivs = []
        for caller in CALLERS:
            seen = {}
            for klass, want in (("pseudo", True), ("intact", False)):
                for (c, s, e, st) in X.calls_for(acc, caller, want):
                    if c not in seen:
                        seen[c] = len(seen)
                    real, _ = resolve(c, reglen.get(caller, {}).get(c), seen[c])
                    if real is None or e > len(contigs[real]) or e - s < 30:
                        continue
                    ivs.append((real, s, e, st, caller, klass))
        if not ivs:
            continue
        stats_["intervals"] += len(ivs)
        comps = L.components(ivs)
        stats_["unjoined"] += sum(1 for c in comps if len(c) == 1)
        for comp in comps:
            if len(comp) < 2:
                continue
            verdict = {}
            for i in comp:
                caller, klass = ivs[i][4], ivs[i][5]
                if caller in verdict and verdict[caller] != klass:
                    verdict[caller] = "ambiguous"
                elif caller not in verdict:
                    verdict[caller] = klass
            if "ambiguous" in verdict.values():
                stats_["ambiguous"] += 1
                continue
            ps = sorted(c for c, v in verdict.items() if v == "pseudo")
            it = sorted(c for c, v in verdict.items() if v == "intact")
            if ps and it:
                kind = "discordant"
            elif it and not ps and len(it) >= 2:
                kind = "concordant_intact"
            else:
                continue
            # analysis_locus.py line 206 takes the LONGEST interval in the
            # component as the locus's representative and records its own
            # coordinates, not the component's min/max span. Reproducing that
            # exactly is what lets a locus be matched back to the validated
            # table; an earlier version used min/max and matched only 46 %.
            longest = max(comp, key=lambda i: ivs[i][2] - ivs[i][1])
            cg, s, e = ivs[longest][0], ivs[longest][1], ivs[longest][2]
            rows.append({"accession": acc, "contig": cg, "start": s, "end": e,
                         "kind": kind,
                         "pseudo_callers": "+".join(ps) or "-",
                         "intact_callers": "+".join(it) or "-"})
        if k % 20 == 0:
            print(f"    {k}/{len(pairs)} pairs, {len(rows)} loci", flush=True)
    return rows, stats_


def g3(t, min_loci=MIN_LOCI):
    d = t[t.kind == "discordant"]
    c = t[t.kind == "concordant_intact"]
    gaps = {}
    for acc in sorted(set(d.accession) & set(c.accession)):
        a = d.loc[d.accession == acc, "omega"]
        b = c.loc[c.accession == acc, "omega"]
        if len(a) >= min_loci and len(b) >= min_loci:
            gaps[acc] = float(np.median(a) - np.median(b))
    nz = [v for v in gaps.values() if v != 0.0]
    if len(nz) < 5:
        return {"n_genomes": len(gaps), "p": None, "median_gap": None}
    pos = sum(1 for v in nz if v > 0)
    return {"n_genomes": len(gaps), "n_tested": len(nz),
            "median_gap": float(np.median(list(gaps.values()))),
            "n_positive": pos,
            "p": float(stats.binomtest(pos, len(nz), 0.5).pvalue)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", type=int, default=0)
    ap.add_argument("--thresholds", default="0.30,0.40,0.50,0.60,0.70")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    base = pd.read_csv(R / "loci.tsv", sep="\t")
    omega = {(r.accession, r.contig, r.start, r.end): r.omega
             for r in base.itertuples()}
    print(f"validated base: {len(base)} loci with omega, "
          f"{int((base.kind == 'discordant').sum())} discordant")
    base_g3 = g3(base)
    print(f"base G3: {base_g3['median_gap']:+.4f} "
          f"{base_g3['n_positive']}/{base_g3['n_tested']} p {base_g3['p']:.3g}")

    pairs = L.build_pairs(probe=a.probe) if hasattr(L, "build_pairs") else None
    if pairs is None:
        # the pair list lives inside analysis_locus.main(); rebuild it from the
        # validated table so the same pairs in the same order are used
        seen = {}
        for r in base.itertuples():
            seen.setdefault(r.accession, (r.genus, r.stratum))
        pairs = [(acc, "", st, gn) for acc, (gn, st) in seen.items()]
    if a.probe:
        pairs = pairs[:a.probe]
    print(f"pairs: {len(pairs)}")

    out = {"note": "POST HOC; isolated outputs; omega reused from the "
                   "validated loci.tsv, so G3 is on the stable subset",
           "base": base_g3, "n_base_loci": int(len(base)), "thresholds": {}}

    for ts in a.thresholds.split(","):
        recip = float(ts)
        t0 = time.time()
        print(f"\n=== reciprocal overlap {recip:.2f}")
        rows, st = build_partition(pairs, recip)
        dt = time.time() - t0
        df = pd.DataFrame(rows)
        key = list(zip(df.accession, df.contig, df.start, df.end))
        df["omega"] = [omega.get(k, np.nan) for k in key]
        matched = df[df.omega.notna()].copy()
        gmap = dict(zip(base.accession, base.genus))
        matched["genus"] = matched.accession.map(gmap)
        r = g3(matched)
        n_new = int(df.omega.isna().sum())
        # compare against the base restricted to the SAME genomes, or a probe
        # over 6 of 99 genomes would look like it had lost 93 % of the base
        accs = set(df.accession)
        base_here = base[base.accession.isin(accs)]
        n_lost = len(base_here) - len(matched)
        out["thresholds"][ts] = {
            "n_loci": int(len(df)),
            "n_discordant": int((df.kind == "discordant").sum()),
            "n_concordant": int((df.kind == "concordant_intact").sum()),
            "n_matched_to_base": int(len(matched)),
            "n_new_without_omega": n_new,
            "n_base_loci_in_these_genomes": int(len(base_here)),
            "n_base_loci_not_present": n_lost,
            "frac_base_preserved": len(matched) / max(len(base_here), 1),
            "intervals": st["intervals"], "unjoined": st["unjoined"],
            "ambiguous": st["ambiguous"],
            "g3_on_matched": r, "seconds": round(dt, 1)}
        p = f"{r['p']:.3g}" if r.get("p") else "n/a"
        mg = f"{r['median_gap']:+.4f}" if r.get("median_gap") else "n/a"
        print(f"  loci {len(df)} ({(df.kind=='discordant').sum()} discordant)"
              f"  matched {len(matched)} of {len(base_here)} base loci here "
              f"({100*len(matched)/max(len(base_here),1):.1f}%)  new {n_new}")
        print(f"  G3 on matched: {mg}  {r.get('n_positive')}/{r.get('n_tested')}"
              f"  p {p}   [{dt:.0f}s]")
        matched.to_csv(OUT / f"loci_recip_{ts}.tsv", sep="\t", index=False)

    (R / "overlap_sensitivity.json").write_text(
        json.dumps(out, indent=2, sort_keys=True) + "\n")
    print("\nwrote results/overlap_sensitivity.json and "
          "results/overlap_sensitivity/*.tsv")
    print("results/loci.tsv and results/locus.json were NOT touched")


if __name__ == "__main__":
    main()
