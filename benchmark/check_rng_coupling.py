#!/usr/bin/env python3
"""Did a shared RNG make PGAP's result move when only Bakta and DFAST changed?

analysis_dnds.py builds ONE module-level generator:

    RNG = np.random.default_rng(20261028)            # line 39

and consumes it inside the per-genome, per-caller loop:

    if len(it) > len(ps):
        idx = RNG.choice(len(it), len(ps), replace=False)   # line 537

So intact controls are down-sampled to match the pseudogene count, from a
single stateful stream shared by every caller and every genome. The classifier
fix changed len(ps) and len(it) for Bakta and DFAST, which changes how far that
stream advances, which changes every draw that comes after.

Caller order is PGAP, Pseudofinder, Bakta, DFAST, so within one genome PGAP
draws first. But genome k's PGAP draw follows genome k-1's Bakta and DFAST
draws. That gives a sharp prediction:

  PGAP's PSEUDOGENE omegas are identical before and after, in every genome.
      They are not sampled and its parser did not change.
  PGAP's INTACT omegas are identical in the FIRST pair only, and differ in
      every later pair, purely from stream position.

If both hold, PGAP's p moving from 0.015 to 0.097 is RNG coupling and not a
parser effect -- and the published result was never robust to anything that
changed another caller's gene counts.

    python3 benchmark/check_rng_coupling.py
"""
import json
import pathlib

import numpy as np
import pandas as pd
from scipy import stats

P = pathlib.Path(__file__).resolve().parent.parent
PRE = P / "results/pre_fix/dnds_pairs.tsv"
POST = P / "results/dnds_pairs.tsv"


def fingerprint(t, caller, call):
    """Per genome, the sorted omega multiset -- a fingerprint of which genes
    were used, without needing gene ids the table does not carry."""
    s = t[(t.caller == caller) & (t.call == call)]
    return {acc: tuple(np.round(np.sort(g.omega.values), 12))
            for acc, g in s.groupby("accession")}


def gaps(t, caller):
    """One median difference per genome: pseudogene minus intact."""
    out = {}
    s = t[t.caller == caller]
    for acc, g in s.groupby("accession"):
        ps = g.loc[g.call == "pseudo", "omega"]
        it = g.loc[g.call == "intact", "omega"]
        if len(ps) >= 3 and len(it) >= 3:
            out[acc] = float(np.median(ps) - np.median(it))
    return out


def sign_test(vals):
    nz = [v for v in vals if v != 0.0]
    if len(nz) < 5:
        return None, len(vals), len(vals) - len(nz), None
    pos = sum(1 for v in nz if v > 0)
    return (float(stats.binomtest(pos, len(nz), 0.5).pvalue),
            len(vals), len(vals) - len(nz), pos)


def main() -> None:
    if not PRE.exists() or not POST.exists():
        raise SystemExit("need both results/pre_fix/dnds_pairs.tsv and "
                         "results/dnds_pairs.tsv")
    pre = pd.read_csv(PRE, sep="\t")
    post = pd.read_csv(POST, sep="\t")
    # pair order: rows are appended as pairs are processed
    order = list(dict.fromkeys(pre.accession))
    print(f"pre {len(pre)} rows, post {len(post)} rows, "
          f"{len(order)} pairs in order\n")

    report = {}
    for caller in ("PGAP", "Pseudofinder", "Bakta", "DFAST"):
        a_ps, b_ps = fingerprint(pre, caller, "pseudo"), fingerprint(post, caller, "pseudo")
        a_it, b_it = fingerprint(pre, caller, "intact"), fingerprint(post, caller, "intact")
        shared = [x for x in order if x in a_ps and x in b_ps]
        same_ps = [x for x in shared if a_ps[x] == b_ps[x]]
        shared_it = [x for x in order if x in a_it and x in b_it]
        same_it = [x for x in shared_it if a_it[x] == b_it[x]]
        print(f"{caller}")
        print(f"  pseudogene omegas identical in {len(same_ps)}/{len(shared)} genomes")
        print(f"  intact     omegas identical in {len(same_it)}/{len(shared_it)} genomes")
        if same_it:
            pos = [order.index(x) + 1 for x in same_it]
            print(f"    identical at pair positions: {sorted(pos)[:12]}")
        report[caller] = {
            "pseudo_identical": len(same_ps), "pseudo_compared": len(shared),
            "intact_identical": len(same_it), "intact_compared": len(shared_it),
            "intact_identical_positions": sorted(order.index(x) + 1
                                                 for x in same_it)}

    print("\n=== PGAP: the per-genome gap, before and after")
    ga, gb = gaps(pre, "PGAP"), gaps(post, "PGAP")
    common = [x for x in order if x in ga and x in gb]
    pa, na, ta, posa = sign_test([ga[x] for x in common])
    pb, nb, tb, posb = sign_test([gb[x] for x in common])
    print(f"  {len(common)} genomes with both")
    print(f"  before: median gap {np.median([ga[x] for x in common]):+.4f}  "
          f"{posa}/{na - ta} positive  sign test p {pa:.4g}")
    print(f"  after : median gap {np.median([gb[x] for x in common]):+.4f}  "
          f"{posb}/{nb - tb} positive  sign test p {pb:.4g}")
    flipped = [x for x in common if (ga[x] > 0) != (gb[x] > 0)]
    print(f"  genomes whose gap changed SIGN: {len(flipped)}")
    for x in flipped[:12]:
        print(f"    {x}  {ga[x]:+.4f} -> {gb[x]:+.4f}  (pair "
              f"{order.index(x) + 1})")
    moved = sum(1 for x in common if abs(ga[x] - gb[x]) > 1e-12)
    print(f"  genomes whose gap changed at all: {moved}/{len(common)}")

    report["PGAP_gap"] = {
        "n_genomes": len(common),
        "before": {"median": float(np.median([ga[x] for x in common])),
                   "p": pa, "n_positive": posa},
        "after": {"median": float(np.median([gb[x] for x in common])),
                  "p": pb, "n_positive": posb},
        "n_sign_flipped": len(flipped), "n_changed": moved,
        "flipped": [{"accession": x, "before": ga[x], "after": gb[x],
                     "pair_position": order.index(x) + 1} for x in flipped]}

    per = P / "benchmark/results/pgap_sampling_comparison.tsv"
    with per.open("w") as fh:
        fh.write("accession\tpair_position\tgap_before\tgap_after\tdelta\t"
                 "sign_flipped\tintact_omegas_identical\n")
        a_it, b_it = fingerprint(pre, "PGAP", "intact"), fingerprint(post, "PGAP", "intact")
        for x in common:
            same = a_it.get(x) == b_it.get(x)
            fh.write(f"{x}\t{order.index(x) + 1}\t{ga[x]:.6f}\t{gb[x]:.6f}\t"
                     f"{gb[x] - ga[x]:+.6f}\t{(ga[x] > 0) != (gb[x] > 0)}\t"
                     f"{same}\n")
    print(f"\nwrote {per.relative_to(P)}")

    verdict = (report["PGAP"]["pseudo_identical"] == report["PGAP"]["pseudo_compared"]
               and report["PGAP"]["intact_identical"] < report["PGAP"]["intact_compared"])
    report["mechanism_confirmed"] = bool(verdict)
    print("\n=== verdict")
    if verdict:
        print("  CONFIRMED. PGAP's pseudogene omegas are unchanged in every")
        print("  genome, so its parser did not change. Its intact controls")
        print("  changed, which can only come from the shared RNG's position.")
    else:
        print("  NOT confirmed by this test; the pattern does not match a pure")
        print("  stream-position effect and something else also moved.")

    (P / "benchmark/results/rng_coupling.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
