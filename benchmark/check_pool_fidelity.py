#!/usr/bin/env python3
"""Does the pool-plus-replay machinery reproduce the original procedure exactly?

A seed sweep is only interpretable if the one thing that differs between runs
is the seed. Two ways my design could have failed that, and both are tested
here against the recomputed run in results/dnds_pairs.tsv.

F1  omega fidelity. The pool computes omega for every gene in one large
    diamond query; the original computes it for a matched subsample in a
    smaller one. If diamond's best hit for a query depended on which other
    queries were present, the pool's omegas would not be the original's. The
    PSEUDOGENE class is never sampled, so the pool's pseudogene omega multiset
    must equal the original's exactly, per genome and caller. That is a direct
    test with no sampling in the way.

F2  procedure fidelity. The original draws

        if len(it) > len(ps):
            idx = RNG.choice(len(it), len(ps), replace=False)

    from ONE global stream seeded 20261028, in pair order, callers in
    CALLERS order, skipping a caller with no pseudogenes BEFORE drawing. It
    draws from the FULL intact list and only afterwards loses genes that
    produce no omega. Replaying that stream and intersecting the drawn indices
    with the pool must reproduce the original's intact omega multiset exactly.

If F1 and F2 both pass, the sweep differs from the published method only in
which stream the indices come from -- which is the whole point.

    python3 benchmark/check_pool_fidelity.py
"""
import ast
import collections
import csv
import json
import pathlib
import sys
import types

import numpy as np

P = pathlib.Path(__file__).resolve().parent.parent
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

# omega_pool.tsv writes omega with :.8f, so its values are truncated at eight
# decimals while dnds_pairs.tsv carries full float precision. An exact multiset
# comparison therefore fails on rounding alone -- which is what it first did,
# for every cell. Compare sorted values with a tolerance well inside the
# written precision and far below any meaningful omega difference.
ATOL = 1e-6


def ms(vals):
    return np.sort(np.asarray(list(vals), dtype=float))


def equal(a, b):
    return len(a) == len(b) and (len(a) == 0 or
                                 bool(np.allclose(a, b, atol=ATOL, rtol=0)))


def main() -> None:
    pool_f = R / "omega_pool.tsv"
    if not pool_f.exists():
        sys.exit("results/omega_pool.tsv not built yet")

    # the original run's omegas, per genome and caller
    orig = collections.defaultdict(lambda: {"pseudo": [], "intact": []})
    order = []
    with (R / "dnds_pairs.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            orig[(r["accession"], r["caller"])][r["call"]].append(float(r["omega"]))
            if r["accession"] not in order:
                order.append(r["accession"])

    # the pool's omegas, keyed by interval index
    pool = collections.defaultdict(lambda: {"pseudo": {}, "intact": {}})
    with pool_f.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pool[(r["accession"], r["caller"])][r["call"]][int(r["idx"])] = \
                float(r["omega"])
    pooled_accs = {a for a, _ in pool}
    complete = len(pooled_accs) == len(order)
    print(f"{len(order)} pairs in the original's order, "
          f"{len(pooled_accs)} present in the pool")
    if not complete:
        # An earlier version reported fidelity_confirmed on a partial pool,
        # which would have licensed a sweep over genomes nobody had checked.
        # Coverage is now part of the verdict.
        print(f"  POOL INCOMPLETE ({len(pooled_accs)}/{len(order)}). The "
              f"overlap is still checked below and is useful as an early read,"
              f"\n  but fidelity is NOT confirmed until every pair is covered.")

    # ---------------------------------------------------------------- F1
    print("\nF1  pseudogene omegas, pool against original (no sampling involved)")
    same = diff = 0
    examples = []
    for (acc, caller), d in sorted(orig.items()):
        if acc not in pooled_accs or not d["pseudo"]:
            continue
        a = ms(d["pseudo"])
        b = ms(pool[(acc, caller)]["pseudo"].values())
        if equal(a, b):
            same += 1
        else:
            diff += 1
            if len(examples) < 5:
                worst = (float(np.max(np.abs(a - b))) if len(a) == len(b)
                         else float("nan"))
                examples.append(f"{acc} {caller}: original {len(a)} values, "
                                f"pool {len(b)}, max |diff| {worst:.3g}")
    print(f"    identical: {same}   different: {diff}")
    for e in examples:
        print(f"      {e}")

    # ---------------------------------------------------------------- F2
    print("\nF2  replaying the single global stream, seed 20261028, in pair order")
    RNG = np.random.default_rng(20261028)
    drawn = {}
    for acc in order:
        for caller in X.CALLERS:
            ps = X.calls_for(acc, caller, True)
            it = X.calls_for(acc, caller, False)
            if not ps:
                continue                      # no draw, exactly as the original
            if len(it) > len(ps):
                idx = RNG.choice(len(it), len(ps), replace=False)
                drawn[(acc, caller)] = set(int(i) for i in idx)
            else:
                drawn[(acc, caller)] = set(range(len(it)))
    print(f"    replayed {len(drawn)} (genome, caller) draws")

    same2 = diff2 = skipped = 0
    examples2 = []
    for (acc, caller), d in sorted(orig.items()):
        if acc not in pooled_accs or not d["intact"]:
            continue
        sel = drawn.get((acc, caller))
        if sel is None:
            skipped += 1
            continue
        have = pool[(acc, caller)]["intact"]
        rebuilt = [v for i, v in have.items() if i in sel]
        a, b = ms(d["intact"]), ms(rebuilt)
        if equal(a, b):
            same2 += 1
        else:
            diff2 += 1
            if len(examples2) < 6:
                worst = (float(np.max(np.abs(a - b))) if len(a) == len(b)
                         else float("nan"))
                examples2.append(
                    f"{acc} {caller}: original {len(a)} values, "
                    f"replay {len(b)}, max |diff| {worst:.3g}")
    print(f"    identical: {same2}   different: {diff2}   no draw recorded: {skipped}")
    for e in examples2:
        print(f"      {e}")

    ok = (diff == 0 and diff2 == 0 and same > 0 and same2 > 0 and complete)
    print("\nverdict")
    if not complete:
        print(f"    NOT CONFIRMED: the pool covers {len(pooled_accs)} of "
              f"{len(order)} pairs.")
        print(f"    On the covered part: F1 {same} identical / {diff} different,"
              f" F2 {same2} / {diff2}.")
    elif ok:
        print("    FIDELITY CONFIRMED. The pool's omegas match the original's")
        print("    where no sampling intervenes, and replaying the global")
        print("    stream reproduces the original's control sets exactly. The")
        print("    sweep therefore varies the seed and nothing else.")
    else:
        print("    FIDELITY NOT CONFIRMED. The sweep's variance would not be")
        print("    comparable to the original's and must not be reported as")
        print("    seed sensitivity until this is resolved.")

    (P / "benchmark/results/pool_fidelity.json").write_text(json.dumps(
        {"pairs_in_order": len(order), "pairs_in_pool": len(pooled_accs),
         "pool_complete": bool(complete),
         "F1_pseudo_identical": same, "F1_pseudo_different": diff,
         "F2_intact_identical": same2, "F2_intact_different": diff2,
         "F2_no_draw": skipped, "fidelity_confirmed": bool(ok),
         "F1_examples": examples, "F2_examples": examples2},
        indent=2, sort_keys=True) + "\n")
    if not ok:
        sys.exit("fidelity not confirmed")


if __name__ == "__main__":
    main()
