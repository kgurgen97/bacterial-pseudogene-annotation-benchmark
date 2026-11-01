#!/usr/bin/env python3
"""Does which genes were SAMPLED change which contig a gene is read from?

analysis_dnds.py resolves a caller's contig name to the assembly's by unique
length, falling back to first-appearance order when the length is ambiguous.
LOCUS_RESULTS.md records that 230,627 of 886,139 lookups -- 26% -- took that
fall-back.

The first-appearance index comes from this, inside the per-caller loop:

    seen = {}
    for tag, lst in (("pseudo", ps), ("intact", it)):      # it is ALREADY SAMPLED
        for i, (c, s, e, st) in enumerate(lst):
            if c not in seen:
                seen[c] = len(seen)

`it` at that point is the RNG-sampled subset. So the order in which contigs are
first encountered -- and therefore the fall-back's answer -- depends on which
intact genes the RNG happened to draw. If so, the published omega of a
PSEUDOGENE can depend on the control sample, because the gene's sequence is
read from a different contig.

That is a third reproducibility hazard, compounding the shared RNG. This script
tests it directly: for each genome and caller, resolve every pseudogene
interval twice -- once with `seen` built from ps plus the sampled intact list,
once from ps plus the full intact list -- and count the disagreements.

    python3 benchmark/check_contig_order_dependence.py [n_genomes]
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


def resolve_all(resolve, reglen, caller, ps, it):
    """Replicate the loop's `seen` bookkeeping and return pseudogene contigs."""
    seen, out, by_order = {}, [], 0
    for tag, lst in (("pseudo", ps), ("intact", it)):
        for i, (c, s, e, st) in enumerate(lst):
            if c not in seen:
                seen[c] = len(seen)
            real, order = resolve(c, reglen.get(c), seen[c])
            if tag == "pseudo":
                out.append(real)
                by_order += bool(order)
    return out, by_order


def main() -> None:
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    order = []
    with (R / "dnds_pairs.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["accession"] not in order:
                order.append(r["accession"])
    accs = order[:n]
    print(f"{len(accs)} genomes, replaying the original stream for the sample\n")

    # the original stream, so the "sampled" arm is the real sampled arm
    RNG = np.random.default_rng(20261028)
    sampled = {}
    for acc in order:
        for caller in X.CALLERS:
            ps = X.calls_for(acc, caller, True)
            it = X.calls_for(acc, caller, False)
            if not ps:
                continue
            if len(it) > len(ps):
                idx = RNG.choice(len(it), len(ps), replace=False)
                sampled[(acc, caller)] = [it[i] for i in idx]
            else:
                sampled[(acc, caller)] = it

    tot = diff = 0
    fallback_s = fallback_f = 0
    per = []
    for acc in accs:
        contigs = X.genome_seqs(acc)
        if not contigs:
            continue
        resolve = X.contig_key(contigs)
        reglen = {
            "Bakta": X.gff_region_lengths(X.D / "n2-bakta" / acc / f"{acc}.gff3"),
            "DFAST": X.gff_region_lengths(X.D / "n2-dfast" / acc / "genome.gff"),
        }
        for caller in X.CALLERS:
            ps = X.calls_for(acc, caller, True)
            it_full = X.calls_for(acc, caller, False)
            it_s = sampled.get((acc, caller))
            if not ps or it_s is None:
                continue
            rl = reglen.get(caller, {})
            a, fa = resolve_all(resolve, rl, caller, ps, it_s)
            b, fb = resolve_all(resolve, rl, caller, ps, it_full)
            d = sum(1 for x, y in zip(a, b) if x != y)
            tot += len(a)
            diff += d
            fallback_s += fa
            fallback_f += fb
            if d:
                per.append({"accession": acc, "caller": caller,
                            "n_pseudo": len(a), "n_differing": d})
    print(f"pseudogene intervals resolved: {tot}")
    print(f"  resolved to a DIFFERENT contig depending on the sample: {diff}"
          f"  ({100 * diff / max(tot, 1):.1f}%)")
    print(f"  order fall-back used, sampled arm: {fallback_s}")
    print(f"  order fall-back used, full arm   : {fallback_f}")
    if per:
        print("\n  worst affected (genome, caller):")
        for r in sorted(per, key=lambda r: -r["n_differing"])[:8]:
            print(f"    {r['accession']:18s} {r['caller']:13s} "
                  f"{r['n_differing']:5d} of {r['n_pseudo']}")
    print("\nverdict")
    if diff:
        print("  CONFIRMED. Which contig a pseudogene's sequence is read from")
        print("  depends on which intact controls the RNG drew, through the")
        print("  order-based contig fall-back. So a pseudogene's omega is not")
        print("  a property of that gene alone, and the omega pool cannot")
        print("  reproduce the original run's values by construction.")
    else:
        print("  not reproduced on these genomes: contig resolution did not")
        print("  depend on the sample here.")
    (P / "benchmark/results/contig_order_dependence.json").write_text(json.dumps(
        {"n_genomes": len(accs), "n_pseudo_resolved": tot,
         "n_differing": diff,
         "frac_differing": diff / max(tot, 1),
         "fallback_sampled": fallback_s, "fallback_full": fallback_f,
         "worst": sorted(per, key=lambda r: -r["n_differing"])[:20]},
        indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
