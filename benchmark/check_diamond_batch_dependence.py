#!/usr/bin/env python3
"""Does DIAMOND's reported best hit depend on which OTHER queries are present?

The contig fall-back was ruled out (check_contig_order_dependence.py: 0 of
17,164 pseudogene intervals resolve differently). That leaves best_hits:

    diamond blastx ... --max-target-seqs 1 --evalue 1e-5
    ...
    if q not in best: best[q] = (s, qstart, qend)     # first line per query wins

then a post-filter on pident and coverage. Two reasons that could be
batch-dependent. DIAMOND partitions queries into blocks and its seed-filtering
heuristics are not guaranteed invariant across different partitions, so the
single reported target for a borderline query can change. And because the
pident/coverage filter is applied AFTER `-k 1`, a query whose one reported
target fails the filter is dropped entirely, even when a slightly lower-scoring
target would have passed -- so the SET of queries with hits can move too.

This takes one genome and runs best_hits three times:
  A  the exact query set the original built (pseudogenes + the sampled controls)
  B  the same query set again, to establish determinism at fixed input
  C  the pool's query set (pseudogenes + ALL controls)

and compares what each returns for the PSEUDOGENE queries, which are identical
in all three.

    python3 benchmark/check_diamond_batch_dependence.py [accession]
"""
import ast
import json
import pathlib
import sys
import tempfile
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


def write_queries(path, contigs, resolve, reglen, acc, use_full, sampled):
    """Rebuild the loop's query FASTA, either as the original or as the pool."""
    meta = {}
    with path.open("w") as fh:
        for caller in X.CALLERS:
            ps = X.calls_for(acc, caller, True)
            it = X.calls_for(acc, caller, False)
            if not ps:
                continue
            it_use = it if use_full else sampled[(acc, caller)]
            seen = {}
            for tag, lst in (("pseudo", ps), ("intact", it_use)):
                for i, (c, s, e, st) in enumerate(lst):
                    if c not in seen:
                        seen[c] = len(seen)
                    real, _ = resolve(c, reglen.get(caller, {}).get(c), seen[c])
                    if real is None:
                        continue
                    if e - s < X.MIN_CODONS * 3 or e > len(contigs[real]):
                        continue
                    name = f"{caller}|{tag}|{i}"
                    fh.write(f">{name}\n{X.sub(contigs[real], s, e, st)}\n")
                    meta[name] = True
    return meta


def main() -> None:
    order, partner = [], {}
    with (R / "dnds_pairs.tsv").open() as fh:
        import csv
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["accession"] not in order:
                order.append(r["accession"])
                partner[r["accession"]] = r["partner"]
    acc = sys.argv[1] if len(sys.argv) > 1 else order[0]
    print(f"genome {acc}, partner {partner[acc]}")

    # replay the original stream up to and including this genome
    RNG = np.random.default_rng(20261028)
    sampled = {}
    for a in order:
        for caller in X.CALLERS:
            ps = X.calls_for(a, caller, True)
            it = X.calls_for(a, caller, False)
            if not ps:
                continue
            if len(it) > len(ps):
                idx = RNG.choice(len(it), len(ps), replace=False)
                sampled[(a, caller)] = [it[i] for i in idx]
            else:
                sampled[(a, caller)] = it
        if a == acc:
            break

    contigs = X.genome_seqs(acc)
    resolve = X.contig_key(contigs)
    reglen = {
        "Bakta": X.gff_region_lengths(X.D / "n2-bakta" / acc / f"{acc}.gff3"),
        "DFAST": X.gff_region_lengths(X.D / "n2-dfast" / acc / "genome.gff"),
    }
    pfaa = X.D / f"n2-pseudofinder/{partner[acc]}/pf_intact.faa"

    runs = {}
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        for label, full in (("A original set", False), ("B original set again", False),
                            ("C pool set", True)):
            qf = td / f"{label.split()[0]}.fna"
            meta = write_queries(qf, contigs, resolve, reglen, acc, full, sampled)
            hits = X.best_hits(qf, pfaa)
            runs[label] = hits
            n_ps = sum(1 for k in meta if "|pseudo|" in k)
            print(f"  {label:22s} queries {len(meta):6d} "
                  f"(pseudo {n_ps})  hits {len(hits)}")

    a, b, c = (runs["A original set"], runs["B original set again"],
               runs["C pool set"])
    pseudo = sorted(k for k in a if "|pseudo|" in k)
    print(f"\n{len(pseudo)} pseudogene queries with a hit in run A")

    def cmp(x, y, lx, ly):
        both = [k for k in pseudo if k in x and k in y]
        same = [k for k in both if x[k] == y[k]]
        only_x = [k for k in pseudo if k in x and k not in y]
        print(f"  {lx} vs {ly}: both {len(both)}, identical {len(same)}, "
              f"in {lx} only {len(only_x)}")
        return len(both), len(same), len(only_x)

    det = cmp(a, b, "A", "B")
    bat = cmp(a, c, "A", "C")

    print("\nverdict")
    if det[0] and det[1] == det[0] and det[2] == 0:
        print("  DIAMOND is deterministic at fixed input (A == B).")
    else:
        print("  DIAMOND is NOT deterministic even at fixed input.")
    if bat[1] < bat[0] or bat[2]:
        print("  But its answer DOES depend on which other queries are present")
        print("  (A != C), so the omega of a given gene is not reproducible")
        print("  across query sets. The pool cannot reproduce the original run's")
        print("  values, and a seed sweep must therefore be computed entirely")
        print("  within one query set.")
    else:
        print("  And its answer does not depend on the other queries present.")

    (P / "benchmark/results/diamond_batch_dependence.json").write_text(
        json.dumps({"accession": acc, "partner": partner[acc],
                    "n_pseudo_queries": len(pseudo),
                    "A_vs_B": {"both": det[0], "identical": det[1],
                               "only_A": det[2]},
                    "A_vs_C": {"both": bat[0], "identical": bat[1],
                               "only_A": bat[2]},
                    "deterministic_at_fixed_input": bool(
                        det[0] and det[1] == det[0] and det[2] == 0),
                    "depends_on_query_set": bool(bat[1] < bat[0] or bat[2])},
                   indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
