#!/usr/bin/env python3
"""Which of Bakta's pseudogene intervals come from which attribute key?

Bakta asserts pseudogenes two ways: `pseudogene=unitary|unknown` on 6,844 CDS
rows corpus-wide, and `pseudo=True` on a further 509. Its own `pseudogenes:`
summary line counts only the first, so the count tables use 6,844 while the
corrected interval analyses use all 7,353 (SPREAD_METRICS.md).

To test whether that choice changes any headline claim, the two definitions
have to be separable *inside* the interval analyses. calls_for returns a plain
list, and the omega pool records each interval by its index in that list, so
this writes the index -> key mapping. Nothing is recomputed: the same intervals
are simply labelled.

    python3 benchmark/map_bakta_definition.py
"""
import ast
import collections
import csv
import json
import pathlib
import re
import types

P = pathlib.Path(__file__).resolve().parent.parent
D = pathlib.Path.home() / "research/datasets"

STRICT = re.compile(r"(?:^|;)\s*pseudogene\s*=", re.I)
LOOSE = re.compile(r"(?:^|;)\s*pseudo\s*=", re.I)


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


def classify(acc):
    """-> {index in calls_for(...,True): 'strict' | 'loose'}

    Rebuilt by walking the GFF in the same order calls_for does, so the index
    agrees by construction rather than by assumption.
    """
    g = D / f"n2-bakta/{acc}/{acc}.gff3"
    out = {}
    if not g.exists():
        return out
    i = 0
    for line in g.read_text(errors="replace").splitlines():
        if line.startswith("##FASTA"):
            break
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] not in ("CDS", "gene", "pseudogene"):
            continue
        try:
            int(f[3].lstrip("<>")), int(f[4].lstrip("<>"))
        except ValueError:
            continue
        if not X.gff_asserts_pseudo(f):
            continue
        out[i] = "strict" if STRICT.search(f[8]) else "loose"
        i += 1
    return out


def main() -> None:
    rows = list(csv.DictReader((P / "results/dnds_pairs.tsv").open(),
                               delimiter="\t"))
    accs = list(dict.fromkeys(r["accession"] for r in rows))
    print(f"{len(accs)} accessions used by the interval analyses")

    out, tally = [], collections.Counter()
    per_genome = []
    for acc in accs:
        m = classify(acc)
        # the mapping must have the same length as what calls_for returns, or
        # the index is not the index the pool recorded
        n_calls = len(X.calls_for(acc, "Bakta", True))
        if len(m) != n_calls:
            raise SystemExit(f"{acc}: mapping has {len(m)} entries but "
                             f"calls_for returns {n_calls}; the index would "
                             f"not line up with the omega pool")
        s = sum(1 for v in m.values() if v == "strict")
        l = len(m) - s
        tally["strict"] += s
        tally["loose"] += l
        per_genome.append({"accession": acc, "strict": s, "loose": l,
                           "total": len(m)})
        for i, v in m.items():
            out.append({"accession": acc, "idx": i, "definition": v})

    with (P / "results/bakta_definition.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t",
                           fieldnames=["accession", "idx", "definition"])
        w.writeheader()
        w.writerows(out)
    with (P / "benchmark/results/bakta_definition_by_genome.tsv").open(
            "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t",
                           fieldnames=["accession", "strict", "loose", "total"])
        w.writeheader()
        w.writerows(per_genome)

    n_loose_genomes = sum(1 for r in per_genome if r["loose"])
    print(f"\nacross these {len(accs)} genomes:")
    print(f"  pseudogene=  (strict)   {tally['strict']}")
    print(f"  pseudo=True  (extra)    {tally['loose']}")
    print(f"  genomes carrying at least one pseudo=True call: "
          f"{n_loose_genomes}/{len(accs)}")
    if tally["strict"] + tally["loose"]:
        print(f"  the inclusive definition is "
              f"{100 * tally['loose'] / (tally['strict'] + tally['loose']):.1f}% "
              f"larger than the strict one here")
    worst = sorted(per_genome, key=lambda r: -r["loose"])[:5]
    print("  genomes most exposed to the choice:")
    for r in worst:
        print(f"    {r['accession']}  strict {r['strict']:4d}  "
              f"pseudo=True {r['loose']:3d}")

    (P / "benchmark/results/bakta_definition.json").write_text(json.dumps(
        {"n_genomes": len(accs), "strict": tally["strict"],
         "loose": tally["loose"],
         "genomes_with_loose": n_loose_genomes}, indent=2, sort_keys=True) + "\n")
    print("\nwrote results/bakta_definition.tsv and "
          "benchmark/results/bakta_definition_by_genome.tsv")


if __name__ == "__main__":
    main()
