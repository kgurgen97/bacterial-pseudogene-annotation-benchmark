#!/usr/bin/env python3
"""Does the corrected classifier agree with each caller's own summary counts?

The fix has an independent oracle. run_dfast_all.sh and run_bakta_all.sh
already counted pseudogenes per genome, with awk, from each caller's own
summary file -- and their headers warn against grepping the GFF for "pseudo",
which is the mistake the interval path made. Those counts are in
results/dfast_counts.tsv and results/bakta_counts.tsv.

So the corrected interval functions must reproduce them:

    len(calls_for(acc, "DFAST", True)) == dfast_counts.dfast_disrupted
    len(calls_for(acc, "Bakta", True)) == bakta_counts.bakta_pseudo

Two code paths written months apart, by different means, from different files.
If they agree across the corpus the fix is right; where they disagree the
difference is reported rather than tolerated.

    python3 benchmark/check_classifier_fix.py [n_genomes]
"""
import ast
import collections
import json
import pathlib
import re
import sys
import types

P = pathlib.Path(__file__).resolve().parent.parent


def borrow(path):
    """Load analysis_dnds.py's definitions without running its top level."""
    tree = ast.parse(path.read_text())
    keep = [n for n in tree.body
            if isinstance(n, (ast.Import, ast.ImportFrom, ast.FunctionDef,
                              ast.ClassDef, ast.Assign, ast.AnnAssign))]
    mod = types.ModuleType("dnds_borrowed")
    mod.__dict__["__file__"] = str(path)
    exec(compile(ast.Module(body=keep, type_ignores=[]), str(path), "exec"),
         mod.__dict__)
    return mod


X = borrow(P / "analysis_dnds.py")


def bakta_decompose(acc):
    """Split Bakta's asserted pseudogenes into the two keys it actually uses.

    Bakta writes `pseudogene=unitary|unknown` and, separately, `pseudo=True`.
    Its own `pseudogenes:` summary line counts ONLY the `pseudogene=` rows on
    CDS features -- it excludes `pseudo=True` and excludes the tRNA rows that
    carry `pseudogene=unknown`. Both keys are genuine assertions, so the
    interval set uses both; this decomposition is what makes the superset
    auditable against the summary rather than merely close to it.
    """
    g = pathlib.Path.home() / f"research/datasets/n2-bakta/{acc}/{acc}.gff3"
    strict = extra = 0
    if not g.exists():
        return strict, extra
    for line in g.read_text(errors="replace").splitlines():
        if line.startswith("##FASTA"):
            break
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] not in ("CDS", "gene", "pseudogene"):
            continue
        a = f[8]
        if re.search(r"(?:^|;)\s*pseudogene\s*=", a, re.I):
            strict += 1
        elif re.search(r"(?:^|;)\s*pseudo\s*=", a, re.I):
            extra += 1
    return strict, extra


def load_counts(name, key):
    out = {}
    p = P / "results" / name
    if not p.exists():
        return out
    lines = p.read_text().splitlines()
    head = lines[0].split("\t")
    i = head.index(key)
    for line in lines[1:]:
        f = line.split("\t")
        if len(f) > i:
            try:
                out[f[0]] = int(f[i])
            except ValueError:
                pass
    return out


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 523
    dfast = load_counts("dfast_counts.tsv", "dfast_disrupted")
    bakta = load_counts("bakta_counts.tsv", "bakta_pseudo")
    print(f"oracles: {len(dfast)} DFAST rows, {len(bakta)} Bakta rows")

    rows, bad = [], collections.Counter()
    extra_total = 0
    accs = sorted(set(dfast) | set(bakta))[:n]
    for acc in accs:
        got_d = len(X.calls_for(acc, "DFAST", True))
        got_b = len(X.calls_for(acc, "Bakta", True))
        strict_b, extra_b = bakta_decompose(acc)
        extra_total += extra_b
        want_d, want_b = dfast.get(acc), bakta.get(acc)

        okd = want_d is None or got_d == want_d
        # Bakta's summary counts only the pseudogene= rows, so THAT is what
        # must match it exactly; the interval set is that plus the pseudo=True
        # rows, and the sum must account for every interval.
        okb = (want_b is None or strict_b == want_b) and got_b == strict_b + extra_b
        if not okd:
            bad["DFAST"] += 1
        if not okb:
            bad["Bakta"] += 1
        rows.append({"accession": acc, "dfast_fixed": got_d,
                     "dfast_oracle": want_d, "dfast_ok": okd,
                     "bakta_fixed": got_b, "bakta_oracle": want_b,
                     "bakta_strict": strict_b, "bakta_pseudo_true": extra_b,
                     "bakta_ok": okb})
        if not (okd and okb):
            print(f"  MISMATCH {acc}: DFAST {got_d} vs {want_d}; "
                  f"Bakta intervals {got_b} = strict {strict_b} "
                  f"(summary {want_b}) + pseudo=True {extra_b}")

    tot_d = sum(r["dfast_fixed"] for r in rows)
    tot_b = sum(r["bakta_fixed"] for r in rows)
    print(f"\n{len(rows)} genomes checked")
    print(f"  DFAST pseudogene intervals, corrected total: {tot_d}"
          f"   mismatches: {bad['DFAST']}")
    print(f"  Bakta pseudogene intervals, corrected total: {tot_b}"
          f"   mismatches: {bad['Bakta']}")
    print(f"    of which pseudogene= rows matching Bakta's own summary: "
          f"{tot_b - extra_total}")
    print(f"    plus pseudo=True rows, which its summary omits:        "
          f"{extra_total}")
    print(f"  for comparison, the defective rule produced 125307 DFAST and "
          f"12068 Bakta")

    res = P / "benchmark/results/classifier_fix.json"
    res.write_text(json.dumps(
        {"n_genomes": len(rows), "dfast_total": tot_d, "bakta_total": tot_b,
         "bakta_pseudo_true": extra_total,
         "bakta_strict": tot_b - extra_total,
         "mismatches": dict(bad)}, indent=2, sort_keys=True) + "\n")
    with (P / "benchmark/results/classifier_fix.tsv").open("w") as fh:
        fh.write("accession\tdfast_fixed\tdfast_oracle\tbakta_fixed\t"
                 "bakta_oracle\tbakta_strict\tbakta_pseudo_true\n")
        for r in rows:
            fh.write(f"{r['accession']}\t{r['dfast_fixed']}\t{r['dfast_oracle']}"
                     f"\t{r['bakta_fixed']}\t{r['bakta_oracle']}"
                     f"\t{r['bakta_strict']}\t{r['bakta_pseudo_true']}\n")
    if sum(bad.values()):
        sys.exit(f"{sum(bad.values())} genomes disagree with the caller's own "
                 f"summary; the fix is not right yet")
    print("\nboth corrected paths reproduce the callers' own summary counts "
          "exactly")


if __name__ == "__main__":
    main()
