#!/usr/bin/env python3
"""Independent audit of the corrected parser, beyond "the totals match".

Reproducing each caller's own summary counts proves the SIZE of the set is
right. It does not prove the intervals are the right intervals. Four things
could still be wrong, and each gets its own check:

  A1  coordinates. DFAST's pseudogene_summary.tsv states a location as
      `sequence1:6528..6896(-)`. If that is parsed wrongly -- off by one, or
      the strand dropped -- every downstream dN/dS is computed on the wrong
      sequence while the counts still match. Checked by looking up the same
      locus id in the GFF and demanding the coordinates agree exactly.

  A2  identifiers. The summary's locus_id must name a feature the GFF actually
      contains, or the intact set is being filtered against ids that do not
      exist and nothing is excluded.

  A3  duplication. If a caller emits both a `gene` and a `CDS` row for one
      locus, the interval extractor counts it twice. The counts could still
      match a summary that counts genes, so this has to be checked directly.

  A4  attribute values. The fixed test keys on `pseudo=` / `pseudogene=`
      existing. If any caller ever writes `pseudo=False`, a key test calls it
      a pseudogene. Every distinct key and value in the corpus is enumerated.

Also checked: no duplicate intervals inside one caller's pseudogene set, and
no interval shared between a caller's pseudogene and intact sets.

    python3 benchmark/audit_fixed_parser.py [n_genomes]
"""
import ast
import collections
import json
import pathlib
import re
import sys
import types

P = pathlib.Path(__file__).resolve().parent.parent
D = pathlib.Path.home() / "research/datasets"


def borrow(path):
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
FAIL = collections.Counter()
NOTE = collections.defaultdict(list)
# genes whose annotated bound runs off a contig edge, reported rather than
# silently tolerated: the interval is truncated at the edge by both paths
FUZZY: list[str] = []


def fail(kind, msg):
    FAIL[kind] += 1
    if len(NOTE[kind]) < 5:
        NOTE[kind].append(msg)


def gff_rows(path):
    if not path.exists():
        return
    with path.open(errors="replace") as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                return
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 9:
                yield f


def attr(f, key):
    for kv in f[8].split(";"):
        if kv.startswith(key + "="):
            return kv.split("=", 1)[1]
    return None


def audit_dfast(acc):
    gff = D / f"n2-dfast/{acc}/genome.gff"
    summ = D / f"n2-dfast/{acc}/pseudogene_summary.tsv"
    if not (gff.exists() and summ.exists()):
        return
    by_id, types_seen, tag_types = {}, collections.Counter(), collections.defaultdict(set)
    for f in gff_rows(gff):
        types_seen[f[2]] += 1
        fid = attr(f, "ID")
        if fid:
            by_id[fid] = f
        tag = attr(f, "locus_tag")
        if tag:
            tag_types[tag].add(f[2])

    # A3: a locus appearing as both gene and CDS would be counted twice
    for tag, ts in tag_types.items():
        if len({t for t in ts if t in ("gene", "CDS")}) > 1:
            fail("A3 duplicate gene/CDS", f"DFAST {acc} {tag}: {sorted(ts)}")

    seen_rows = 0
    for line in summ.read_text(errors="replace").splitlines()[1:]:
        f = line.split("\t")
        if len(f) < 4 or f[3] not in X.DFAST_DISRUPTED:
            continue
        seen_rows += 1
        lid, loc = f[0], f[2]
        m = re.match(r"^([^:]+):<?(\d+)\.\.>?(\d+)\(([+-])\)", loc)
        if not m:
            fail("A1 unparsable location", f"DFAST {acc} {lid}: {loc!r}")
            continue
        # A2: the id must exist in the GFF
        row = by_id.get(lid)
        if row is None:
            fail("A2 id missing from gff", f"DFAST {acc} {lid}")
            continue
        # A1: and its coordinates must agree. A gene running off a contig edge
        # gets a fuzzy bound -- the GFF writes `<1` where the summary writes
        # `1`. Both parser paths strip `<` and `>`, so the comparison must too;
        # comparing the raw strings flagged the only 2 such genes in 45,976 as
        # mismatches when nothing was wrong.
        g_start, g_end = row[3].lstrip("<>"), row[4].lstrip("<>")
        if (row[0] != m.group(1) or g_start != m.group(2)
                or g_end != m.group(3) or row[6] != m.group(4)):
            fail("A1 coordinate mismatch",
                 f"DFAST {acc} {lid}: summary {m.group(1)}:{m.group(2)}-"
                 f"{m.group(3)}({m.group(4)}) vs gff {row[0]}:{row[3]}-"
                 f"{row[4]}({row[6]})")
        if row[3].startswith(("<", ">")) or row[4].startswith(("<", ">")):
            FUZZY.append(f"{acc} {lid} {row[0]}:{row[3]}-{row[4]}")

    # the parser must return exactly the rows we just counted
    ids, ivs = X.dfast_pseudo(acc)
    if len(ivs) != seen_rows:
        fail("A1 count drift", f"DFAST {acc}: parser {len(ivs)} vs audit {seen_rows}")
    # no duplicate intervals, and pseudo and intact must not overlap
    if len(set(ivs)) != len(ivs):
        fail("duplicate intervals", f"DFAST {acc}: {len(ivs) - len(set(ivs))} repeats")
    intact = X.dfast_intervals(acc, False)
    if set(ivs) & set(intact):
        fail("pseudo in intact set",
             f"DFAST {acc}: {len(set(ivs) & set(intact))} intervals on both sides")
    return types_seen


def audit_bakta(acc):
    gff = D / f"n2-bakta/{acc}/{acc}.gff3"
    if not gff.exists():
        return
    keys = collections.Counter()
    tag_types = collections.defaultdict(set)
    types_seen = collections.Counter()
    for f in gff_rows(gff):
        types_seen[f[2]] += 1
        tag = attr(f, "locus_tag")
        if tag:
            tag_types[tag].add(f[2])
        for kv in f[8].split(";"):
            k = kv.split("=", 1)[0].strip().lower()
            if k.startswith("pseudo"):
                keys[(f[2], kv.strip())] += 1
    for tag, ts in tag_types.items():
        if len({t for t in ts if t in ("gene", "CDS")}) > 1:
            fail("A3 duplicate gene/CDS", f"Bakta {acc} {tag}: {sorted(ts)}")
    ps = X.calls_for(acc, "Bakta", True)
    it = X.calls_for(acc, "Bakta", False)
    if len(set(ps)) != len(ps):
        fail("duplicate intervals", f"Bakta {acc}: {len(ps) - len(set(ps))} repeats")
    if set(ps) & set(it):
        fail("pseudo in intact set",
             f"Bakta {acc}: {len(set(ps) & set(it))} intervals on both sides")
    return keys, types_seen


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 523
    accs = sorted(p.name for p in (D / "n2-dfast").iterdir() if p.is_dir())[:n]
    dfast_types, bakta_types = collections.Counter(), collections.Counter()
    attr_values = collections.Counter()

    for i, acc in enumerate(accs, 1):
        t = audit_dfast(acc)
        if t:
            dfast_types.update(t)
        r = audit_bakta(acc)
        if r:
            k, t2 = r
            attr_values.update(k)
            bakta_types.update(t2)
        if i % 100 == 0:
            print(f"  {i}/{len(accs)} genomes", flush=True)

    print(f"\n{len(accs)} genomes audited\n")
    print("A4  every pseudo* attribute in the corpus, by feature type and value")
    for (ftype, kv), c in sorted(attr_values.items(), key=lambda x: -x[1]):
        used = "COUNTED" if ftype in ("CDS", "gene", "pseudogene") else "ignored"
        print(f"    {ftype:18s} {kv:28s} {c:6d}   {used}")
    bad_vals = [kv for (_, kv), _ in attr_values.items()
                if re.search(r"=\s*(false|no|0)\s*$", kv, re.I)]
    print(f"    negated values that a key test would wrongly count: "
          f"{bad_vals if bad_vals else 'none'}")

    print("\nA3  feature types present (a 'gene' type alongside CDS risks "
          "double counting)")
    print(f"    DFAST: {dict(dfast_types.most_common(8))}")
    print(f"    Bakta: {dict(bakta_types.most_common(8))}")

    print(f"\nA1  pseudogenes with a fuzzy bound (gene runs off a contig edge): "
          f"{len(FUZZY)}")
    for x in FUZZY[:5]:
        print(f"    {x}")
    if FUZZY:
        print("    both the summary parser and the GFF reader strip < and >, so "
              "these are\n    truncated at the contig edge consistently")

    print("\nresults")
    if not FAIL:
        print("    all checks pass: coordinates agree, every summary id exists "
              "in the GFF,\n    no gene/CDS duplication, no duplicate or "
              "double-sided intervals")
    for kind, c in FAIL.most_common():
        print(f"    FAIL {kind}: {c}")
        for m in NOTE[kind]:
            print(f"         {m}")

    out = P / "benchmark/results/parser_audit.json"
    out.write_text(json.dumps(
        {"n_genomes": len(accs), "failures": dict(FAIL),
         "examples": {k: v for k, v in NOTE.items()},
         "fuzzy_bound_pseudogenes": FUZZY,
         "bakta_pseudo_attributes": {f"{t}|{kv}": c
                                     for (t, kv), c in attr_values.items()},
         "dfast_feature_types": dict(dfast_types),
         "bakta_feature_types": dict(bakta_types),
         "negated_values": bad_vals}, indent=2, sort_keys=True) + "\n")
    print(f"\nwrote {out.relative_to(P)}")
    if FAIL:
        sys.exit(f"{sum(FAIL.values())} audit failures")


if __name__ == "__main__":
    main()
