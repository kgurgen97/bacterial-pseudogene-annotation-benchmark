#!/usr/bin/env python3
"""What did the classifier fix change? A leaf-by-leaf diff of every result json.

results/pre_fix/ holds the numbers the write-ups published; results/ holds the
recomputation. This compares them without knowing anything about their shape,
so a changed verdict cannot hide in a nested structure nobody thought to check.

Verdicts and p-values are reported first because those are what the prose
quotes.

    python3 benchmark/compare_pre_post_fix.py
"""
import json
import pathlib

P = pathlib.Path(__file__).resolve().parent.parent
PRE = P / "results/pre_fix"
POST = P / "results"

FILES = ["dnds.json", "locus.json", "clustering.json", "dfast.json",
         "frame.json"]


def leaves(obj, path=""):
    """Every scalar in a nested structure, keyed by its path."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from leaves(v, f"{path}.{k}" if path else str(k))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from leaves(v, f"{path}[{i}]")
    else:
        yield path, obj


def fmt(v):
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def changed(a, b):
    if isinstance(a, float) and isinstance(b, float):
        if a == b:
            return False
        if a == 0 or b == 0:
            return abs(a - b) > 1e-12
        return abs(a - b) / max(abs(a), abs(b)) > 1e-9
    return a != b


def main() -> None:
    summary = []
    for name in FILES:
        pre, post = PRE / name, POST / name
        print(f"\n{'=' * 70}\n{name}")
        if not pre.exists():
            print("  no pre-fix copy saved")
            continue
        if not post.exists():
            print("  not recomputed yet")
            continue
        a = dict(leaves(json.loads(pre.read_text())))
        b = dict(leaves(json.loads(post.read_text())))

        verdict_rows, number_rows, gone = [], [], []
        for k in sorted(set(a) | set(b)):
            if k not in a:
                gone.append((k, "(new)", fmt(b[k])))
                continue
            if k not in b:
                gone.append((k, fmt(a[k]), "(removed)"))
                continue
            if not changed(a[k], b[k]):
                continue
            row = (k, fmt(a[k]), fmt(b[k]))
            if "verdict" in k.lower():
                verdict_rows.append(row)
            else:
                number_rows.append(row)

        if verdict_rows:
            print("  VERDICTS THAT CHANGED")
            for k, x, y in verdict_rows:
                print(f"    {k}\n       was {x}\n       now {y}")
                summary.append((name, k, x, y))
        else:
            print("  no verdict changed")

        shown = [r for r in number_rows
                 if any(t in r[0].lower() for t in
                        ("p_value", "pvalue", ".p", "median", "n_", "frac",
                         "effect", "total", "count"))]
        rest = [r for r in number_rows if r not in shown]
        if shown:
            print(f"  headline numbers that changed ({len(shown)}):")
            for k, x, y in shown[:40]:
                print(f"    {k:54s} {x:>14s} -> {y:>14s}")
            if len(shown) > 40:
                print(f"    ... and {len(shown) - 40} more")
        if rest:
            print(f"  other values changed: {len(rest)}")
        if gone:
            print(f"  keys added or removed: {len(gone)}")
            for k, x, y in gone[:10]:
                print(f"    {k:54s} {x:>14s} -> {y:>14s}")

    print(f"\n{'=' * 70}\nVERDICT CHANGES, ALL FILES")
    if not summary:
        print("  none -- every verdict survived the fix")
    for name, k, x, y in summary:
        print(f"  {name:16s} {k:34s} {x} -> {y}")

    out = P / "benchmark/results/pre_post_fix.json"
    out.write_text(json.dumps(
        [{"file": n, "key": k, "before": x, "after": y}
         for n, k, x, y in summary], indent=2) + "\n")
    print(f"\nwrote {out.relative_to(P)}")


if __name__ == "__main__":
    main()
