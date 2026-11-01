#!/usr/bin/env python3
"""Build the machine-readable result registry.

Walks every committed result JSON and emits one row per scalar leaf, so that
any figure quoted in the manuscript can be located by key rather than by
memory. Reporting only: nothing is recomputed.

    python3 analysis_result_registry.py
"""
from __future__ import annotations

import hashlib
import json
import pathlib

import pandas as pd

P = pathlib.Path(__file__).resolve().parent
R = P / "results"
OUT = R / "result_registry.tsv"
MAXDEPTH = 6


def walk(obj, prefix, rows, src, depth=0):
    if depth > MAXDEPTH:
        return
    if isinstance(obj, dict):
        for k, v in obj.items():
            walk(v, f"{prefix}.{k}" if prefix else str(k), rows, src, depth + 1)
    elif isinstance(obj, list):
        if obj and all(isinstance(x, (int, float)) for x in obj):
            rows.append((src, prefix, "list", json.dumps(obj)))
        else:
            for i, v in enumerate(obj[:40]):
                walk(v, f"{prefix}[{i}]", rows, src, depth + 1)
    elif isinstance(obj, (int, float, bool, str)) or obj is None:
        t = ("bool" if isinstance(obj, bool) else
             "number" if isinstance(obj, (int, float)) else
             "null" if obj is None else "text")
        s = str(obj)
        rows.append((src, prefix, t, s if len(s) <= 300 else s[:297] + "..."))


def main() -> None:
    rows, files = [], []
    # both result trees. A first version walked only results/ and so could not
    # resolve any figure sourced from benchmark/results/, which is where the
    # simulation, parser-audit and spread outputs live.
    sources = sorted(R.glob("*.json")) + sorted((P / "benchmark" / "results")
                                                .glob("*.json"))
    for f in sources:
        try:
            d = json.loads(f.read_text())
        except Exception as e:                    # keep going; record the fault
            files.append((f.name, f.stat().st_size, "UNPARSEABLE", str(e)[:80]))
            continue
        before = len(rows)
        walk(d, "", rows, f.name)
        files.append((f.name, f.stat().st_size,
                      hashlib.sha256(f.read_bytes()).hexdigest()[:16],
                      f"{len(rows) - before} leaves"))
    reg = pd.DataFrame(rows, columns=["source_file", "key", "type", "value"])
    reg.to_csv(OUT, sep="\t", index=False)
    man = pd.DataFrame(files, columns=["file", "bytes", "sha256_16", "note"])
    man.to_csv(R / "result_registry_files.tsv", sep="\t", index=False)
    print(f"{len(reg)} entries from {len(files)} result files")
    print(f"wrote {OUT.relative_to(P)} and result_registry_files.tsv")
    bad = man[man.sha256_16 == "UNPARSEABLE"]
    if len(bad):
        print("UNPARSEABLE:")
        print(bad.to_string(index=False))
    # a few spot values the manuscript depends on, printed for eyeballing
    want = [
        ("locus.json", "n_loci"), ("locus.json", "n_discordant"),
        ("rbh_audited_all.json", "n_loci"),
        ("reference_sensitivity_all.json", "coverage.n_loci"),
        ("frame2.json",), ("review_response.json", "R1_genus_level.p"),
    ]
    print("\nspot checks")
    for w in want:
        sel = reg[reg.source_file == w[0]]
        if len(w) > 1:
            sel = sel[sel.key == w[1]]
            for _, r in sel.iterrows():
                print(f"  {r.source_file:34s} {r.key:28s} {r.value}")
        else:
            print(f"  {w[0]:34s} {len(sel)} entries")


if __name__ == "__main__":
    main()
