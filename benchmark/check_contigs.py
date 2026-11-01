#!/usr/bin/env python3
"""Did the caller keep the input's contig names? If not, the run has failed.

Plan section 5. This benchmark matches a lesion's known coordinate inside each
caller's own output, which only works if the caller's sequence identifiers are
the ones it was given. Bakta renames to contig_1.. unless --keep-contig-headers
is passed; DFAST renames to sequence1.. unless --use_original_name t is. The
parent study lived with the consequences: 26% of 886,176 interval lookups had
to be resolved by position in the file rather than by a unique length.

So this is a gate, not a diagnostic. A mismatch is reported as a failed run and
the assembly is dropped, rather than reconciled by length or order.

    python3 check_contigs.py --calls DIR --assemblies DIR
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

from Bio import SeqIO

HERE = pathlib.Path(__file__).resolve().parent


def input_contigs(fna: pathlib.Path) -> list[str]:
    return [r.id for r in SeqIO.parse(str(fna), "fasta")]


def gff_seqids(gff: pathlib.Path) -> list[str]:
    """Sequence ids in first-appearance order, from column 1 and ##sequence-region."""
    seen, out = set(), []
    with gff.open(errors="replace") as fh:
        for line in fh:
            if line.startswith("##sequence-region"):
                p = line.split()
                if len(p) > 1 and p[1] not in seen:
                    seen.add(p[1])
                    out.append(p[1])
            elif line.startswith("#"):
                continue
            else:
                p = line.split("\t")
                if len(p) > 8 and p[0] not in seen:
                    seen.add(p[0])
                    out.append(p[0])
    return out


OUTPUTS = {
    "bakta": lambda d, i: d / f"{i}.gff3",
    "dfast": lambda d, i: d / "genome.gff",
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True)
    ap.add_argument("--assemblies", required=True)
    a = ap.parse_args()
    calls, asm = pathlib.Path(a.calls), pathlib.Path(a.assemblies)

    rows, bad = [], 0
    for caller_dir in sorted(p for p in calls.iterdir()
                             if p.is_dir() and not p.name.startswith(".")):
        caller = caller_dir.name
        if caller not in OUTPUTS:
            print(f"  ?  {caller}: no output rule, skipped")
            continue
        for run in sorted(p for p in caller_dir.iterdir() if p.is_dir()):
            gff = OUTPUTS[caller](run, run.name)
            fna = asm / f"{run.name}.fna"
            if not gff.exists():
                print(f"  -  {caller}/{run.name}: no output yet")
                continue
            if not fna.exists():
                print(f"  !! {caller}/{run.name}: input assembly missing")
                bad += 1
                continue
            want, got = input_contigs(fna), gff_seqids(gff)
            missing = [c for c in want if c not in got]
            extra = [c for c in got if c not in want]
            ok = not missing and not extra
            rows.append({"caller": caller, "assembly": run.name, "ok": ok,
                         "n_input": len(want), "n_output": len(got),
                         "missing_from_output": missing[:6],
                         "unexpected_in_output": extra[:6]})
            if ok:
                print(f"  ok {caller}/{run.name}: all {len(want)} contig names "
                      f"preserved verbatim")
            else:
                bad += 1
                print(f"  !! {caller}/{run.name}: FAILED")
                print(f"       input  {want[:4]}")
                print(f"       output {got[:4]}")
                if extra:
                    print(f"       the caller renamed contigs; "
                          f"{len(extra)} names are not in the input")

    out = HERE / "results/contig_check.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rows, indent=2, sort_keys=True) + "\n")
    print(f"\n{len(rows)} runs checked, {bad} failed")
    if bad:
        sys.exit("contig names were not preserved; those runs are failures and "
                 "must not be reconciled by length or order")


if __name__ == "__main__":
    main()
