#!/usr/bin/env python3
"""Pre- and post-flight checks for the 92-target RBH expansion.

Asserts RBH_EXPANSION_PLAN.md's hash, and the hashes of the two plans it adopts
by reference, so an expanded run cannot quietly proceed under an edited plan.
Then records a manifest of the pilot evidence, which the same script re-verifies
afterwards: the expansion must leave every pilot file byte-identical, and that
is checked rather than asserted in prose.

    python3 analysis_rbh_expansion_check.py before
    python3 analysis_rbh_expansion_check.py after
"""
from __future__ import annotations

import hashlib
import json
import pathlib
import subprocess
import sys

import pandas as pd

P = pathlib.Path(__file__).resolve().parent
R = P / "results"
MANIFEST = R / "rbh_expansion_manifest.json"

PLANS = {
    "RBH_EXPANSION_PLAN.md":
        "c32bcc1a434168779b0b0cb00e4373be2bc7139b72274489ad0342c83d5a129a",
    "RBH_STRATIFIED_PLAN.md":
        "1fadc6db0ee4310b8004a9c488deb262d9d0815687423200edafa3e7fc77d4ef",
    "FORWARD_AMBIGUITY_PLAN.md":
        "e626f6d4386c4e824bd0872ead39ddee0206ef30c60ad6753a644365edb49698",
    "FRAME2_PLAN.md":
        None,                      # recorded, not pinned here
}

# Everything the expansion must not touch. DATA and EVIDENCE only.
#
# RBH_ORTHOLOGY_VALIDATION.md was in this list on the first pass and the
# post-flight check duly flagged it, because writing the expansion up is
# exactly what must change it. Protecting the write-up against the run that it
# documents is a category error: the guard exists so the run cannot disturb the
# numbers it is being compared against, and prose under version control has its
# own history. It is listed under RELEASED below so the removal is visible
# rather than silent.
PROTECTED = [
    "results/loci.tsv", "results/locus.json", "MANUSCRIPT.md",
    "results/reference_sensitivity/per_locus.tsv",
    "results/reference_sensitivity/per_locus_all.tsv",
    "results/rbh3/per_locus_rbh_audited.tsv",
    "results/rbh3/input_snapshot.tsv",
    "results/rbh3/matching_grid.tsv",
    "results/rbh_audited.json",
    "results/rbh_stratified.json",
    "results/rbh_stratified/per_genome.tsv",
    "results/forward_ambiguity.json",
    "results/forward_ambiguity/runner_up_homology.tsv",
]

# Recorded by an earlier pre-flight, then deliberately released. A manifest
# entry that is no longer protected is reported, never silently ignored.
RELEASED = {
    "RBH_ORTHOLOGY_VALIDATION.md":
        "the write-up of this very run; protecting it against itself was a "
        "category error, and git carries its history",
}

INPUT = "results/reference_sensitivity/per_locus_all.tsv"
EXPECT_LOCI = 36535
EXPECT_PAIRS = 92


def fsha(p: pathlib.Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def host_free_gb():
    out = subprocess.run(["df", "--output=avail", "-BG", "/mnt/c"],
                         capture_output=True, text=True, check=True).stdout
    return int("".join(c for c in out.splitlines()[-1] if c.isdigit()))


def before():
    print("PRE-FLIGHT")
    free = host_free_gb()
    print(f"  host C: free {free} GB (from /mnt/c)")
    if free < 10:
        sys.exit(f"  refusing: {free} GB < 10 GB start floor")

    for name, want in PLANS.items():
        p = P / name
        if not p.exists():
            sys.exit(f"  missing plan {name}")
        got = fsha(p)
        if want is None:
            print(f"  {name:28s} {got[:16]}…  (recorded)")
            continue
        if got != want:
            sys.exit(f"  FROZEN PLAN CHANGED: {name}\n"
                     f"    expected {want}\n    found    {got}")
        print(f"  {name:28s} {got[:16]}…  OK")

    src = P / INPUT
    df = pd.read_csv(src, sep="\t")
    pairs = df[["accession", "partner"]].drop_duplicates()
    print(f"  input {INPUT}")
    print(f"    sha256 {fsha(src)[:16]}…")
    print(f"    {len(df)} rows, {len(pairs)} pairs, {df.genus.nunique()} genera,"
          f" {df.accession.nunique()} genomes")
    if len(df) != EXPECT_LOCI or len(pairs) != EXPECT_PAIRS:
        sys.exit(f"  expected {EXPECT_LOCI} rows and {EXPECT_PAIRS} pairs")
    reps = df.groupby("genus").accession.nunique()
    multi = int((reps > 1).sum())
    print(f"    genera with more than one genome: {multi} of {len(reps)}"
          f"  (max {int(reps.max())} genomes in a genus)")
    if multi < 1:
        sys.exit("  no genus has repeated genomes; expansion would add nothing")

    man = {p_: fsha(P / p_) for p_ in PROTECTED if (P / p_).exists()}
    MANIFEST.write_text(json.dumps({"input_sha256": fsha(src),
                                    "protected": man}, indent=1,
                                   sort_keys=True))
    print(f"  recorded {len(man)} protected files -> "
          f"{MANIFEST.relative_to(P)}")
    print("  PRE-FLIGHT OK")


def after():
    print("POST-FLIGHT")
    if not MANIFEST.exists():
        sys.exit("  no manifest; run 'before' first")
    man = json.loads(MANIFEST.read_text())
    bad, skipped = [], []
    for p_, want in man["protected"].items():
        if p_ not in PROTECTED:
            skipped.append(p_)
            continue
        f = P / p_
        if not f.exists():
            bad.append((p_, "MISSING"))
        elif fsha(f) != want:
            bad.append((p_, "MODIFIED"))
    for p_ in skipped:
        why = RELEASED.get(p_, "no reason recorded — investigate")
        print(f"  RELEASED, not checked: {p_}\n    {why}")
    src = P / INPUT
    if fsha(src) != man["input_sha256"]:
        bad.append((INPUT, "INPUT CHANGED DURING THE RUN"))
    if bad:
        for p_, why in bad:
            print(f"  {why}: {p_}")
        sys.exit(f"  {len(bad)} protected file(s) disturbed")
    n = len(man["protected"]) - len(skipped)
    print(f"  all {n} protected data files byte-identical")
    print(f"  input unchanged across the run")
    print(f"  host C: free {host_free_gb()} GB")
    print("  POST-FLIGHT OK")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "before"
    {"before": before, "after": after}.get(mode, before)()
