#!/usr/bin/env python3
"""Does the parent study's pseudogene rule match product names?

analysis_dnds.py decides a GFF feature is a pseudogene with:

    is_pseudo = ("pseudo" in f[8].lower()) or f[2] == "pseudogene"

That is a substring search over the whole attributes column, which holds
product names, gene names, Dbxrefs and free-text notes. Reading one Bakta and
one DFAST output showed the attribute strings contain "pseudouridine",
"pseudoknot", "Pseudomonadota" and "pseudotuberculosis" -- all of them
properties of perfectly intact genes.

This counts, on the parent study's own caller output, how many features the
rule calls pseudogenes and how many of those carry an actual pseudogene
assertion rather than merely the letters "pseudo" in a name. It changes no
result by itself; it measures the size of a possible error.

    python3 check_substring_classifier.py [n_genomes]
"""
import collections
import csv
import json
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
D = pathlib.Path.home() / "research/datasets"

# What a real pseudogene assertion looks like in each format.
REAL = re.compile(r"(?:^|;)\s*(?:pseudo|pseudogene)\s*=", re.I)


def rows(p):
    """GFF3 data lines, stopping at an appended ##FASTA block."""
    with p.open(errors="replace") as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                return
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 9:
                yield f


def audit(p):
    """Replicate the parent rule, then split its hits into real and spurious."""
    flagged = real = 0
    why = collections.Counter()
    examples = []
    for f in rows(p):
        if f[2] not in ("CDS", "gene", "pseudogene"):
            continue
        attrs = f[8]
        hit = ("pseudo" in attrs.lower()) or f[2] == "pseudogene"
        if not hit:
            continue
        flagged += 1
        if REAL.search(attrs) or f[2] == "pseudogene":
            real += 1
            continue
        # not a pseudogene assertion: find which word triggered it
        for m in re.finditer(r"\w*pseudo\w*", attrs, re.I):
            why[m.group(0).lower()] += 1
        if len(examples) < 3:
            examples.append(f"{f[0]} {f[2]} {f[3]}-{f[4]} :: {attrs[:150]}")
    return {"flagged": flagged, "real": real, "spurious": flagged - real,
            "trigger_words": dict(why.most_common(8)), "examples": examples}


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 25
    out, totals = [], collections.defaultdict(lambda: collections.Counter())
    words = collections.defaultdict(collections.Counter)
    shown = collections.defaultdict(list)

    accs = sorted(p.name for p in (D / "n2-dfast").iterdir() if p.is_dir())[:n]
    for acc in accs:
        for caller, path in (("DFAST", D / f"n2-dfast/{acc}/genome.gff"),
                             ("Bakta", next((D / "n2-bakta" / acc).glob("*.gff3"),
                                            D / f"n2-bakta/{acc}/none"))):
            if not path.exists():
                continue
            a = audit(path)
            a.update(accession=acc, caller=caller)
            out.append(a)
            totals[caller]["flagged"] += a["flagged"]
            totals[caller]["real"] += a["real"]
            totals[caller]["spurious"] += a["spurious"]
            words[caller].update(a["trigger_words"])
            if a["spurious"] and len(shown[caller]) < 3:
                shown[caller] += a["examples"][:1]

    print(f"{len(accs)} genomes from the parent study's own caller output\n")
    for caller in ("DFAST", "Bakta"):
        t = totals[caller]
        if not t["flagged"]:
            print(f"{caller}: no output found")
            continue
        pct = 100 * t["spurious"] / t["flagged"]
        print(f"{caller}")
        print(f"  features the parent rule calls pseudogenes : {t['flagged']}")
        print(f"  of those, carrying a real pseudo= assertion: {t['real']}")
        print(f"  of those, matched only by a name or note   : {t['spurious']}"
              f"  ({pct:.1f}%)")
        if words[caller]:
            print("  words that triggered the spurious hits:")
            for w, c in words[caller].most_common(8):
                print(f"    {w:28s} {c}")
        for e in shown[caller]:
            print(f"    e.g. {e}")
        print()

    res = ROOT / "benchmark/results/substring_classifier.json"
    res.parent.mkdir(parents=True, exist_ok=True)
    res.write_text(json.dumps(
        {"n_genomes": len(accs),
         "totals": {k: dict(v) for k, v in totals.items()},
         "trigger_words": {k: dict(v) for k, v in words.items()}},
        indent=2, sort_keys=True) + "\n")
    with (ROOT / "benchmark/results/substring_classifier.tsv").open(
            "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t",
                           fieldnames=["accession", "caller", "flagged", "real",
                                       "spurious"], extrasaction="ignore")
        w.writeheader()
        w.writerows(out)
    print(f"wrote benchmark/results/substring_classifier.{{json,tsv}}")


if __name__ == "__main__":
    main()
