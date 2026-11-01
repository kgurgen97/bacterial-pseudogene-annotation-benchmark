#!/usr/bin/env python3
"""Which numbers in MANUSCRIPT.md are now stale? Read-only audit.

MANUSCRIPT.md is NOT modified. This locates every claim whose value changed
under the classifier correction, the RNG finding, or the new frame analysis, so
a revision proposal can be written against a list rather than an impression.

For each stale figure: the old value, the corrected value, where it appears, and
why it moved.

    python3 benchmark/audit_manuscript_claims.py
"""
import json
import pathlib
import re

P = pathlib.Path(__file__).resolve().parent.parent
MAN = P / "MANUSCRIPT.md"

J = {n: json.loads((P / f"results/{n}.json").read_text())
     for n in ("locus", "dfast", "dnds", "clustering", "frame", "frame2",
               "caller_comparison", "locus_sensitivity")}

# (label, regexes that match the OLD value, corrected value, why)
CHECKS = [
    ("discordant loci", [r"7[ ,]?565"], "4 886",
     "classifier fix: DFAST's 125,307 name-match calls removed"),
    ("scored loci", [r"103[ ,]?386"], "103 224", "classifier fix"),
    ("concordant-intact loci", [r"95[ ,]?821"], "98 338", "classifier fix"),
    ("G3 per-genome median", [r"\+0\.022\b"], "+0.031", "classifier fix"),
    ("G3 p-value", [r"4\.03\s*[x×]\s*10|4\s*[x×]\s*10.{0,3}12|4\.03e-12"],
     "6.16e-12", "classifier fix"),
    ("G3 genome count", [r"over\s+96\s+genomes"], "92 genomes", "classifier fix"),
    ("Pseudofinder-only loci", [r"1[ ,]?955"], "1 742", "classifier fix"),
    ("DFAST-only loci", [r"2[ ,]?639"], "540",
     "classifier fix: DFAST's calls were product names"),
    ("Pseudofinder-only gap", [r"\+0\.035\b"], "+0.0379", "classifier fix"),
    ("DFAST-only gap", [r"\+0\.006\b"], "+0.0389",
     "classifier fix: sign and significance both change"),
    ("DFAST-only p", [r"p\s*=?\s*0\.06\b"], "3.31e-06", "classifier fix"),
    ("L3 multi-caller fraction", [r"0\.1\s*%\s*of loci"], "0.2 %", "classifier fix"),
    ("total intervals", [r"886[ ,]?176"], "886 139", "classifier fix"),
    ("G2 best-worst spread", [r"spread 0\.035|0\.035\s*\(Pseudofinder"], "0.015",
     "classifier fix: DFAST moves from +0.002 to +0.023"),
    ("DFAST dN/dS gap", [r"\+0\.002\b"], "+0.023", "classifier fix"),
    ("PGAP clustered sign test", [r"0\.015\b"], "0.0967 and seed-dependent",
     "shared-RNG coupling; the published value was one stream position"),
    ("Bakta clustered sign test", [r"0\.0023\b|0\.00228"], "0.113",
     "classifier fix: Bakta drops to 11 usable genomes"),
]

# claims that are now INVALID rather than merely stale
INVALID = [
    ("frame analysis conclusions",
     [r"three-quarters of the separation", r"F4\b.*CONFIRMED",
      r"the objection was the strongest one raised"],
     "FRAME_RESULTS.md is VOID: F4 fired (DFAST pre gap +0.018 against a "
     "+0.010 bar). FRAME2_RESULTS.md replaces the question and finds a "
     "frameshift inflates omega by +0.0315 on its own."),
    ("DFAST characterised as making uncorroborated calls",
     [r"DFAST still does not", r"DFAST.{0,40}does not",
      r"uncorroborated"],
     "Inverted. DFAST's corrected unique calls are corroborated: +0.0389, "
     "40/48 genomes, p = 3.3e-06."),
    ("four independent callers",
     [r"four callers", r"all four", r"four annotation"],
     "Pseudofinder post-processes PGAP's annotation, so this is three "
     "independent annotators plus one dependent post-processor."),
]


def main() -> None:
    text = MAN.read_text()
    lines = text.splitlines()
    print(f"MANUSCRIPT.md: {len(lines)} lines, {len(text)} chars")
    print("NOT MODIFIED by this script.\n")

    found, missing, report = 0, [], {"stale": [], "invalid": []}
    print("=" * 76)
    print("STALE NUMBERS: a value that changed under a corrected analysis")
    for label, pats, new, why in CHECKS:
        hits = []
        for i, l in enumerate(lines, 1):
            for p in pats:
                if re.search(p, l, re.I):
                    hits.append((i, l.strip()[:110]))
                    break
        if hits:
            found += 1
            print(f"\n  {label}  ->  corrected: {new}")
            print(f"    why: {why}")
            for i, l in hits[:4]:
                print(f"    line {i}: {l}")
            if len(hits) > 4:
                print(f"    ... and {len(hits) - 4} more occurrences")
            report["stale"].append({"label": label, "corrected": new,
                                    "why": why,
                                    "lines": [i for i, _ in hits]})
        else:
            missing.append(label)

    print("\n" + "=" * 76)
    print("INVALID CLAIMS: not a number but a conclusion that no longer holds")
    for label, pats, why in INVALID:
        hits = []
        for i, l in enumerate(lines, 1):
            for p in pats:
                if re.search(p, l, re.I):
                    hits.append((i, l.strip()[:110]))
                    break
        if hits:
            print(f"\n  {label}  ({len(hits)} occurrences)")
            print(f"    {why}")
            for i, l in hits[:5]:
                print(f"    line {i}: {l}")
            if len(hits) > 5:
                print(f"    ... and {len(hits) - 5} more")
            report["invalid"].append({"label": label, "why": why,
                                      "n_occurrences": len(hits),
                                      "lines": [i for i, _ in hits]})

    print("\n" + "=" * 76)
    print(f"stale figures located in the text : {found} of {len(CHECKS)}")
    if missing:
        print("  not found in the manuscript (may be phrased differently,")
        print("  or may genuinely not appear):")
        for m in missing:
            print(f"    {m}")

    # the corrected values, for the proposal to quote
    print("\ncorrected values, from the json:")
    l = J["locus"]
    print(f"  locus: {l['n_loci']} loci, {l['n_discordant']} discordant, "
          f"{l['n_concordant_intact']} concordant-intact")
    print(f"  G3   : {l['verdicts']['G3']['value']}")
    s = J["dfast"]["sole_caller_vs_baseline"]
    for c in ("Pseudofinder", "DFAST"):
        r = s[c]
        print(f"  {c:13s} sole {r['n_loci']} loci, {r['median_diff']:+.4f}, "
              f"{r['n_positive']}/{r['n_tested']}, p {r['p']:.3g}")
    f2 = J["frame2"]["arms"]["frameshift_high"]
    print(f"  FRAME2 frameshift inflation: median {f2['median_delta']:+.4f} "
          f"(n {f2['n']})")
    m = J["caller_comparison"]["M1_matched_pairs"]["Pseudofinder vs DFAST"]
    print(f"  matched PF vs DFAST: p {m['p']:.3g} over {m['n_genomes']} genomes")
    ls = J["locus_sensitivity"]["summary"]
    print(f"  G3 sensitivity: {ls['n_positive_significant']}/"
          f"{ls['n_testable']} checks positive and significant")

    (P / "benchmark/results/manuscript_audit.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("\nwrote benchmark/results/manuscript_audit.json")


if __name__ == "__main__":
    main()
