#!/usr/bin/env python3
"""How far does Pseudofinder still reach into G3 after the structural fix?

POST HOC. Not registered. PF_INPUT_DEPENDENCE.md established that Pseudofinder
post-processes PGAP's annotation, so the two cannot disagree where PGAP called a
pseudogene. Two residual dependencies are NOT covered by that finding and are
tested here.

  D1  REFERENCE SELECTION. Every omega in this project is computed against the
      PARTNER genome's pf_intact.ffn / pf_intact.faa - which is Pseudofinder's
      INTACT set for that partner. Any partner gene Pseudofinder called a
      pseudogene is therefore absent from the orthologue reference, so it can
      never serve as the comparison for a query gene. That is a Pseudofinder
      filter sitting underneath every dN/dS value in the study, including
      PGAP's, Bakta's and DFAST's.

  D2  THE BASELINE. A concordant-intact locus needs two or more callers to
      agree the gene is intact, and Pseudofinder is one of the four. If
      Pseudofinder's agreement is load-bearing for the baseline, the control
      side of G3 is partly Pseudofinder's construct too.

Both are quantified, and D2 is tested directly by rebuilding G3 with
Pseudofinder excluded from the baseline entirely.

    python3 analysis_pf_residual_dependence.py
"""
import collections
import csv
import json
import pathlib

import numpy as np
import pandas as pd
from scipy import stats

P = pathlib.Path(__file__).resolve().parent
R = P / "results"
D = pathlib.Path.home() / "research/datasets"
MIN_LOCI = 3


def g3(t, min_loci=MIN_LOCI):
    d = t[t.kind == "discordant"]
    c = t[t.kind == "concordant_intact"]
    gaps = {}
    for acc in sorted(set(d.accession) & set(c.accession)):
        a = d.loc[d.accession == acc, "omega"]
        b = c.loc[c.accession == acc, "omega"]
        if len(a) >= min_loci and len(b) >= min_loci:
            gaps[acc] = float(np.median(a) - np.median(b))
    nz = [v for v in gaps.values() if v != 0.0]
    if len(nz) < 5:
        return {"n_genomes": len(gaps), "p": None, "median_gap": None,
                "note": f"only {len(nz)} non-tied genomes"}
    pos = sum(1 for v in nz if v > 0)
    return {"n_genomes": len(gaps), "n_tested": len(nz),
            "median_gap": float(np.median(list(gaps.values()))),
            "n_positive": pos,
            "p": float(stats.binomtest(pos, len(nz), 0.5).pvalue)}


def count_fasta(p):
    if not p.exists():
        return 0
    n = 0
    with p.open(errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def main() -> None:
    out = {"note": "POST HOC; residual Pseudofinder dependence in G3"}

    # ---------------------------------------------------------------- D1
    print("D1  the orthologue reference is Pseudofinder's intact set")
    partners = {}
    with (R / "dnds_pairs.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            partners.setdefault(r["accession"], r["partner"])
    uniq = sorted(set(partners.values()))
    rows = []
    for pt in uniq:
        intact = count_fasta(D / f"n2-pseudofinder/{pt}/pf_intact.ffn")
        pseud = 0
        g = D / f"n2-pseudofinder/{pt}/pf_pseudos.gff"
        if g.exists():
            pseud = sum(1 for line in g.read_text(errors="replace").splitlines()
                        if line and not line.startswith("#"))
        rows.append({"partner": pt, "reference_genes": intact,
                     "excluded_as_pseudogene": pseud,
                     "frac_excluded": pseud / (intact + pseud) if intact + pseud else None})
    df = pd.DataFrame(rows)
    tot_i = int(df.reference_genes.sum())
    tot_p = int(df.excluded_as_pseudogene.sum())
    print(f"  {len(uniq)} distinct partner genomes used as references")
    print(f"  genes available as orthologues : {tot_i}")
    print(f"  genes Pseudofinder withheld    : {tot_p} "
          f"({100 * tot_p / (tot_i + tot_p):.1f}% of the partner gene sets)")
    print(f"  per partner, withheld fraction: median "
          f"{df.frac_excluded.median():.3f}, max {df.frac_excluded.max():.3f}")
    print("  every dN/dS in the study - PGAP's, Bakta's and DFAST's included -")
    print("  is computed against a reference Pseudofinder filtered.")
    out["D1_reference"] = {
        "n_partners": len(uniq), "reference_genes": tot_i,
        "withheld_as_pseudogene": tot_p,
        "frac_withheld": tot_p / (tot_i + tot_p),
        "per_partner_median_frac": float(df.frac_excluded.median()),
        "per_partner_max_frac": float(df.frac_excluded.max())}
    df.to_csv(R / "pf_reference_exposure.tsv", sep="\t", index=False)

    # ---------------------------------------------------------------- D2
    print("\nD2  how load-bearing is Pseudofinder in the concordant-intact baseline?")
    t = pd.read_csv(R / "loci.tsv", sep="\t")
    base = t[t.kind == "concordant_intact"].copy()
    sets = base.intact_callers.fillna("").apply(
        lambda s: set(x for x in str(s).split("+") if x and x != "-"))
    base["n_callers"] = sets.apply(len)
    base["has_pf"] = sets.apply(lambda s: "Pseudofinder" in s)
    n_all4 = int((base.n_callers == 4).sum())
    n_pf = int(base.has_pf.sum())
    # loci that would FALL BELOW the two-caller bar without Pseudofinder
    would_fail = int(((base.n_callers == 2) & base.has_pf).sum())
    print(f"  baseline loci {len(base)}")
    print(f"    all four callers agree intact : {n_all4} "
          f"({100 * n_all4 / len(base):.1f}%)")
    print(f"    Pseudofinder among them       : {n_pf} "
          f"({100 * n_pf / len(base):.1f}%)")
    print(f"    rest on exactly two, one being Pseudofinder: {would_fail} "
          f"({100 * would_fail / len(base):.1f}%)")

    # rebuild G3 with Pseudofinder removed from the baseline's intact set
    keep_idx = base.index[~((base.n_callers == 2) & base.has_pf)]
    t_nopf = pd.concat([t[t.kind == "discordant"], base.loc[keep_idx]])
    r_full = g3(t)
    r_nopf = g3(t_nopf)
    # and the strictest version: baseline must hold WITHOUT Pseudofinder at all
    strict_idx = base.index[(base.n_callers - base.has_pf.astype(int)) >= 2]
    t_strict = pd.concat([t[t.kind == "discordant"], base.loc[strict_idx]])
    r_strict = g3(t_strict)
    print(f"\n  G3, baseline as published            "
          f"{r_full['median_gap']:+.4f}  {r_full['n_positive']}/{r_full['n_tested']}"
          f"  p {r_full['p']:.3g}")
    print(f"  G3, dropping PF-dependent 2-caller   "
          f"{r_nopf['median_gap']:+.4f}  {r_nopf['n_positive']}/{r_nopf['n_tested']}"
          f"  p {r_nopf['p']:.3g}   (-{len(base) - len(keep_idx)} baseline loci)")
    print(f"  G3, baseline must hold WITHOUT PF    "
          f"{r_strict['median_gap']:+.4f}  "
          f"{r_strict['n_positive']}/{r_strict['n_tested']}"
          f"  p {r_strict['p']:.3g}   (-{len(base) - len(strict_idx)} baseline loci)")
    out["D2_baseline"] = {
        "n_baseline": int(len(base)), "n_all_four_agree": n_all4,
        "n_with_pseudofinder": n_pf,
        "n_two_caller_dependent_on_pf": would_fail,
        "g3_as_published": r_full,
        "g3_without_pf_dependent_pairs": r_nopf,
        "g3_baseline_without_pf_at_all": r_strict}

    print("\nverdict")
    if (r_strict.get("p") is not None and r_strict["p"] < 0.05
            and r_strict["median_gap"] > 0):
        print("  D2 CLEARED. G3 holds when the baseline is required to stand")
        print("  without Pseudofinder at all. The control side of G3 is not")
        print("  Pseudofinder's construct.")
        out["D2_verdict"] = "CLEARED"
    else:
        print("  D2 NOT CLEARED: the baseline depends on Pseudofinder.")
        out["D2_verdict"] = "NOT_CLEARED"
    print("  D1 CANNOT BE CLEARED with these data. The reference is")
    print("  Pseudofinder-filtered by construction, and testing it would need")
    print("  omega recomputed against an independent orthologue set.")
    out["D1_verdict"] = "UNRESOLVED_BY_DESIGN"

    (R / "pf_residual_dependence.json").write_text(
        json.dumps(out, indent=2, sort_keys=True) + "\n")
    print("\nwrote results/pf_residual_dependence.json")


if __name__ == "__main__":
    main()
