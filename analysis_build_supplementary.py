#!/usr/bin/env python3
"""Assemble the supplementary tables from committed outputs.

Each table is written to supplementary/ with a header row, and an index
records what each one is and where it came from. Nothing is recomputed.

    python3 analysis_build_supplementary.py
"""
from __future__ import annotations

import hashlib
import json
import pathlib

import pandas as pd

P = pathlib.Path(__file__).resolve().parent
R = P / "results"
B = P / "benchmark" / "results"
SUP = P / "supplementary"


def write(df, name, title, source):
    SUP.mkdir(exist_ok=True)
    out = SUP / f"{name}.tsv"
    df.to_csv(out, sep="\t", index=False)
    sha = hashlib.sha256(out.read_bytes()).hexdigest()[:16]
    print(f"  {name:6s} {len(df):7d} rows  {title}")
    return {"table": name, "title": title, "rows": int(len(df)),
            "file": f"supplementary/{name}.tsv", "source": source,
            "sha256_16": sha}


def main() -> None:
    idx = []
    print("building supplementary tables")

    # S1 per-genome counts and fractions, four configurations
    cor = pd.read_csv(R / "corpus523.tsv", sep="\t")
    keep = ["accession", "stratum", "size_bp", "gc_percent",
            "pgap_pseudo", "pgap_cds", "pf_pseudo", "pf_intact",
            "bakta_pseudo", "bakta_cds", "dfast_disrupted", "dfast_cds",
            "frac_PGAP", "frac_Pseudofinder", "frac_Bakta", "frac_DFAST"]
    idx.append(write(cor[[c for c in keep if c in cor.columns]], "S1",
                     "Per-genome pseudogene counts and fractions, four "
                     "configurations, 523 genomes", "results/corpus523.tsv"))

    # S2 per-genus RBH outcome
    idx.append(write(pd.read_csv(R / "rbh3_all" / "by_genus.tsv", sep="\t"),
                     "S2", "Per-genus reciprocal-best-hit outcome, 28 genera",
                     "results/rbh3_all/by_genus.tsv"))

    # S3 the literature re-tests
    lit = json.loads((R / "literature_retest.json").read_text())
    rows = []
    for key, t in lit["C1"]["per_test"].items():
        parts = key.split("|")
        rows.append({"test": key, "level": parts[0],
                     "caller": parts[1] if len(parts) > 1 else "",
                     "metric": parts[2] if len(parts) > 2 else "",
                     "rho": t.get("rho"), "p": t.get("p"),
                     "partial_rho": t.get("partial_rho"),
                     "partial_p": t.get("partial_p"), "n": t.get("n"),
                     "verdict": t.get("verdict")})
    idx.append(write(pd.DataFrame(rows), "S3",
                     "Genome-reduction re-test: every caller x metric x "
                     "stratification test", "results/literature_retest.json"))

    # S4 direct between-configuration contrasts
    parts = []
    for m in ("fraction", "count", "per_Mb"):
        f = R / "caller_effect_contrast" / f"contrasts_{m}.tsv"
        if f.exists():
            parts.append(pd.read_csv(f, sep="\t"))
    if parts:
        idx.append(write(pd.concat(parts), "S4",
                         "Direct paired contrasts of effect size between "
                         "configurations, with cluster-bootstrap intervals",
                         "results/caller_effect_contrast/"))

    # S5 seed sensitivity
    seed = json.loads((R / "seed_sensitivity.json").read_text())
    rows = []
    for caller, arms in seed["per_caller"].items():
        for arm, blk in arms.items():
            rf = blk.get("rng_free", {})
            ac = blk.get("across_seeds", {})
            rows.append({
                "caller": caller, "arm": arm,
                "median_gap_original_rules":
                    rf.get("original_rules", {}).get("median_gap"),
                "sign_p_one_sided":
                    rf.get("original_rules", {}).get("sign_p_one_sided"),
                "genus_cluster_ci95": str(rf.get("genus_cluster_ci95")),
                "n_seeds": ac.get("n_seeds"),
                "seeds_p_below_0.05":
                    ac.get("original_rules", {}).get("n_seeds_p_below_0.05")})
    idx.append(write(pd.DataFrame(rows), "S5",
                     "Seed sensitivity across 40 fixed seeds, four "
                     "configurations", "results/seed_sensitivity.json"))

    # S6 FRAME2 arms
    f2 = json.loads((R / "frame2.json").read_text())
    rows = [{"arm": k, **v} for k, v in f2["arms"].items()]
    idx.append(write(pd.DataFrame(rows), "S6",
                     "Exact same-gene control: change in omega by lesion type",
                     "results/frame2.json"))

    # S7 lesion classification and controls
    idx.append(write(pd.read_csv(R / "lesion_extended" /
                                 "per_locus_extended.tsv", sep="\t"),
                     "S7", "Caller-independent lesion classification on "
                           "extended intervals, per locus",
                     "results/lesion_extended/per_locus_extended.tsv"))
    idx.append(write(pd.read_csv(R / "lesion_extended" / "controls.tsv",
                                 sep="\t"), "S7b",
                     "Lesion-classifier controls: artificial lesion and "
                     "unmodified arms",
                     "results/lesion_extended/controls.tsv"))

    # S8 the gate register, every gate including the failures
    gates = []
    for src, label in ((R / "rbh_audited_all.json", "RBH, 92 targets"),
                       (R / "lesion_extended.json", "lesion, extended"),
                       (R / "forward_ambiguity_all.json", "forward ambiguity"),
                       (R / "rbh_stratified_all.json", "stratified retention")):
        if not src.exists():
            continue
        d = json.loads(src.read_text())
        for key in ("gates", "gates_technical", "gates_scientific", "gate"):
            blk = d.get(key)
            if isinstance(blk, dict):
                for g, v in blk.items():
                    if isinstance(v, dict) and "pass" in v:
                        gates.append({"analysis": label, "block": key,
                                      "gate": g, "value": v.get("value"),
                                      "threshold": v.get("threshold"),
                                      "pass": v.get("pass")})
    idx.append(write(pd.DataFrame(gates), "S8",
                     "Gate register: every declared gate and whether it "
                     "passed, failures included", "results/*.json"))

    # S9 provenance: result registry files with digests
    idx.append(write(pd.read_csv(R / "result_registry_files.tsv", sep="\t"),
                     "S9", "Provenance: every result file with size and "
                           "SHA-256 prefix",
                     "results/result_registry_files.tsv"))

    man = pd.DataFrame(idx)
    man.to_csv(SUP / "INDEX.tsv", sep="\t", index=False)
    lines = ["# Supplementary materials", "",
             "Assembled by `analysis_build_supplementary.py` from committed "
             "outputs. Nothing here is recomputed.", "",
             "| table | title | rows | file |", "|---|---|---|---|"]
    for r in idx:
        lines.append(f"| **{r['table']}** | {r['title']} | {r['rows']:,} | "
                     f"`{r['file']}` |")
    lines += ["", "Full per-file provenance, including SHA-256 prefixes, is in "
                  "`supplementary/INDEX.tsv`."]
    (SUP / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\n{len(idx)} tables written to supplementary/ with INDEX.tsv "
          f"and README.md")


if __name__ == "__main__":
    main()
