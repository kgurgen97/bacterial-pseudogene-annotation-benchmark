#!/usr/bin/env python3
"""Per-genome and per-genus tables for RBH_ORTHOLOGY_VALIDATION.md.

Reporting only. Reads results/rbh3/per_locus_rbh_audited.tsv, computes nothing
new about orthology and recomputes no omega; it reshapes what
analysis_rbh_audited.py already decided into the tables the write-up needs.

    python3 analysis_rbh_tables.py
"""
from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd

P = pathlib.Path(__file__).resolve().parent
OUT = P / "results" / "rbh3"
MIN_LOCI = 3

ARMS = [("omega", "E1 all"), ("omega_hiconf", "E2 RBH"),
        ("omega_standard", "E2 std"), ("omega_syn1", "E2 syn1"),
        ("omega_syn_collinear", "E2 collin"), ("omega_excluded", "E3 excl")]


def gap(sub, col):
    a = sub.loc[sub.kind == "discordant", col].dropna()
    b = sub.loc[sub.kind == "concordant_intact", col].dropna()
    if len(a) < MIN_LOCI or len(b) < MIN_LOCI:
        return None
    return float(np.median(a) - np.median(b))


def table(df, by):
    rows = []
    for key, sub in df.groupby(by):
        d = sub[sub.kind == "discordant"]
        c = sub[sub.kind == "concordant_intact"]
        r = {by: key, "n_dis": len(d), "n_con": len(c),
             "ret_dis": float(d.hiconf.mean()) if len(d) else None,
             "ret_con": float(c.hiconf.mean()) if len(c) else None,
             "fwdamb_dis": float(pd.to_numeric(
                 d.forward_ambiguous, errors="coerce").mean()) if len(d) else None,
             "cds_dis": float(d.has_annotated_cds.mean()) if len(d) else None}
        r["ret_imbalance_pp"] = (
            100 * (r["ret_con"] - r["ret_dis"])
            if r["ret_dis"] is not None and r["ret_con"] is not None else None)
        for col, lab in ARMS:
            r[lab] = gap(sub, col)
        rows.append(r)
    return pd.DataFrame(rows)


def md(t, by, label):
    print(f"\n### {label}\n")
    head = (f"| {by} | n dis | n con | retained dis | retained con | "
            f"imbalance pp | fwd-amb dis | has CDS dis | "
            + " | ".join(lab for _, lab in ARMS) + " |")
    print(head)
    print("|" + "---|" * (8 + len(ARMS)))
    for _, r in t.iterrows():
        cells = [str(r[by]), f"{r.n_dis:.0f}", f"{r.n_con:.0f}",
                 f"{100 * r.ret_dis:.1f} %" if pd.notna(r.ret_dis) else "—",
                 f"{100 * r.ret_con:.1f} %" if pd.notna(r.ret_con) else "—",
                 f"{r.ret_imbalance_pp:+.1f}" if pd.notna(r.ret_imbalance_pp) else "—",
                 f"{100 * r.fwdamb_dis:.1f} %" if pd.notna(r.fwdamb_dis) else "—",
                 f"{100 * r.cds_dis:.1f} %" if pd.notna(r.cds_dis) else "—"]
        for _, lab in ARMS:
            cells.append(f"{r[lab]:+.4f}" if pd.notna(r[lab]) else "—")
        print("| " + " | ".join(cells) + " |")


def main():
    df = pd.read_csv(OUT / "per_locus_rbh_audited.tsv", sep="\t")
    print(f"{len(df)} loci, {df.accession.nunique()} genomes, "
          f"{df.genus.nunique()} genera")
    for by, label in (("accession", "Per genome"), ("genus", "Per genus")):
        t = table(df, by).sort_values(by)
        t.to_csv(OUT / f"by_{by}.tsv", sep="\t", index=False)
        md(t, by, label)
    # how often does each arm keep the sign, across units?
    print("\n### Sign of the gap across genomes, per arm\n")
    t = table(df, "accession")
    print("| arm | genomes with a gap | positive | negative |")
    print("|---|---|---|---|")
    for _, lab in ARMS:
        v = t[lab].dropna()
        print(f"| {lab} | {len(v)} | {int((v > 0).sum())} | "
              f"{int((v < 0).sum())} |")


if __name__ == "__main__":
    main()
