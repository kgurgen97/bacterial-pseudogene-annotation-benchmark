#!/usr/bin/env python3
"""The forward-ambiguity criterion, corrected, and the asymmetry diagnosed.

Executes FORWARD_AMBIGUITY_PLAN.md, whose SHA-256 is asserted below.

THE DEFECT BEING FIXED. analysis_rbh_audited.py scored "forward ambiguity
<= 10%" POOLED: 5.46%, PASS. Per class it is 1.9% for concordant-intact loci
and 13.9% for discordant ones - a sevenfold difference - so the pooled pass
came from baseline loci outnumbering discordant loci four to one. The criterion
measured the wrong thing in the wrong way.

WHAT REPLACES IT, AND WHY THAT IS NOT CHEATING. A rate threshold is a proxy
whose relationship to bias is unknown. Omega is known for every ambiguous
locus - ambiguity is a recorded property, not a reason a locus is missing - so
the CONSEQUENCE of ambiguity is directly measurable and a proxy is unnecessary.
Replacing a failing proxy with a measurable consequence is otherwise
indistinguishable from choosing the criterion that passes, so both verdicts are
printed side by side everywhere, and the proxy's class-stratified failure is
stated rather than dropped.

THE QUESTION THE NUMBERS CANNOT ANSWER ON THEIR OWN. Why are discordant loci
seven times more ambiguous? Either they sit in multi-copy families where several
partner proteins are genuinely near-equivalent - transposases and insertion
sequences are both pseudogene hotspots and large families - or they are shorter
and more degraded, so the top two bitscores land close together by chance. The
first is real biology and no threshold removes it; the second is a power
problem. The plan declares both tests and what each outcome would mean, before
seeing either.

    python3 analysis_forward_ambiguity.py [--input T] [--out-json J]
"""
from __future__ import annotations

import argparse
import ast
import collections
import hashlib
import json
import pathlib
import subprocess
import sys
import tempfile
import types

import numpy as np
import pandas as pd
from scipy import stats

P = pathlib.Path(__file__).resolve().parent
R = P / "results"
OUT = R / "forward_ambiguity"
D = pathlib.Path.home() / "research/datasets"
REFDIR = D / "n2-refseq-reference"          # READ ONLY

PLAN = P / "FORWARD_AMBIGUITY_PLAN.md"
PLAN_SHA = "e626f6d4386c4e824bd0872ead39ddee0206ef30c60ad6753a644365edb49698"

CONSEQUENCE_TOL = 0.005
PROXY_THRESHOLD = 0.10          # the old criterion, scored per class for the record
PARALOGUE_REL_BITSCORE = 0.20   # runner-up counts as homologous above this share
IPW_COVS = ["log_span", "forward_bitscore", "log_candidates"]
THREADS = 4
FLOOR_START_GB = 10


def borrow(path):
    tree = ast.parse(path.read_text())
    keep = [n for n in tree.body
            if isinstance(n, (ast.Import, ast.ImportFrom, ast.FunctionDef,
                              ast.ClassDef))
            or (isinstance(n, ast.Assign)
                and all(isinstance(t, ast.Name) and t.id.isupper()
                        for t in n.targets))]
    m = types.ModuleType(path.stem + "_b")
    m.__dict__["__file__"] = str(path)
    exec(compile(ast.Module(body=keep, type_ignores=[]), str(path), "exec"),
         m.__dict__)
    return m


# the frozen weighting machinery, reused rather than rewritten
S = borrow(P / "analysis_rbh_stratified.py")


def host_free_gb():
    out = subprocess.run(["df", "--output=avail", "-BG", "/mnt/c"],
                         capture_output=True, text=True, check=True).stdout
    return int("".join(c for c in out.splitlines()[-1] if c.isdigit()))


def partner_faa(partner):
    faa, sent = REFDIR / f"{partner}.faa", REFDIR / f"{partner}.json"
    if not (faa.exists() and sent.exists()):
        return None
    try:
        if json.loads(sent.read_text()).get("version") != 2:
            return None
    except Exception:
        return None
    return faa


def read_fasta(p):
    out, name, buf = {}, None, []
    with open(p) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


def self_and_pair_scores(pairs, faa, workdir):
    """-> {(best, second): (pair_bitscore, best_self_bitscore)}.

    Aligns each ambiguous locus's best subject against its runner-up, and
    against itself, so the pair score can be read as a share of the self score.
    A runner-up that is a strong hit of the best is a paralogue, which is the
    plan's H1 test.
    """
    prot = read_fasta(faa)
    want = [(b, s) for b, s in pairs if b in prot and s in prot]
    if not want:
        return {}
    workdir.mkdir(parents=True, exist_ok=True)
    subj = {s for _, s in want} | {b for b, _ in want}
    dbf = workdir / "db.faa"
    with dbf.open("w") as fh:
        for k in sorted(subj):
            fh.write(f">{k}\n{prot[k]}\n")
    qf = workdir / "q.faa"
    with qf.open("w") as fh:
        for k in sorted({b for b, _ in want}):
            fh.write(f">{k}\n{prot[k]}\n")
    dmnd = workdir / "db.dmnd"
    subprocess.run(["diamond", "makedb", "--in", str(dbf), "-d", str(dmnd),
                    "--quiet", "--threads", str(THREADS)], check=True)
    out = workdir / "hits.tsv"
    subprocess.run(["diamond", "blastp", "-q", str(qf), "-d", str(dmnd),
                    "-o", str(out), "--quiet", "--threads", str(THREADS),
                    "--max-target-seqs", "500", "--evalue", "10",
                    "--outfmt", "6", "qseqid", "sseqid", "bitscore"],
                   check=True)
    best = {}
    for line in out.read_text().splitlines():
        f = line.split("\t")
        if len(f) >= 3:
            k = (f[0], f[1])
            v = float(f[2])
            if v > best.get(k, -1.0):
                best[k] = v
    res = {}
    for b, s in want:
        res[(b, s)] = (best.get((b, s)), best.get((b, b)))
    return res


def rate_table(df, by=None):
    rows = []
    groups = [("all", df)] if by is None else list(df.groupby(by))
    for key, g in groups:
        row = {"key": key}
        for kind in ("discordant", "concordant_intact"):
            h = g[g.kind == kind]
            amb = pd.to_numeric(h.forward_ambiguous, errors="coerce")
            row[f"n_{kind}"] = int(len(h))
            row[f"rate_{kind}"] = float(amb.mean()) if len(h) else None
        a, b = row.get("rate_discordant"), row.get("rate_concordant_intact")
        row["difference_pp"] = (100 * (a - b)
                                if a is not None and b is not None else None)
        row["ratio"] = (a / b if a is not None and b not in (None, 0) else None)
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input",
                    default="results/rbh3/per_locus_rbh_audited.tsv")
    ap.add_argument("--out-json", default="results/forward_ambiguity.json")
    ap.add_argument("--out-dir", default=None)   # isolate an expanded run
    a = ap.parse_args()
    global OUT
    if a.out_dir:
        OUT = P / a.out_dir

    got = hashlib.sha256(PLAN.read_bytes()).hexdigest()
    print(f"plan sha256: {got[:16]}…")
    if got != PLAN_SHA:
        sys.exit(f"FROZEN PLAN CHANGED\n  expected {PLAN_SHA}\n  found    {got}")
    free0 = host_free_gb()
    print(f"host C: free {free0} GB (from /mnt/c)")
    if free0 < FLOOR_START_GB:
        sys.exit(f"refusing to start: {free0} < {FLOOR_START_GB} GB")
    OUT.mkdir(parents=True, exist_ok=True)

    src = P / a.input
    in_sha = hashlib.sha256(src.read_bytes()).hexdigest()
    print(f"input {src.relative_to(P)}  sha256 {in_sha[:16]}…")
    df = pd.read_csv(src, sep="\t")
    df["log_span"] = np.log10(df.span.clip(lower=1))
    df["log_candidates"] = np.log10(
        pd.to_numeric(df.n_forward_candidates, errors="coerce").clip(lower=1))
    df["amb"] = pd.to_numeric(df.forward_ambiguous, errors="coerce")
    genus_of = dict(zip(df.accession, df.genus))
    print(f"{len(df)} loci, {df.accession.nunique()} genomes, "
          f"{df.genus.nunique()} genera")

    # ============================================= 1. rates, per class
    pooled = float(df.amb.mean(skipna=True))
    t_all = rate_table(df)
    r_dis = float(t_all.rate_discordant.iloc[0])
    r_con = float(t_all.rate_concordant_intact.iloc[0])
    print("\nPART 3 — RATES, reported with no threshold attached")
    print(f"  pooled                    {100 * pooled:6.2f} %   "
          f"(the old criterion scored this and passed)")
    print(f"  concordant intact         {100 * r_con:6.2f} %")
    print(f"  discordant                {100 * r_dis:6.2f} %   "
          f"ratio {r_dis / r_con:.1f}x")
    print(f"  PROXY CRITERION, per class, threshold {PROXY_THRESHOLD:.0%}: "
          f"{'PASS' if max(r_dis, r_con) <= PROXY_THRESHOLD else 'FAIL'}"
          f"  — it fails, and that is reported rather than dropped")

    for by, name in (("accession", "per_genome"), ("genus", "per_genus")):
        t = rate_table(df, by).sort_values("rate_discordant", ascending=False)
        t.to_csv(OUT / f"rates_{name}.tsv", sep="\t", index=False)
    t_gen = rate_table(df, "accession").sort_values("rate_discordant",
                                                    ascending=False)
    print("\n  per genome, worst first")
    for _, r in t_gen.iterrows():
        print(f"    {r.key:20s} discordant {100 * r.rate_discordant:6.2f} %  "
              f"baseline {100 * r.rate_concordant_intact:5.2f} %  "
              f"({r.n_discordant:.0f} / {r.n_concordant_intact:.0f} loci)")

    # ============================================= 2. consequence on G3
    arms, gaps_by_arm = [], {}
    for label, mask in (("A0 all scored loci  (TRUTH)", None),
                        ("A1 forward-unambiguous only", df.amb == 0),
                        ("A2 forward-AMBIGUOUS only", df.amb == 1)):
        g = {}
        for acc, sub in df.groupby("accession"):
            m = None if mask is None else mask.loc[sub.index]
            v = S.gap_unweighted(sub, mask=m)
            if v is not None:
                g[acc] = v
        gaps_by_arm[label] = g
        arms.append(S.combine(g, genus_of, label))

    # A3: IPW recovering the full universe from unambiguous loci
    g3, diags = {}, []
    for acc, sub in df.groupby("accession"):
        wm, ok = {}, True
        for kind in ("discordant", "concordant_intact"):
            cell = sub[sub.kind == kind].copy()
            cell["hiconf"] = (cell.amb == 0).astype(int)   # "retained" = clear
            w, d = S.ipw_weights(cell, covs=IPW_COVS)
            d.update({"accession": acc, "kind": kind})
            diags.append(d)
            if w is None:
                ok = False
                continue
            keep = cell[cell.hiconf == 1]
            wm[kind] = (keep.omega.to_numpy(float), w)
        if ok and len(wm) == 2:
            md = S.wmedian(*wm["discordant"])
            mb = S.wmedian(*wm["concordant_intact"])
            if md is not None and mb is not None:
                g3[acc] = md - mb
    gaps_by_arm["A3 IPW on unambiguous loci"] = g3
    arms.append(S.combine(g3, genus_of, "A3 IPW on unambiguous loci"))
    pd.DataFrame(diags).to_csv(OUT / "ipw_diagnostics.tsv", sep="\t",
                               index=False)

    print("\nPARTS 1 and 2 — CONSEQUENCE")
    for A in arms:
        if A.get("median_gap") is None:
            print(f"  {A['label']:34s} {A.get('note')}")
            continue
        ci = A.get("genome_cluster_ci95") or [None, None]
        line = (f"  {A['label']:34s} {A['median_gap']:+.4f}  "
                f"{A['n_positive']}/{A['n_tested']}  p {A['p']:<9.3g} "
                f"({A['n_genomes']} genomes)")
        if ci[0] is not None:
            line += f"  genome CI [{ci[0]:+.4f}, {ci[1]:+.4f}]"
        print(line)

    truth = arms[0]["median_gap"]
    a0, a1 = gaps_by_arm["A0 all scored loci  (TRUTH)"], \
        gaps_by_arm["A1 forward-unambiguous only"]
    shared = sorted(set(a0) & set(a1))
    same_set_truth = float(np.median([a0[k] for k in shared]))
    consequence = arms[1]["median_gap"] - same_set_truth
    flips = [k for k in shared if np.sign(a0[k]) != np.sign(a1[k])]
    print(f"\n  A1 minus A0 on the same {len(shared)} genomes: "
          f"{consequence:+.5f}  (tolerance {CONSEQUENCE_TOL})")
    print(f"  genomes changing sign between A0 and A1: {len(flips)} "
          f"{flips if flips else ''}")
    per_g = {k: a1[k] - a0[k] for k in shared}
    print(f"  per-genome change: median {np.median(list(per_g.values())):+.5f}, "
          f"max |change| {max(abs(v) for v in per_g.values()):.5f}")

    # ============================================= 3. mechanism
    print("\nPART 5 — MECHANISM")
    fit = df.dropna(subset=IPW_COVS + ["amb"]).copy()
    fit["is_dis"] = (fit.kind == "discordant").astype(float)
    y = fit.amb.to_numpy(float)
    b_un = S.logistic_fit(np.column_stack(
        [np.ones(len(fit)), fit.is_dis.to_numpy(float)]), y)
    Xa = np.column_stack([np.ones(len(fit)), fit.is_dis.to_numpy(float)]
                         + [fit[c].to_numpy(float) for c in IPW_COVS])
    b_adj = S.logistic_fit(Xa, y)
    mech = {}
    if b_un is not None and b_adj is not None:
        mech["class_logodds_unadjusted"] = float(b_un[1])
        mech["class_logodds_adjusted"] = float(b_adj[1])
        mech["class_or_unadjusted"] = float(np.exp(b_un[1]))
        mech["class_or_adjusted"] = float(np.exp(b_adj[1]))
        mech["adjusted_coefs"] = {c: float(v) for c, v in
                                  zip(["intercept", "is_discordant"] + IPW_COVS,
                                      b_adj)}
        print(f"  discordant odds ratio for ambiguity, unadjusted "
              f"{mech['class_or_unadjusted']:6.2f}")
        print(f"  the same, adjusted for span, bitscore, candidate count "
              f"{mech['class_or_adjusted']:6.2f}")
        for c in IPW_COVS:
            print(f"    coefficient on {c:20s} "
                  f"{mech['adjusted_coefs'][c]:+.3f}")

    cands = {}
    for kind, g in df.groupby("kind"):
        n = pd.to_numeric(g.n_forward_candidates, errors="coerce")
        cands[kind] = {"median": float(n.median()),
                       "p90": float(n.quantile(0.9)),
                       "median_ambiguous": float(
                           n[g.amb == 1].median()) if (g.amb == 1).any() else None}
        print(f"  forward candidates, {kind:18s} median {cands[kind]['median']:.0f}"
              f"  90th pct {cands[kind]['p90']:.0f}"
              f"  median among ambiguous "
              f"{cands[kind]['median_ambiguous']}")

    # runner-up homology, the H1 test
    para_rows = []
    amb = df[(df.amb == 1) & df.forward_second.notna()]
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        for (acc, partner), g in amb.groupby(["accession", "partner"]):
            faa = partner_faa(partner)
            if faa is None:
                continue
            pairs = list({(r.forward_subject, r.forward_second)
                          for r in g.itertuples()})
            sc = self_and_pair_scores(pairs, faa, td / acc)
            for r in g.itertuples():
                pair, self_ = sc.get((r.forward_subject, r.forward_second),
                                     (None, None))
                rel = (pair / self_ if pair and self_ else None)
                para_rows.append({
                    "accession": acc, "kind": r.kind,
                    "best": r.forward_subject, "second": r.forward_second,
                    "pair_bitscore": pair, "best_self_bitscore": self_,
                    "relative_bitscore": rel,
                    "homologous": (int(rel >= PARALOGUE_REL_BITSCORE)
                                   if rel is not None else None)})
    para = pd.DataFrame(para_rows)
    paralogy = {}
    if len(para):
        para.to_csv(OUT / "runner_up_homology.tsv", sep="\t", index=False)
        print(f"\n  runner-up homology on {len(para)} ambiguous loci "
              f"(threshold: pair bitscore >= "
              f"{PARALOGUE_REL_BITSCORE:.0%} of the best's self score)")
        for kind, g in para.groupby("kind"):
            h = pd.to_numeric(g.homologous, errors="coerce").dropna()
            rel = pd.to_numeric(g.relative_bitscore, errors="coerce").dropna()
            paralogy[kind] = {
                "n": int(len(g)), "n_decided": int(len(h)),
                "frac_homologous": float(h.mean()) if len(h) else None,
                "median_relative_bitscore": float(rel.median())
                if len(rel) else None}
            print(f"    {kind:18s} n {len(g):5d}  homologous "
                  f"{100 * (h.mean() if len(h) else float('nan')):5.1f} %  "
                  f"median relative bitscore "
                  f"{rel.median() if len(rel) else float('nan'):.3f}")

    # ============================================= verdicts
    gate = {
        "consequence": {"value": consequence, "tolerance": CONSEQUENCE_TOL,
                        "pass": bool(abs(consequence) <= CONSEQUENCE_TOL)},
        "no_genome_flips": {"value": len(flips), "pass": bool(not flips)},
    }
    proxy = {"threshold": PROXY_THRESHOLD,
             "rate_discordant": r_dis, "rate_concordant_intact": r_con,
             "pass": bool(max(r_dis, r_con) <= PROXY_THRESHOLD)}

    dis_or_survives = (mech.get("class_or_adjusted") is not None
                       and mech["class_or_adjusted"] >= 1.5)
    homol = paralogy.get("discordant", {}).get("frac_homologous")
    if dis_or_survives and homol is not None and homol >= 0.5:
        verdict = "H1_paralogy"
    elif not dis_or_survives and (homol is None or homol < 0.5):
        verdict = "H2_weak_signal"
    else:
        verdict = "neither_H1_nor_H2_asymmetry_unexplained"

    print("\nPLAN SECTION 5.1 — the declared reading of the mechanism")
    print(f"  {verdict}")
    print("\nPLAN SECTION 6 — gate")
    for k, v in gate.items():
        print(f"  {k:22s} {str(v['value'])[:9]:>9}  -> "
              f"{'PASS' if v['pass'] else 'FAIL'}")
    print(f"  {'proxy, per class':22s} "
          f"{100 * max(r_dis, r_con):8.2f} %  -> "
          f"{'PASS' if proxy['pass'] else 'FAIL'}  (reported, not the criterion)")
    all_pass = all(v["pass"] for v in gate.values())
    print(f"  overall, as declared: {'PASS' if all_pass else 'FAIL'}")

    res = {"note": "executes FORWARD_AMBIGUITY_PLAN.md; no omega recomputed",
           "plan_sha256": got, "input": str(src.relative_to(P)),
           "input_sha256": in_sha, "n_loci": int(len(df)),
           "rates": {"pooled": pooled, "discordant": r_dis,
                     "concordant_intact": r_con,
                     "ratio": r_dis / r_con if r_con else None},
           "proxy_criterion_per_class": proxy,
           "arms": arms,
           "consequence_same_genome_set": {
               "n_genomes": len(shared), "A0": same_set_truth,
               "A1": arms[1]["median_gap"], "difference": consequence,
               "genomes_changing_sign": flips,
               "per_genome_median_change": float(
                   np.median(list(per_g.values()))),
               "per_genome_max_abs_change": float(
                   max(abs(v) for v in per_g.values()))},
           "mechanism": mech, "candidate_counts": cands,
           "runner_up_homology": paralogy,
           "paralogue_rel_bitscore_threshold": PARALOGUE_REL_BITSCORE,
           "verdict_per_plan_section_5_1": verdict,
           "gate": gate, "all_pass": all_pass,
           "caveat":
               "The criterion was changed from a rate proxy to a measured "
               "consequence. That is legitimate only because omega is known "
               "for ambiguous loci, and the proxy's class-stratified failure "
               "is reported alongside the consequence everywhere.",
           "host_free_gb_start": free0, "host_free_gb_end": host_free_gb()}
    (P / a.out_json).write_text(json.dumps(res, indent=1, default=str,
                                           sort_keys=True))
    print(f"\nwrote {a.out_json} and {OUT}/*.tsv")


if __name__ == "__main__":
    main()
