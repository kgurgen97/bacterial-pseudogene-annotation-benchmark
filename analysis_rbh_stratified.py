#!/usr/bin/env python3
"""Genome-stratified handling of differential RBH retention.

Executes RBH_STRATIFIED_PLAN.md, whose SHA-256 is asserted below. The plan was
frozen before these estimators ran and is never edited because of a result.

THE ONE THING THAT MAKES THIS DIFFERENT FROM EVERY OTHER SELECTION CORRECTION.
Omega is known for every locus, including the ones RBH excluded: exclusion
removes a locus from an ARM, not from the data. So the quantity the correction
is supposed to recover - the full-universe gap S0 - is already measured, and the
correction can be CHECKED against it instead of justified by an assumption. If
inverse-probability weighting cannot recover a number we already know, it cannot
be trusted on a number we do not, and that is reported as a failure rather than
smoothed over.

Reporting only in the sense that matters: no omega is recomputed, no orthologue
is reassigned, no locus is reclassified. The input is the committed per-locus
evidence from analysis_rbh_audited.py.

    python3 analysis_rbh_stratified.py [--input T] [--out-json J]
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import pathlib
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy import stats

P = pathlib.Path(__file__).resolve().parent
R = P / "results"
OUT = R / "rbh_stratified"

PLAN = P / "RBH_STRATIFIED_PLAN.md"
PLAN_SHA = "1fadc6db0ee4310b8004a9c488deb262d9d0815687423200edafa3e7fc77d4ef"

COVS = ["log_span", "forward_coverage"]     # NOT pident: it is a mediator
MIN_LOCI = 3                                # the project's own gap rule
MIN_RETAINED_FOR_MODEL = 10
MIN_EXCLUDED_FOR_MODEL = 3
WEIGHT_CAP = 10.0
ESS_FLOOR = 0.30
RIDGE = 1e-6
IMBALANCE_CUT = 0.10                        # S4's restriction
CALIPER = 0.10                              # S3, per the plan's disclosure
RECOVERY_TOL = 0.005                        # the gate
FLOOR_START_GB, FLOOR_ABORT_GB = 10, 8
BOOT = 2000
SEED = 20261031


def host_free_gb():
    out = subprocess.run(["df", "--output=avail", "-BG", "/mnt/c"],
                         capture_output=True, text=True, check=True).stdout
    return int("".join(c for c in out.splitlines()[-1] if c.isdigit()))


# --------------------------------------------------------------- estimators
def wmedian(v, w):
    """Weighted median: the smallest value whose cumulative weight reaches half.

    With equal weights this is the lower of the two central order statistics for
    an even count, which is numpy's 'lower' convention rather than its default
    midpoint; the unit tests pin that down so the IPW and unweighted arms are
    never compared across two different definitions.
    """
    v = np.asarray(v, float)
    w = np.asarray(w, float)
    keep = np.isfinite(v) & np.isfinite(w) & (w > 0)
    v, w = v[keep], w[keep]
    if not len(v):
        return None
    o = np.argsort(v)
    v, w = v[o], w[o]
    c = np.cumsum(w)
    return float(v[np.searchsorted(c, c[-1] / 2.0)])


def logistic_fit(Xm, y, ridge=RIDGE, iters=50):
    """Newton-Raphson logistic regression. -> coefficients, or None.

    Ridge only on the slopes, never the intercept, and only for numerical
    stability; it is 1e-6 and cannot shrink anything materially.
    """
    Xm = np.asarray(Xm, float)
    y = np.asarray(y, float)
    n, p = Xm.shape
    b = np.zeros(p)
    pen = np.eye(p) * ridge
    pen[0, 0] = 0.0
    for _ in range(iters):
        eta = np.clip(Xm @ b, -30, 30)
        mu = 1.0 / (1.0 + np.exp(-eta))
        Wd = np.clip(mu * (1 - mu), 1e-9, None)
        g = Xm.T @ (y - mu) - pen @ b
        H = -(Xm.T * Wd) @ Xm - pen
        try:
            step = np.linalg.solve(H, g)
        except np.linalg.LinAlgError:
            return None
        b_new = b - step
        if not np.all(np.isfinite(b_new)):
            return None
        if np.max(np.abs(b_new - b)) < 1e-8:
            return b_new
        b = b_new
    return b


def ipw_weights(cell, covs=COVS):
    """-> (weights for the retained rows, diagnostics) within one cell.

    A cell is one (genome, locus class). Retention is modelled from the
    covariates the plan declares; the weight is 1/p-hat, stabilised to mean 1
    inside the cell and capped per the plan.
    """
    y = (cell.hiconf == 1).to_numpy(float)
    n_ret, n_exc = int(y.sum()), int((1 - y).sum())
    diag = {"n_retained": n_ret, "n_excluded": n_exc}
    if n_ret < MIN_RETAINED_FOR_MODEL or n_exc < MIN_EXCLUDED_FOR_MODEL:
        diag["status"] = "too_few_to_model"
        return None, diag
    Xm = np.column_stack([np.ones(len(cell))]
                         + [cell[c].to_numpy(float) for c in covs])
    if not np.all(np.isfinite(Xm)):
        diag["status"] = "non_finite_covariate"
        return None, diag
    b = logistic_fit(Xm, y)
    if b is None:
        diag["status"] = "fit_failed"
        return None, diag
    phat = 1.0 / (1.0 + np.exp(-np.clip(Xm @ b, -30, 30)))
    ret = y == 1
    p_ret = np.clip(phat[ret], 1e-6, 1.0)
    w = 1.0 / p_ret
    w = w / w.mean()                      # stabilise to mean 1
    n_clipped = int((w > WEIGHT_CAP).sum())
    w = np.minimum(w, WEIGHT_CAP)
    ess = float((w.sum() ** 2) / (w ** 2).sum())
    diag.update({"status": "ok", "max_weight": float(w.max()),
                 "n_clipped": n_clipped, "ess": ess,
                 "ess_frac": ess / n_ret,
                 "weight_dominated": bool(ess / n_ret < ESS_FLOOR),
                 "coef": {c: float(v) for c, v in
                          zip(["intercept"] + list(covs), b)},
                 "min_phat": float(p_ret.min())})
    return w, diag


def smd_w(a, wa, b, wb):
    """Weighted standardised mean difference, pooled weighted SD."""
    def mv(x, w):
        x, w = np.asarray(x, float), np.asarray(w, float)
        m = np.average(x, weights=w)
        return m, np.average((x - m) ** 2, weights=w)
    ma, va = mv(a, wa)
    mb, vb = mv(b, wb)
    s = np.sqrt((va + vb) / 2)
    return 0.0 if s == 0 else float((ma - mb) / s)


def combine(gaps, genus_of, label, boot=BOOT, seed=SEED):
    """Median of per-genome gaps, with genome- and genus-clustered intervals."""
    out = {"label": label, "n_genomes": len(gaps)}
    if not gaps:
        out["note"] = "no genome scoreable"
        return out
    v = list(gaps.values())
    nz = [x for x in v if x != 0]
    pos = sum(1 for x in nz if x > 0)
    out.update({"median_gap": float(np.median(v)),
                "n_positive": pos, "n_tested": len(nz),
                "p": float(stats.binomtest(pos, len(nz), 0.5).pvalue)
                if nz else None})
    g = np.random.default_rng(seed)
    for unit in ("genome", "genus"):
        buckets = collections.defaultdict(list)
        for acc, val in gaps.items():
            buckets[acc if unit == "genome" else genus_of.get(acc, "?")].append(val)
        keys = list(buckets)
        if len(keys) < 3:
            out[f"{unit}_cluster_ci95"] = [None, None]
            continue
        meds = [float(np.median([x for k in g.choice(keys, len(keys),
                                                     replace=True)
                                 for x in buckets[k]])) for _ in range(boot)]
        out[f"{unit}_cluster_ci95"] = [float(np.percentile(meds, 2.5)),
                                       float(np.percentile(meds, 97.5))]
    return out


def gap_unweighted(sub, col="omega", mask=None):
    s = sub if mask is None else sub[mask]
    a = s.loc[s.kind == "discordant", col].dropna()
    b = s.loc[s.kind == "concordant_intact", col].dropna()
    if len(a) < MIN_LOCI or len(b) < MIN_LOCI:
        return None
    return float(np.median(a) - np.median(b))


def match_within(sub, caliper=CALIPER, seed=SEED):
    """S3: 1:1 nearest neighbour on the declared covariates, with a caliper."""
    g = np.random.default_rng(seed)
    keep = sub[sub.hiconf == 1].dropna(subset=COVS + ["omega"])
    dis = keep[keep.kind == "discordant"]
    con = keep[keep.kind == "concordant_intact"]
    if len(dis) < MIN_LOCI or len(con) < MIN_LOCI:
        return None, 0, 0
    Z = keep[COVS]
    mu, sd = Z.mean(), Z.std().replace(0, 1)
    zd = ((dis[COVS] - mu) / sd).to_numpy()
    zc = ((con[COVS] - mu) / sd).to_numpy()
    taken, pairs = set(), []
    for i in g.permutation(len(zd)):
        d = np.sqrt(((zc - zd[i]) ** 2).sum(axis=1))
        if taken:
            d[list(taken)] = np.inf
        j = int(np.argmin(d))
        if not np.isfinite(d[j]) or d[j] > caliper:
            continue
        taken.add(j)
        pairs.append((float(dis.omega.iloc[i]), float(con.omega.iloc[j])))
    if len(pairs) < MIN_LOCI:
        return None, len(pairs), len(dis)
    a = np.array([p[0] for p in pairs])
    b = np.array([p[1] for p in pairs])
    return float(np.median(a) - np.median(b)), len(pairs), len(dis)


# ----------------------------------------------------------------- main
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input",
                    default="results/rbh3/per_locus_rbh_audited.tsv")
    ap.add_argument("--out-json", default="results/rbh_stratified.json")
    ap.add_argument("--out-dir", default=None)   # isolate an expanded run
    a = ap.parse_args()
    global OUT
    if a.out_dir:
        OUT = P / a.out_dir

    got = hashlib.sha256(PLAN.read_bytes()).hexdigest()
    print(f"plan sha256: {got[:16]}…")
    if got != PLAN_SHA:
        sys.exit(f"FROZEN PLAN CHANGED\n  expected {PLAN_SHA}\n  found    {got}\n"
                 "  a frozen plan is never edited because of a result; write a "
                 "new, separately labelled plan instead")
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
    genus_of = dict(zip(df.accession, df.genus))
    print(f"{len(df)} loci, {df.accession.nunique()} genomes, "
          f"{df.genus.nunique()} genera")

    # ---------------------------------------------- per genome, every arm
    rows, diags, s0, s1, s2, s3, s4 = [], [], {}, {}, {}, {}, {}
    for acc, sub in df.groupby("accession"):
        g0 = gap_unweighted(sub)
        g1 = gap_unweighted(sub, mask=(sub.hiconf == 1))
        if g0 is not None:
            s0[acc] = g0
        if g1 is not None:
            s1[acc] = g1

        # S2: IPW, fitted per class inside this genome
        wm, ok, cell_diag = {}, True, {}
        for kind in ("discordant", "concordant_intact"):
            cell = sub[sub.kind == kind]
            w, d = ipw_weights(cell)
            d.update({"accession": acc, "kind": kind})
            cell_diag[kind] = d
            diags.append(d)
            if w is None:
                ok = False
                continue
            ret = cell[cell.hiconf == 1]
            wm[kind] = (ret.omega.to_numpy(float), w, ret)
        if ok and len(wm) == 2:
            md = wmedian(*wm["discordant"][:2])
            mb = wmedian(*wm["concordant_intact"][:2])
            if md is not None and mb is not None:
                s2[acc] = md - mb

        g3, n_pairs, n_dis_elig = match_within(sub)
        if g3 is not None:
            s3[acc] = g3

        d_ret = float(sub.loc[sub.kind == "discordant", "hiconf"].mean())
        b_ret = float(sub.loc[sub.kind == "concordant_intact", "hiconf"].mean())
        imb = b_ret - d_ret
        if abs(imb) <= IMBALANCE_CUT and g1 is not None:
            s4[acc] = g1

        rows.append({
            "accession": acc, "genus": genus_of[acc],
            "n_dis": int((sub.kind == "discordant").sum()),
            "n_con": int((sub.kind == "concordant_intact").sum()),
            "ret_dis": d_ret, "ret_con": b_ret, "imbalance": imb,
            "S0_all": g0, "S1_retained": g1, "S2_ipw": s2.get(acc),
            "S3_matched": g3, "S3_pairs": n_pairs,
            "S3_dis_eligible": n_dis_elig,
            "in_S4": int(abs(imb) <= IMBALANCE_CUT),
            "ipw_ess_frac_dis": cell_diag.get("discordant", {}).get("ess_frac"),
            "ipw_ess_frac_con": cell_diag.get("concordant_intact", {})
                                         .get("ess_frac"),
            "ipw_weight_dominated": int(any(
                cell_diag.get(k, {}).get("weight_dominated") for k in cell_diag)),
        })

    per = pd.DataFrame(rows).sort_values("imbalance", ascending=False)
    per.to_csv(OUT / "per_genome.tsv", sep="\t", index=False)
    pd.DataFrame(diags).to_csv(OUT / "ipw_diagnostics.tsv", sep="\t",
                               index=False)

    arms = [combine(s0, genus_of, "S0 all scored loci  (TRUTH)"),
            combine(s1, genus_of, "S1 RBH-retained, unweighted (biased)"),
            combine(s2, genus_of, "S2 IPW on retained loci"),
            combine(s3, genus_of, f"S3 matched, caliper {CALIPER}"),
            combine(s4, genus_of, f"S4 genomes with imbalance <= "
                                  f"{IMBALANCE_CUT:.0%}")]

    print("\nPER-GENOME TABLE (sorted by retention imbalance)")
    print(f"  {'genome':20s} {'imb':>7} {'S0':>9} {'S1':>9} {'S2':>9} "
          f"{'S3':>9} {'ESSd':>6} S4")
    for _, r in per.iterrows():
        def f(v):
            return f"{v:+.4f}" if pd.notna(v) else "     —"
        print(f"  {r.accession:20s} {100 * r.imbalance:+6.1f} {f(r.S0_all)} "
              f"{f(r.S1_retained)} {f(r.S2_ipw)} {f(r.S3_matched)} "
              f"{(f'{100 * r.ipw_ess_frac_dis:5.1f}' if pd.notna(r.ipw_ess_frac_dis) else '    —')} "
              f"{'yes' if r.in_S4 else ''}")

    print("\nARMS")
    for A in arms:
        if A.get("median_gap") is None:
            print(f"  {A['label']:44s} {A.get('note')}")
            continue
        ci = A.get("genome_cluster_ci95") or [None, None]
        line = (f"  {A['label']:44s} {A['median_gap']:+.4f}  "
                f"{A['n_positive']}/{A['n_tested']}  p {A['p']:<9.3g} "
                f"({A['n_genomes']} genomes)")
        if ci[0] is not None:
            line += f"  genome CI [{ci[0]:+.4f}, {ci[1]:+.4f}]"
        print(line)

    truth = arms[0].get("median_gap")
    bias = {}
    for A in arms[1:]:
        if A.get("median_gap") is not None and truth is not None:
            bias[A["label"]] = A["median_gap"] - truth
    print(f"\nBIAS AGAINST S0 = {truth:+.4f} (the known full-universe gap)")
    for k, v in bias.items():
        print(f"  {k:44s} {v:+.4f}")

    # ---- ADDED AFTER THE PLAN WAS FROZEN, as a correction to its own gate.
    # S2 can only be fitted where a genome has excluded loci to model, so it
    # covers 9 of 12 genomes while S0 covers 12. The plan compares them
    # directly, which is comparing two different genome sets; that was an
    # oversight in the plan, and the plan is not edited to hide it. Both the
    # declared comparison and the genome-set-matched one are reported.
    matched_sets = {}
    for name, arm_gaps in (("S1", s1), ("S2", s2), ("S3", s3), ("S4", s4)):
        keys = sorted(set(arm_gaps) & set(s0))
        if len(keys) < 3:
            continue
        t = float(np.median([s0[k] for k in keys]))
        e = float(np.median([arm_gaps[k] for k in keys]))
        matched_sets[name] = {"n_genomes": len(keys), "S0_same_genomes": t,
                              "arm": e, "bias": e - t,
                              "genomes": keys}
    print("\nSAME-GENOME-SET COMPARISON (the plan compared different sets)")
    for name, m in matched_sets.items():
        print(f"  {name}  {m['n_genomes']:2d} genomes   "
              f"S0 on those {m['S0_same_genomes']:+.4f}   "
              f"arm {m['arm']:+.4f}   bias {m['bias']:+.4f}")

    # ---- paired, genome by genome: does IPW move S1 toward S0?
    shared = sorted(set(s0) & set(s1) & set(s2))
    paired = {"n_genomes": len(shared)}
    if len(shared) >= 5:
        d1 = np.array([abs(s1[k] - s0[k]) for k in shared])
        d2 = np.array([abs(s2[k] - s0[k]) for k in shared])
        closer = int((d2 < d1).sum())
        paired.update({
            "median_abs_error_S1": float(np.median(d1)),
            "median_abs_error_S2": float(np.median(d2)),
            "n_genomes_ipw_closer": closer,
            "sign_test_p": float(stats.binomtest(
                closer, len(shared), 0.5).pvalue)})
        print(f"\nPER-GENOME RECOVERY, {len(shared)} genomes")
        print(f"  median |S1 - S0| {paired['median_abs_error_S1']:.4f}   "
              f"median |S2 - S0| {paired['median_abs_error_S2']:.4f}")
        print(f"  IPW closer to the truth in {closer}/{len(shared)} genomes, "
              f"sign test p {paired['sign_test_p']:.3g}")

    # ---- the interpreting diagnostic: does the gap track the imbalance?
    sub = per.dropna(subset=["S1_retained"])
    rho = p_rho = None
    if len(sub) >= 5:
        sp = stats.spearmanr(sub.imbalance, sub.S1_retained)
        rho, p_rho = float(sp.statistic), float(sp.pvalue)
        print(f"\nIMBALANCE vs RETAINED GAP across {len(sub)} genomes: "
              f"Spearman rho {rho:+.3f}, p {p_rho:.3g}")

    # ---- diagnostics for the weighting, two of them, measuring DIFFERENT things
    #
    # The plan asked for covariate balance between the two classes. That is the
    # right target for MATCHING and the WRONG target for IPW: inverse-
    # probability weighting is meant to restore each class's own full-sample
    # covariate distribution, not to make the two classes resemble each other.
    # Between-class differences in span and coverage are real properties of
    # discordant and intact genes, not confounding to be weighted away. Both
    # diagnostics are therefore reported, with the IPW-appropriate one added
    # after the plan was frozen and labelled as such.
    bal, recov = [], []
    for acc, g in df.groupby("accession"):
        cells = {}
        for kind in ("discordant", "concordant_intact"):
            cell = g[g.kind == kind]
            w, _ = ipw_weights(cell)
            if w is None:
                continue
            ret = cell[cell.hiconf == 1]
            cells[kind] = (cell, ret, w)
            for c in COVS:
                # does weighting restore THIS class's full-cell distribution?
                recov.append({
                    "accession": acc, "kind": kind, "covariate": c,
                    "smd_retained_vs_full": smd_w(
                        ret[c], np.ones(len(ret)), cell[c],
                        np.ones(len(cell))),
                    "smd_ipw_vs_full": smd_w(
                        ret[c], w, cell[c], np.ones(len(cell)))})
        if len(cells) == 2:
            _, rd, wd = cells["discordant"]
            _, rc, wc = cells["concordant_intact"]
            for c in COVS:
                bal.append({
                    "accession": acc, "covariate": c,
                    "smd_unweighted": smd_w(rd[c], np.ones(len(rd)),
                                            rc[c], np.ones(len(rc))),
                    "smd_ipw": smd_w(rd[c], wd, rc[c], wc)})

    balance, recovery = {}, {}
    if bal:
        bt = pd.DataFrame(bal)
        bt.to_csv(OUT / "ipw_balance.tsv", sep="\t", index=False)
        for c, s in bt.groupby("covariate"):
            balance[c] = {
                "median_abs_smd_unweighted": float(
                    np.nanmedian(np.abs(s.smd_unweighted))),
                "median_abs_smd_ipw": float(np.nanmedian(np.abs(s.smd_ipw)))}
        print("\nBETWEEN-CLASS BALANCE (the plan's diagnostic; the wrong target "
              "for IPW)")
        for c, v in balance.items():
            print(f"  {c:20s} unweighted {v['median_abs_smd_unweighted']:.3f}"
                  f"   IPW {v['median_abs_smd_ipw']:.3f}")
    if recov:
        rt = pd.DataFrame(recov)
        rt.to_csv(OUT / "ipw_distribution_recovery.tsv", sep="\t", index=False)
        print("\nDISTRIBUTION RECOVERY (added after the plan; the right target "
              "for IPW)")
        print("  median |SMD| of the retained sample against its own FULL cell")
        for (kind, c), s in rt.groupby(["kind", "covariate"]):
            u = float(np.nanmedian(np.abs(s.smd_retained_vs_full)))
            w = float(np.nanmedian(np.abs(s.smd_ipw_vs_full)))
            recovery[f"{kind}|{c}"] = {"unweighted": u, "ipw": w,
                                       "improved": bool(w < u)}
            print(f"  {kind:18s} {c:18s} unweighted {u:.3f}  IPW {w:.3f}  "
                  f"{'better' if w < u else 'WORSE'}")

    # ---------------------------------------------- the plan's own verdict
    n_dom = int(per.ipw_weight_dominated.sum())
    dom_frac = n_dom / len(per) if len(per) else 1.0
    recovered = (truth is not None and arms[2].get("median_gap") is not None
                 and abs(arms[2]["median_gap"] - truth) <= RECOVERY_TOL)
    gaps_all = [A["median_gap"] for A in arms[1:]
                if A.get("median_gap") is not None]
    sign_stable = bool(gaps_all and (all(v > 0 for v in gaps_all)
                                     or all(v < 0 for v in gaps_all)))
    spread_ok = bool(gaps_all and (max(gaps_all) - min(gaps_all)) <= 0.010)

    verdict = ("S2_recovers_S0" if recovered else
               "S2_between_S1_and_S0" if (
                   truth is not None and arms[2].get("median_gap") is not None
                   and arms[1].get("median_gap") is not None
                   and min(truth, arms[1]["median_gap"])
                   <= arms[2]["median_gap"]
                   <= max(truth, arms[1]["median_gap"]))
               else "S2_fails_to_recover_S0")

    gate = {
        "retention_imbalance_handled": {
            "recovery_within_tol": recovered,
            "tolerance": RECOVERY_TOL,
            "weight_dominated_genome_frac": dom_frac,
            "weight_dominated_threshold": 0.25,
            "pass": bool(recovered and dom_frac <= 0.25)},
        "sign_stable_across_arms": {"value": sign_stable, "pass": sign_stable},
        "arms_within_0.010": {"value": spread_ok, "pass": spread_ok},
    }
    print("\nPLAN SECTION 4 — the declared reading of this result")
    print(f"  {verdict}")
    print("\nPLAN SECTION 6 — gate")
    for k, v in gate.items():
        print(f"  {k:34s} -> {'PASS' if v['pass'] else 'FAIL'}")
    all_pass = all(v["pass"] for v in gate.values())
    print(f"  overall: {'PASS' if all_pass else 'FAIL'}")

    res = {"note": "executes RBH_STRATIFIED_PLAN.md; no omega recomputed",
           "plan_sha256": got, "input": str(src.relative_to(P)),
           "input_sha256": in_sha,
           "n_loci": int(len(df)), "n_genomes": int(df.accession.nunique()),
           "params": {"covariates": COVS, "weight_cap": WEIGHT_CAP,
                      "ess_floor": ESS_FLOOR, "caliper": CALIPER,
                      "imbalance_cut": IMBALANCE_CUT,
                      "recovery_tol": RECOVERY_TOL, "boot": BOOT,
                      "seed": SEED},
           "arms": arms, "bias_against_truth": bias,
           "same_genome_set_comparison": matched_sets,
           "paired_recovery": paired,
           "imbalance_vs_gap": {"spearman_rho": rho, "p": p_rho},
           "between_class_balance_plan_diagnostic": balance,
           "distribution_recovery_added_after_plan": recovery,
           "n_weight_dominated_genomes": n_dom,
           "verdict_per_plan_section_4": verdict,
           "gate": gate, "all_pass": all_pass,
           "selection_warning":
               "Recovering S0 from retained loci shows that retention is "
               "explainable by span and coverage. It does not show that the "
               "orthologue assignments RBH kept are correct, and it is not "
               "biological corroboration of G3.",
           "host_free_gb_start": free0, "host_free_gb_end": host_free_gb()}
    (P / a.out_json).write_text(json.dumps(res, indent=1, default=str,
                                           sort_keys=True))
    print(f"\nwrote {a.out_json}")
    print(f"wrote {OUT / 'per_genome.tsv'}, {OUT / 'ipw_diagnostics.tsv'}, "
          f"{OUT / 'ipw_balance.tsv'}")


if __name__ == "__main__":
    main()
