#!/usr/bin/env python3
"""Do the configurations differ in effect MAGNITUDE, or only in verdict?

Executes CALLER_EFFECT_CONTRAST_PLAN.md, whose SHA-256 is asserted below.

THE ERROR THIS CORRECTS. The manuscript concluded that annotation choice
changes a downstream result because PGAP and Pseudofinder gave p < 0.05 while
Bakta and DFAST did not. Two correlations on opposite sides of an arbitrary
threshold can be statistically indistinguishable from each other; comparing
verdicts is not comparing effects. The contrast the claim requires is the
PAIRED difference of correlations on the same genomes, with an interval.

No annotation is re-run, no genome downloaded, no dN/dS recomputed. Input is
the existing 363-genome re-test table.

    python3 analysis_caller_effect_contrast.py
"""
from __future__ import annotations

import collections
import hashlib
import itertools
import json
import pathlib
import re
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy import stats

P = pathlib.Path(__file__).resolve().parent
R = P / "results"
OUT = R / "caller_effect_contrast"

PLAN = P / "CALLER_EFFECT_CONTRAST_PLAN.md"
PLAN_SHA = "0d1d64582b5da461f7f0fa8387d6d386ed7884f58785f3e5a2d717bbd1aa3870"

B = 5000
SEED = 20261101
ALPHA = 0.05
CALLERS = ["PGAP", "Pseudofinder", "Bakta", "DFAST"]
# Pseudofinder post-processes PGAP output: pairs involving it are not two
# independent measurements and are quarantined from the headline reading.
DEPENDENT = {"Pseudofinder"}
METRICS = {"fraction": {"PGAP": "frac_pgap", "Pseudofinder": "frac_pfinder",
                        "Bakta": "frac_bakta", "DFAST": "frac_dfast"},
           "count": {"PGAP": "cnt_pgap", "Pseudofinder": "cnt_pfinder",
                     "Bakta": "cnt_bakta", "DFAST": "cnt_dfast"},
           "per Mb": {"PGAP": "per_mb_pgap", "Pseudofinder": "per_mb_pfinder",
                      "Bakta": "per_mb_bakta", "DFAST": "per_mb_dfast"}}
PRIMARY_METRIC = "fraction"


def host_free_gb():
    out = subprocess.run(["df", "--output=avail", "-BG", "/mnt/c"],
                         capture_output=True, text=True, check=True).stdout
    return int("".join(c for c in out.splitlines()[-1] if c.isdigit()))


def genus_of(organism):
    """Same rule analysis_dnds.py uses: first word, skipping 'Candidatus'."""
    if not isinstance(organism, str):
        return None
    w = re.sub(r"[\[\]]", "", organism).split()
    if w and w[0] == "Candidatus":
        w = w[1:]
    return w[0] if w else None


def spearman(x, y):
    if len(x) < 4:
        return np.nan
    return float(stats.spearmanr(x, y).statistic)


def partial_spearman(x, y, z):
    """Spearman of x,y controlling z, via ranks and residualisation."""
    if len(x) < 5:
        return np.nan
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))
    def resid(a, b):
        b1 = np.column_stack([np.ones(len(b)), b])
        coef, *_ = np.linalg.lstsq(b1, a, rcond=None)
        return a - b1 @ coef
    ex, ey = resid(rx, rz), resid(ry, rz)
    if ex.std() == 0 or ey.std() == 0:
        return np.nan
    return float(np.corrcoef(ex, ey)[0, 1])


def holm(pvals):
    """Holm-Bonferroni adjusted p-values, order preserved."""
    idx = np.argsort(pvals)
    m = len(pvals)
    adj = np.empty(m)
    run = 0.0
    for rank, i in enumerate(idx):
        val = (m - rank) * pvals[i]
        run = max(run, val)
        adj[i] = min(1.0, run)
    return adj


def main() -> None:
    got = hashlib.sha256(PLAN.read_bytes()).hexdigest()
    print(f"plan sha256: {got[:16]}…")
    if got != PLAN_SHA:
        sys.exit(f"FROZEN PLAN CHANGED\n  expected {PLAN_SHA}\n  found    {got}")
    free = host_free_gb()
    print(f"host C: free {free} GB (from /mnt/c)")
    OUT.mkdir(parents=True, exist_ok=True)

    lit = pd.read_csv(R / "literature_retest.tsv", sep="\t")
    # DEVIATION D1 from the frozen plan, section 2. The plan names
    # results/corpus523.tsv and parsing genus from its organism string. That
    # file carries `organism` for only the 160 EXPANSION genomes; for the 363
    # genomes of this re-test it is empty, so the plan as written yields zero
    # eligible genomes. results/corpus_traits.tsv covers exactly these 363 and
    # carries a parsed `genus` column already. The plan is NOT edited; the
    # deviation is recorded here and in the output. This is the same
    # lineage-source mismatch PHYLO_RESULTS.md documented for its own plan.
    traits = pd.read_csv(R / "corpus_traits.tsv", sep="\t")[["accession",
                                                             "genus"]]
    df = lit.merge(traits, on="accession", how="left")
    n_from_traits = int(df.genus.notna().sum())
    # fall back to parsing corpus523's organism for anything traits misses
    if df.genus.isna().any():
        corpus = pd.read_csv(R / "corpus523.tsv", sep="\t")[["accession",
                                                             "organism"]]
        df = df.merge(corpus, on="accession", how="left")
        df["genus"] = df.genus.fillna(df.organism.map(genus_of))
    deviation = {
        "id": "D1",
        "plan_section": "2",
        "plan_said": "genus parsed from results/corpus523.tsv organism string",
        "problem": "that column is populated for the 160 expansion genomes "
                   "only; for these 363 it is empty, so the plan as written "
                   "yields zero eligible genomes",
        "done_instead": "genus taken from results/corpus_traits.tsv, which "
                        "covers exactly these 363 and carries a parsed genus "
                        "column; corpus523 organism parsing kept as a "
                        "fallback",
        "n_resolved_from_corpus_traits": n_from_traits}
    print(f"DEVIATION D1: genus from corpus_traits.tsv "
          f"({n_from_traits}/{len(df)} resolved); the plan named corpus523, "
          f"whose organism column is empty for these genomes")
    need = ["size_bp", "gc", "genus"] + [c for m in METRICS.values()
                                         for c in m.values()]
    before = len(df)
    df = df.dropna(subset=need)
    print(f"{before} genomes in the table, {len(df)} eligible "
          f"(all four measures, size, GC and a resolvable genus)")
    print(f"{df.genus.nunique()} genera; largest holds "
          f"{df.genus.value_counts().iloc[0]} genomes")
    in_sha = hashlib.sha256((R / "literature_retest.tsv").read_bytes()).hexdigest()

    size = df.size_bp.to_numpy(float)
    gc = df.gc.to_numpy(float)
    genera = df.genus.to_numpy()
    by_genus = collections.defaultdict(list)
    for i, g in enumerate(genera):
        by_genus[g].append(i)
    gkeys = list(by_genus)

    results = {"note": "paired contrast of correlations between annotation "
                       "configurations, same genomes",
               "plan_sha256": got, "input_sha256": in_sha,
               "n_genomes": int(len(df)), "n_genera": int(len(gkeys)),
               "B": B, "seed": SEED, "deviations": [deviation],
               "metrics": {}}

    for metric, cols in METRICS.items():
        vals = {c: df[cols[c]].to_numpy(float) for c in CALLERS}

        # ---- verdict table: each configuration on its own
        verdict = []
        for c in CALLERS:
            r = spearman(vals[c], size)
            pr = float(stats.spearmanr(vals[c], size).pvalue)
            pa = partial_spearman(vals[c], size, gc)
            verdict.append({"metric": metric, "caller": c, "rho": r, "p": pr,
                            "significant": bool(pr < ALPHA),
                            "partial_rho": pa})
        vt = pd.DataFrame(verdict)

        # ---- paired cluster bootstrap over genera
        rng = np.random.default_rng(SEED)
        boot_rho = {c: np.empty(B) for c in CALLERS}
        boot_par = {c: np.empty(B) for c in CALLERS}
        for b in range(B):
            pick = rng.choice(gkeys, len(gkeys), replace=True)
            idx = np.concatenate([by_genus[g] for g in pick])
            s, z = size[idx], gc[idx]
            for c in CALLERS:
                v = vals[c][idx]
                boot_rho[c][b] = spearman(v, s)
                boot_par[c][b] = partial_spearman(v, s, z)

        rows = []
        for a, bq in itertools.combinations(CALLERS, 2):
            for kind, store in (("rho", boot_rho), ("partial_rho", boot_par)):
                d = store[a] - store[bq]
                d = d[np.isfinite(d)]
                if not len(d):
                    continue
                point = (vt.set_index("caller").loc[a, kind]
                         - vt.set_index("caller").loc[bq, kind])
                lo, hi = np.percentile(d, [2.5, 97.5])
                p = 2 * min((d <= 0).mean(), (d >= 0).mean())
                p = max(p, 1.0 / B)
                rows.append({
                    "metric": metric, "statistic": kind, "a": a, "b": bq,
                    "delta": float(point), "lo": float(lo), "hi": float(hi),
                    "p_boot": float(p),
                    "independent": bool(a not in DEPENDENT
                                        and bq not in DEPENDENT),
                    "ci_excludes_zero": bool(lo > 0 or hi < 0)})
        ct = pd.DataFrame(rows)
        # Holm within metric and statistic, separately for the two blocks
        ct["p_holm"] = np.nan
        for (stat, indep), grp in ct.groupby(["statistic", "independent"]):
            ct.loc[grp.index, "p_holm"] = holm(grp.p_boot.to_numpy())
        ct["significant_adj"] = ct.p_holm < ALPHA

        # verdict-vs-magnitude labelling, declared in the plan
        sig = dict(zip(vt.caller, vt.significant))
        ct["verdict_differs"] = [sig[r.a] != sig[r.b] for r in ct.itertuples()]
        ct["label"] = np.where(
            ct.ci_excludes_zero, "magnitudes differ",
            np.where(ct.verdict_differs,
                     "verdict differs, magnitudes not distinguishable",
                     "no difference detected"))

        vt.to_csv(OUT / f"verdicts_{metric.replace(' ', '_')}.tsv", sep="\t",
                  index=False)
        ct.to_csv(OUT / f"contrasts_{metric.replace(' ', '_')}.tsv", sep="\t",
                  index=False)

        print(f"\n{'=' * 74}\nMETRIC: {metric}"
              f"{'   (PRIMARY)' if metric == PRIMARY_METRIC else ''}")
        print("=" * 74)
        print("  VERDICT TABLE — what the manuscript currently reports")
        for _, r in vt.iterrows():
            print(f"    {r.caller:14s} rho {r.rho:+.4f}  p {r.p:<10.3g}"
                  f"  {'SIGNIFICANT' if r.significant else 'n.s.':12s}"
                  f"  partial {r.partial_rho:+.4f}")
        print("\n  CONTRAST TABLE — what the claim requires (rho difference)")
        for blk, lab in ((True, "INDEPENDENT pairs (headline)"),
                         (False, "pairs involving Pseudofinder "
                                 "(NOT independent of PGAP)")):
            sel = ct[(ct.statistic == "rho") & (ct.independent == blk)]
            print(f"    {lab}")
            for _, r in sel.iterrows():
                print(f"      {r.a:12s} - {r.b:12s} "
                      f"delta {r.delta:+.4f}  CI [{r.lo:+.4f}, {r.hi:+.4f}]"
                      f"  p {r.p_boot:<8.4f} Holm {r.p_holm:<8.4f}"
                      f"  {r.label}")
        results["metrics"][metric] = {
            "verdicts": vt.to_dict("records"),
            "contrasts": ct.to_dict("records")}

    # ---------------------------------------------- the declared reading
    prim = pd.DataFrame(results["metrics"][PRIMARY_METRIC]["contrasts"])
    ind = prim[(prim.statistic == "rho") & (prim.independent)]
    any_adj = bool((ind.p_holm < ALPHA).any())
    any_raw = bool(ind.ci_excludes_zero.any())
    if any_adj:
        reading = "magnitudes_differ_after_adjustment"
        says = ("configurations differ in effect MAGNITUDE for the pairs "
                "named; the claim stands in that form")
    elif any_raw:
        reading = "suggestive_unadjusted_only"
        says = ("evidence is SUGGESTIVE only; report unadjusted and adjusted "
                "side by side and do not claim a magnitude difference")
    else:
        reading = "retreat_verdict_only"
        says = ("RETREAT: the verdict differs while the magnitudes are NOT "
                "distinguishable. Title, abstract, section 9.3 and the "
                "Conclusions must be corrected")
    print(f"\n{'=' * 74}")
    print(f"PLAN SECTION 7 — declared reading on the primary metric "
          f"({PRIMARY_METRIC}, independent pairs)")
    print(f"  {reading}")
    print(f"  {says}")
    results["primary_metric"] = PRIMARY_METRIC
    results["declared_reading"] = reading
    results["manuscript_must_say"] = says
    results["n_independent_pairs_ci_excludes_zero"] = int(
        ind.ci_excludes_zero.sum())
    results["n_independent_pairs_holm_significant"] = int(
        (ind.p_holm < ALPHA).sum())

    (R / "caller_effect_contrast.json").write_text(
        json.dumps(results, indent=1, default=str, sort_keys=True))
    print(f"\nwrote results/caller_effect_contrast.json and {OUT}/*.tsv")


if __name__ == "__main__":
    main()
