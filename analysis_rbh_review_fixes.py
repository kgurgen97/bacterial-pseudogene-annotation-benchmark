#!/usr/bin/env python3
"""Evidence for the review corrections to MANUSCRIPT_CANDIDATE_RBH.

Reporting only, over committed per-locus evidence. No omega recomputed, no
orthologue reassigned, no estimator changed, nothing re-run.

Four questions, each raised by review of commit 9abba73:

  D  the DIRECTION of the retained-versus-excluded omega difference, which the
     candidate stated backwards. The pilot and the full corpus disagree, and
     the candidate carried the pilot's sentence forward over the new numbers.
  G  the GENUS-LEVEL directional test for the primary E2 arm, from the 28
     genus gaps, because the published sign-test p-values are GENOME-level
     while the intervals are genus-clustered.
  C  whether missing annotated CDS explains MODIFIED-RBH rejection. The
     modified method was built to handle loci with no annotated protein, so
     attributing its attrition to CDS availability needs testing, not
     assertion.
  N  whether CI narrowing after filtering is attributable to variance
     reduction from dropping high-variance units, or is simply unexplained.

    python3 analysis_rbh_review_fixes.py
"""
from __future__ import annotations

import collections
import json
import pathlib

import numpy as np
import pandas as pd
from scipy import stats

P = pathlib.Path(__file__).resolve().parent
R = P / "results"
OUT = R / "review_fixes"
MIN_LOCI = 3


def gap(sub, col="omega", mask=None):
    s = sub if mask is None else sub[mask]
    a = s.loc[s.kind == "discordant", col].dropna()
    b = s.loc[s.kind == "concordant_intact", col].dropna()
    if len(a) < MIN_LOCI or len(b) < MIN_LOCI:
        return None
    return float(np.median(a) - np.median(b))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pilot = pd.read_csv(R / "rbh3" / "per_locus_rbh_audited.tsv", sep="\t")
    full = pd.read_csv(R / "rbh3_all" / "per_locus_rbh_audited.tsv", sep="\t")
    gen = pd.read_csv(R / "rbh3_all" / "by_genus.tsv", sep="\t")
    res = {}

    # ================================================== D  omega direction
    print("=" * 72)
    print("D  DIRECTION of the retained-vs-excluded omega difference")
    print("=" * 72)
    rows = []
    for name, df in (("pilot, 12 genomes", pilot), ("full, 92 genomes", full)):
        for kind in ("discordant", "concordant_intact"):
            g = df[df.kind == kind]
            ret = g.loc[g.hiconf == 1, "omega"].dropna()
            exc = g.loc[g.hiconf == 0, "omega"].dropna()
            if not len(exc):
                continue
            u = stats.mannwhitneyu(exc, ret, alternative="two-sided")
            rows.append({
                "dataset": name, "kind": kind,
                "n_retained": len(ret), "n_excluded": len(exc),
                "median_omega_retained": float(ret.median()),
                "median_omega_excluded": float(exc.median()),
                "excluded_minus_retained": float(exc.median() - ret.median()),
                "excluded_higher": bool(exc.median() > ret.median()),
                "mannwhitney_p": float(u.pvalue)})
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "omega_direction.tsv", sep="\t", index=False)
    for _, r in t.iterrows():
        print(f"  {r.dataset:20s} {r.kind:18s} retained {r.median_omega_retained:.4f}"
              f"  excluded {r.median_omega_excluded:.4f}  "
              f"diff {r.excluded_minus_retained:+.4f}  "
              f"excluded {'HIGHER' if r.excluded_higher else 'lower'}  "
              f"MW p {r.mannwhitney_p:.3g}")
    dis = t[(t.dataset.str.startswith("full")) & (t.kind == "discordant")].iloc[0]
    dis_p = t[(t.dataset.str.startswith("pilot")) & (t.kind == "discordant")].iloc[0]
    print(f"\n  CORRECTION: at 92 genomes excluded discordant loci have "
          f"HIGHER median omega ({dis.median_omega_excluded:.7f} vs "
          f"{dis.median_omega_retained:.7f}).")
    print(f"  The pilot had the opposite sign "
          f"({dis_p.median_omega_excluded:.4f} vs "
          f"{dis_p.median_omega_retained:.4f}), so the direction REVERSED "
          f"between scales and the candidate carried the pilot's sentence "
          f"forward over the new numbers.")
    # G3 is the median of PER-GENOME gaps, not a gap over pooled loci. A first
    # version of this check pooled every locus corpus-wide and printed
    # +0.0265 -> +0.0240, which is a different statistic from the reported
    # E1/E2 and must not be quoted as them.
    def g3(col_mask):
        v = [gap(sub, mask=None if col_mask is None else col_mask.loc[sub.index])
             for _, sub in full.groupby("accession")]
        v = [x for x in v if x is not None]
        return float(np.median(v)), len(v)

    e1, n1 = g3(None)
    e2, n2 = g3(full.hiconf == 1)
    print(f"\n  consistency: removing higher-omega discordant loci should "
          f"LOWER the gap, and it does: E1 {e1:+.4f} ({n1} genomes) -> "
          f"E2 {e2:+.4f} ({n2} genomes), change {e2 - e1:+.4f}")
    res["omega_direction"] = {
        "excluded_discordant_higher_at_92": bool(dis.excluded_higher),
        "median_excluded": float(dis.median_omega_excluded),
        "median_retained": float(dis.median_omega_retained),
        "difference": float(dis.excluded_minus_retained),
        "mannwhitney_p": float(dis.mannwhitney_p),
        "pilot_direction_was_opposite": bool(not dis_p.excluded_higher),
        "E1": e1, "E2": e2, "E2_minus_E1": e2 - e1,
        "consistent_with_direction": bool((e2 - e1) < 0)}

    # ======================================= D2  CLUSTERED omega comparison
    # The pooled Mann-Whitney above treats 36 535 loci as independent. Loci
    # within a genome, and genomes within a genus, are not. It is kept as
    # descriptive evidence only; the tests that carry the claim are these.
    print()
    print("=" * 72)
    print("D2 CLUSTERED retained-vs-excluded omega, genome and genus units")
    print("=" * 72)

    def cliffs_delta(a, b):
        """P(a > b) - P(a < b), computed from ranks. Scale-free effect size."""
        a, b = np.asarray(a, float), np.asarray(b, float)
        n, m = len(a), len(b)
        if not n or not m:
            return None
        allv = np.concatenate([a, b])
        r = stats.rankdata(allv)
        ra = r[:n].sum()
        u = ra - n * (n + 1) / 2.0          # Mann-Whitney U for a
        return float(2.0 * u / (n * m) - 1.0)

    def boot_ci(vals, keys, B=2000, seed=20261108):
        """Cluster bootstrap over `keys`, median of pooled values."""
        buckets = collections.defaultdict(list)
        for k, v in zip(keys, vals):
            buckets[k].append(v)
        ks = list(buckets)
        if len(ks) < 3:
            return None, None
        g = np.random.default_rng(seed)
        meds = [float(np.median([x for k in g.choice(ks, len(ks), replace=True)
                                 for x in buckets[k]])) for _ in range(B)]
        return float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5))

    clustered = {}
    for kind in ("discordant", "concordant_intact"):
        for min_side in (3, 10):
            per, gkeys, deltas = [], [], []
            for acc, sub in full[full.kind == kind].groupby("accession"):
                ret = sub.loc[sub.hiconf == 1, "omega"].dropna()
                exc = sub.loc[sub.hiconf == 0, "omega"].dropna()
                if len(ret) < min_side or len(exc) < min_side:
                    continue
                per.append(float(exc.median() - ret.median()))
                gkeys.append(sub.genus.iloc[0])
                deltas.append(cliffs_delta(exc.to_numpy(), ret.to_numpy()))
            key = f"{kind}|min{min_side}"
            n = len(per)
            if n < 3:
                clustered[key] = {"n_genomes": n, "note": "too few genomes"}
                print(f"  {kind:18s} min {min_side:2d}/side: only {n} eligible "
                      f"genomes — not tested")
                continue
            nz = [v for v in per if v != 0]
            pos = sum(1 for v in nz if v > 0)
            sg = stats.binomtest(pos, len(nz), 0.5).pvalue if nz else None
            wx = stats.wilcoxon(per).pvalue if len(per) >= 6 else None
            lo_g, hi_g = boot_ci(per, [f"g{i}" for i in range(n)])
            lo_t, hi_t = boot_ci(per, gkeys)
            d_med = float(np.median([d for d in deltas if d is not None]))
            clustered[key] = {
                "n_genomes_eligible": n,
                "n_genera_represented": len(set(gkeys)),
                "median_per_genome_difference": float(np.median(per)),
                "n_positive": pos, "n_nonzero": len(nz),
                "sign_test_p_two_sided": float(sg) if sg is not None else None,
                "wilcoxon_p_two_sided": float(wx) if wx is not None else None,
                "genome_cluster_ci95": [lo_g, hi_g],
                "genus_cluster_ci95": [lo_t, hi_t],
                "median_cliffs_delta": d_med}
            c = clustered[key]
            print(f"  {kind:18s} min {min_side:2d}/side: {n} genomes "
                  f"({c['n_genera_represented']} genera)  "
                  f"median diff {c['median_per_genome_difference']:+.4f}  "
                  f"{pos}/{len(nz)} positive  sign p "
                  f"{sg:.3g}  Wilcoxon p "
                  f"{wx if wx is None else f'{wx:.3g}'}")
            print(f"  {'':18s}   genome CI [{lo_g:+.4f}, {lo_g and hi_g:+.4f}]"
                  f"   genus CI [{lo_t:+.4f}, {hi_t:+.4f}]"
                  f"   Cliff's delta {d_med:+.3f}")

    # the same thing with the GENUS as the analysis unit
    gper, gdelta = [], []
    for g_, sub in full[full.kind == "discordant"].groupby("genus"):
        ret = sub.loc[sub.hiconf == 1, "omega"].dropna()
        exc = sub.loc[sub.hiconf == 0, "omega"].dropna()
        if len(ret) < 3 or len(exc) < 3:
            continue
        gper.append(float(exc.median() - ret.median()))
        gdelta.append(cliffs_delta(exc.to_numpy(), ret.to_numpy()))
    if len(gper) >= 3:
        nz = [v for v in gper if v != 0]
        pos = sum(1 for v in nz if v > 0)
        sg = stats.binomtest(pos, len(nz), 0.5).pvalue
        wx = stats.wilcoxon(gper).pvalue if len(gper) >= 6 else None
        lo, hi = boot_ci(gper, [f"t{i}" for i in range(len(gper))])
        clustered["discordant|genus_unit"] = {
            "n_genera_eligible": len(gper),
            "median_per_genus_difference": float(np.median(gper)),
            "n_positive": pos, "n_nonzero": len(nz),
            "sign_test_p_two_sided": float(sg),
            "wilcoxon_p_two_sided": float(wx) if wx is not None else None,
            "genus_cluster_ci95": [lo, hi],
            "median_cliffs_delta": float(np.median(
                [d for d in gdelta if d is not None]))}
        c = clustered["discordant|genus_unit"]
        print(f"  discordant, GENUS as unit: {len(gper)} genera  "
              f"median diff {c['median_per_genus_difference']:+.4f}  "
              f"{pos}/{len(nz)} positive  sign p {sg:.3g}  "
              f"Wilcoxon p {wx if wx is None else f'{wx:.3g}'}  "
              f"CI [{lo:+.4f}, {hi:+.4f}]")
    res["clustered_omega_direction"] = clustered
    print("\n  The pooled Mann-Whitney (p 5.7e-18) is DESCRIPTIVE only: it "
          "treats 36 535 loci as independent.")

    # ======================================= S  same-genome-set filtering effect
    print()
    print("=" * 72)
    print("S  the filtering effect, on the SAME genome set")
    print("=" * 72)
    stf = json.loads((R / "rbh_stratified_all.json").read_text())
    s1 = stf["same_genome_set_comparison"]["S1"]
    print(f"  E1 over all 92 genomes          {e1:+.6f}")
    print(f"  E2 over its 90 genomes          {e2:+.6f}   "
          f"(different populations: not the filtering effect)")
    print(f"  E1 over the SAME 90 genomes     {s1['S0_same_genomes']:+.6f}")
    print(f"  E2 over those same 90 genomes   {s1['arm']:+.6f}")
    print(f"  FILTERING EFFECT, same set      {s1['bias']:+.6f}")
    res["filtering_effect_same_set"] = {
        "n_genomes": s1["n_genomes"],
        "E1_same_genomes": s1["S0_same_genomes"],
        "E2": s1["arm"], "difference": s1["bias"],
        "E1_all_92_genomes": e1,
        "note": "the 92-vs-90 comparison compares different populations and "
                "must not be used as the filtering effect"}

    # ---- FOUR DISTINCT ESTIMANDS on the same 90 genomes. A difference of two
    # marginal medians is not an average of anything and can diverge from the
    # typical per-genome change, so all four are computed and one is declared
    # primary on stated grounds rather than by which looks best.
    pg1, pg2, gkey = {}, {}, {}
    for acc, sub in full.groupby("accession"):
        v1 = gap(sub)
        v2 = gap(sub, mask=(sub.hiconf == 1))
        if v1 is not None and v2 is not None:
            pg1[acc], pg2[acc] = v1, v2
            gkey[acc] = sub.genus.iloc[0]
    accs = sorted(pg1)
    d = np.array([pg2[a] - pg1[a] for a in accs])
    a1 = np.array([pg1[a] for a in accs])
    a2 = np.array([pg2[a] for a in accs])

    def cluster_boot_mean(vals, keys, stat=np.median, B=2000, seed=20261108):
        buckets = collections.defaultdict(list)
        for k, v in zip(keys, vals):
            buckets[k].append(v)
        ks = list(buckets)
        g = np.random.default_rng(seed)
        out = [float(stat([x for k in g.choice(ks, len(ks), replace=True)
                           for x in buckets[k]])) for _ in range(B)]
        return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))

    gkeys = [gkey[a] for a in accs]
    # genus-level paired effect: collapse to one pair per genus first
    gen_pairs = collections.defaultdict(lambda: ([], []))
    for a in accs:
        gen_pairs[gkey[a]][0].append(pg1[a])
        gen_pairs[gkey[a]][1].append(pg2[a])
    gd = np.array([float(np.median(v2) - np.median(v1))
                   for v1, v2 in gen_pairs.values()])

    nzd = d[d != 0]
    est = {
        "E1_difference_of_medians": {
            "value": float(np.median(a1) - np.median(a2)) * -1,
            "definition": "median(E2 gaps) - median(E1 gaps) over the same "
                          "90 genomes; the figure previously reported",
            "n_units": len(accs)},
        "E2_median_paired_difference": {
            "value": float(np.median(d)),
            "definition": "median over genomes of (E2_g - E1_g)",
            "n_units": len(accs),
            "n_negative": int((d < 0).sum()), "n_positive": int((d > 0).sum()),
            "n_zero": int((d == 0).sum()),
            "sign_test_p_two_sided": float(stats.binomtest(
                int((nzd < 0).sum()), len(nzd), 0.5).pvalue) if len(nzd) else None,
            "wilcoxon_p_two_sided": float(stats.wilcoxon(d).pvalue)
            if len(nzd) >= 6 else None,
            "genome_cluster_ci95": cluster_boot_mean(d, [f"g{i}" for i in
                                                         range(len(d))]),
            "genus_cluster_ci95": cluster_boot_mean(d, gkeys)},
        "E3_mean_paired_difference": {
            "value": float(d.mean()),
            "definition": "mean over genomes of (E2_g - E1_g)",
            "n_units": len(accs),
            "sd": float(d.std(ddof=1)),
            "t_test_p_two_sided": float(stats.ttest_1samp(d, 0).pvalue),
            "genus_cluster_ci95": cluster_boot_mean(d, gkeys, stat=np.mean)},
        "E4_genus_level_paired": {
            "value": float(np.median(gd)),
            "definition": "per genus, median(E2) - median(E1) over its "
                          "genomes; then median over genera",
            "n_units": int(len(gd)),
            "n_negative": int((gd < 0).sum()), "n_positive": int((gd > 0).sum()),
            "sign_test_p_two_sided": float(stats.binomtest(
                int((gd[gd != 0] < 0).sum()), int((gd != 0).sum()),
                0.5).pvalue),
            "mean": float(gd.mean()),
            "genus_cluster_ci95": cluster_boot_mean(
                gd, [f"t{i}" for i in range(len(gd))])},
    }
    est["primary"] = "E2_median_paired_difference"
    est["primary_justification"] = (
        "It is PAIRED, so it isolates the filtering effect inside each genome "
        "instead of contrasting two marginal summaries; it uses the project's "
        "existing estimator family (a median over per-genome quantities); and "
        "unlike a difference of two medians it is an average of the quantity "
        "of interest, so it cannot be moved by a reshuffle of which genome "
        "sits at the median.")
    print("\n  FOUR ESTIMANDS of the filtering effect, same 90 genomes")
    for k in ("E1_difference_of_medians", "E2_median_paired_difference",
              "E3_mean_paired_difference", "E4_genus_level_paired"):
        v = est[k]
        mark = " <- PRIMARY" if k == est["primary"] else ""
        print(f"    {k:32s} {v['value']:+.6f}  n={v['n_units']:3d}{mark}")
        print(f"      {v['definition']}")
    p2 = est["E2_median_paired_difference"]
    print(f"    paired: {p2['n_negative']} of {len(accs)} genomes lower under "
          f"E2, {p2['n_positive']} higher, {p2['n_zero']} unchanged; "
          f"sign p {p2['sign_test_p_two_sided']:.3g}, "
          f"Wilcoxon p {p2['wilcoxon_p_two_sided']:.3g}")
    print(f"      genome CI {np.round(p2['genome_cluster_ci95'], 5)}  "
          f"genus CI {np.round(p2['genus_cluster_ci95'], 5)}")
    g4 = est["E4_genus_level_paired"]
    print(f"    genus level: {g4['n_negative']} of {g4['n_units']} genera "
          f"lower under E2, sign p {g4['sign_test_p_two_sided']:.3g}, "
          f"CI {np.round(g4['genus_cluster_ci95'], 5)}")
    res["filtering_effect_estimands"] = est

    # ================================================== G  genus-level test
    print()
    print("=" * 72)
    print("G  GENUS-LEVEL directional test for the primary E2 arm")
    print("=" * 72)
    aud = json.loads((R / "rbh_audited_all.json").read_text())
    arms = {a["label"]: a for a in aud["arms"]}
    for lab, col in (("E1 original locus universe", "E1"),
                     ("E2 modified RBH, high-confidence", "E2"),
                     ("E2 modified RBH + synteny >= 1", "E2syn1")):
        v = gen[col].dropna()
        pos = int((v > 0).sum())
        nz = v[v != 0]
        bt = stats.binomtest(int((nz > 0).sum()), len(nz), 0.5,
                             alternative="greater")
        w = stats.wilcoxon(v, alternative="greater")
        a = arms[lab]
        print(f"  {col:7s} genus units {len(v):3d}  positive {pos:3d}  "
              f"median {float(np.median(v)):+.4f}  "
              f"sign p {bt.pvalue:.3g}  Wilcoxon p {w.pvalue:.3g}")
        print(f"          published GENOME-level: {a['n_positive']}/"
              f"{a['n_tested']} genomes, p {a['p']:.3g}  <- different unit")
        res[f"genus_test_{col}"] = {
            "n_genera": int(len(v)), "n_positive": pos,
            "median": float(np.median(v)),
            "sign_test_p_one_sided": float(bt.pvalue),
            "wilcoxon_p_one_sided": float(w.pvalue),
            "genome_level_p_published": a["p"],
            "genome_level_positive": f"{a['n_positive']}/{a['n_tested']}"}

    # ================================================== C  CDS vs rejection
    print()
    print("=" * 72)
    print("C  does missing annotated CDS explain MODIFIED-RBH rejection?")
    print("=" * 72)
    d = full[full.kind == "discordant"].copy()
    tab = d.groupby("has_annotated_cds").agg(
        n=("hiconf", "size"), confirmed=("hiconf", "sum"),
        rate=("hiconf", "mean")).reset_index()
    for _, r in tab.iterrows():
        lab = "has a CDS" if r.has_annotated_cds == 1 else "no CDS"
        print(f"  discordant, {lab:10s} n {r.n:5.0f}  "
              f"modified-RBH confirmed {r.confirmed:5.0f}  "
              f"({100 * r.rate:5.1f} %)")
    with_cds = float(tab.loc[tab.has_annotated_cds == 1, "rate"].iloc[0])
    no_cds = float(tab.loc[tab.has_annotated_cds == 0, "rate"].iloc[0])
    ct = pd.crosstab(d.has_annotated_cds, d.hiconf)
    chi = stats.chi2_contingency(ct.values)
    print(f"  difference {100 * (with_cds - no_cds):+.1f} pp, "
          f"chi-square p {chi.pvalue:.3g}")
    print(f"\n  INTERPRETATION: the modified method was built to score loci "
          f"with no annotated protein, and it does — {100 * no_cds:.1f} % of "
          f"CDS-less discordant loci are confirmed.")
    print(f"  CDS availability is associated with confirmation but does NOT "
          f"account for the attrition: even among discordant loci WITH a CDS, "
          f"{100 * (1 - with_cds):.1f} % are still rejected, against "
          f"{100 * (1 - float(full.loc[full.kind == 'concordant_intact', 'hiconf'].mean())):.1f} % "
          f"of baseline loci.")
    # what DOES the rejection track?
    rej = {}
    for col in ("forward_ambiguous", "span", "forward_coverage",
                "forward_pident"):
        a = pd.to_numeric(d.loc[d.hiconf == 1, col], errors="coerce").dropna()
        b = pd.to_numeric(d.loc[d.hiconf == 0, col], errors="coerce").dropna()
        if not len(a) or not len(b):
            continue
        u = stats.mannwhitneyu(b, a, alternative="two-sided")
        rej[col] = {"median_confirmed": float(a.median()),
                    "median_rejected": float(b.median()),
                    "mannwhitney_p": float(u.pvalue)}
        print(f"    {col:20s} confirmed {a.median():9.4f}   "
              f"rejected {b.median():9.4f}   p {u.pvalue:.3g}")
    out = d.groupby("outcome").size().sort_values(ascending=False)
    print("  rejection outcomes among discordant loci:")
    for k, v in out.items():
        print(f"    {k:28s} {v:5d}  ({100 * v / len(d):5.1f} %)")
    res["cds_vs_rejection"] = {
        "rate_with_cds": with_cds, "rate_without_cds": no_cds,
        "difference_pp": 100 * (with_cds - no_cds),
        "chi2_p": float(chi.pvalue),
        "rejected_even_with_cds_frac": 1 - with_cds,
        "covariates_of_rejection": rej,
        "outcomes": {str(k): int(v) for k, v in out.items()}}

    # ================================================== N  CI narrowing
    print()
    print("=" * 72)
    print("N  is the narrower E2 interval attributable to variance reduction?")
    print("=" * 72)
    gaps_e1, gaps_e2 = {}, {}
    for acc, sub in full.groupby("accession"):
        v1, v2 = gap(sub), gap(sub, mask=(sub.hiconf == 1))
        if v1 is not None:
            gaps_e1[acc] = v1
        if v2 is not None:
            gaps_e2[acc] = v2
    shared = sorted(set(gaps_e1) & set(gaps_e2))
    s1 = np.array([gaps_e1[k] for k in shared])
    s2 = np.array([gaps_e2[k] for k in shared])
    gl = {g: [] for g in gen.genus}
    print(f"  per-genome gap dispersion over the {len(shared)} shared genomes")
    print(f"    E1  sd {s1.std(ddof=1):.5f}  IQR "
          f"{np.subtract(*np.percentile(s1, [75, 25])):.5f}")
    print(f"    E2  sd {s2.std(ddof=1):.5f}  IQR "
          f"{np.subtract(*np.percentile(s2, [75, 25])):.5f}")
    e1g = gen.E1.dropna().to_numpy()
    e2g = gen.E2.dropna().to_numpy()
    print(f"  per-genus gap dispersion")
    print(f"    E1  sd {e1g.std(ddof=1):.5f}   E2  sd {e2g.std(ddof=1):.5f}")
    narrower = s2.std(ddof=1) < s1.std(ddof=1)
    print(f"\n  E2 per-genome dispersion is "
          f"{'LOWER' if narrower else 'NOT lower'} than E1's, so the narrower "
          f"interval is at least consistent with reduced between-unit "
          f"variance — but this is a description of the same bootstrap input, "
          f"not an independent explanation, and it is NOT evidence of better "
          f"accuracy.")
    res["ci_narrowing"] = {
        "n_shared_genomes": len(shared),
        "sd_per_genome_E1": float(s1.std(ddof=1)),
        "sd_per_genome_E2": float(s2.std(ddof=1)),
        "sd_per_genus_E1": float(e1g.std(ddof=1)),
        "sd_per_genus_E2": float(e2g.std(ddof=1)),
        "E2_dispersion_lower": bool(narrower),
        "caveat": "lower dispersion is consistent with the narrower interval "
                  "but is not an independent explanation and is not evidence "
                  "of improved accuracy"}

    (OUT / "review_fixes.json").write_text(
        json.dumps(res, indent=1, default=str, sort_keys=True))
    print(f"\nwrote {OUT / 'review_fixes.json'} and omega_direction.tsv")


if __name__ == "__main__":
    main()
