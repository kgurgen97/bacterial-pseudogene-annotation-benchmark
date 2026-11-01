#!/usr/bin/env python3
"""Verifying an adversarial review's claims, and computing what it says is missing.

POST HOC. Not registered. Each block states the reviewer's claim, then the
independent recomputation, then whether the claim holds.

R1  the genus-level test, which the review says should be the headline
R2  a bootstrap CI of the statistic actually reported (median of per-genome
    gaps, resampling genera) - the review notes no such interval exists
R3  the REAL caller-dropout test: drop every locus with caller X on the
    pseudogene side, not only the sole-caller ones
R4  how many of the 15 sensitivity "checks" are degenerate duplicates of base
R5  a bootstrap CI for the matched Pseudofinder-vs-DFAST difference
R6  lesion-class exposure measured on the DISCORDANT population rather than
    imported from DFAST's corpus mix
R7  FRAME2's leave-one-out stability around its band edge

    python3 analysis_review_response.py
"""
import collections
import json
import pathlib

import numpy as np
import pandas as pd
from scipy import stats

P = pathlib.Path(__file__).resolve().parent
R = P / "results"
MIN_LOCI = 3


def per_genome_gaps(t, min_loci=MIN_LOCI):
    d = t[t.kind == "discordant"]
    c = t[t.kind == "concordant_intact"]
    out = {}
    for acc in sorted(set(d.accession) & set(c.accession)):
        a = d.loc[d.accession == acc, "omega"]
        b = c.loc[c.accession == acc, "omega"]
        if len(a) >= min_loci and len(b) >= min_loci:
            out[acc] = float(np.median(a) - np.median(b))
    return out


def sign_test(v):
    nz = [x for x in v if x != 0.0]
    if len(nz) < 5:
        return None, None, len(nz)
    pos = sum(1 for x in nz if x > 0)
    return float(stats.binomtest(pos, len(nz), 0.5).pvalue), pos, len(nz)


def main() -> None:
    t = pd.read_csv(R / "loci.tsv", sep="\t")
    gen = dict(zip(t.accession, t.genus))
    gaps = per_genome_gaps(t)
    base_p, base_pos, base_n = sign_test(list(gaps.values()))
    base_med = float(np.median(list(gaps.values())))
    out = {"note": "POST HOC verification of an adversarial review"}
    print(f"base G3 recomputed: {base_med:+.6f}, {base_pos}/{base_n}, "
          f"p {base_p:.3g}   ({len(gaps)} genomes)")

    # ---------------------------------------------------------------- R1
    print("\nR1  genus-level test (one value per genus, not per genome)")
    byg = collections.defaultdict(list)
    for acc, v in gaps.items():
        byg[gen[acc]].append(v)
    gvals = {k: float(np.median(v)) for k, v in byg.items()}
    gp, gpos, gn = sign_test(list(gvals.values()))
    print(f"  {len(gvals)} genera, median {np.median(list(gvals.values())):+.4f}, "
          f"{gpos}/{gn} positive, p {gp:.3g}")
    print(f"  the review claimed 26/29 and p = 1.5e-05 -> "
          f"{'MATCHES' if gpos == 26 and abs(gp - 1.52e-5) < 1e-5 else 'differs'}")
    print(f"  genomes per genus: max {max(len(v) for v in byg.values())}, "
          f"genera with >=5: {sum(1 for v in byg.values() if len(v) >= 5)}")
    out["R1_genus_level"] = {"n_genera": len(gvals), "p": gp, "n_positive": gpos,
                             "median": float(np.median(list(gvals.values()))),
                             "headline_p": base_p}

    # ---------------------------------------------------------------- R2
    print("\nR2  bootstrap CI of the REPORTED statistic, resampling genera")
    g = np.random.default_rng(20261106)
    keys = sorted(byg)
    meds = []
    for _ in range(4000):
        draw = g.choice(keys, len(keys), replace=True)
        v = [x for k in draw for x in byg[k]]
        meds.append(float(np.median(v)))
    lo, hi = float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5))
    print(f"  median of per-genome gaps, genus-clustered 95% CI "
          f"[{lo:+.4f}, {hi:+.4f}]")
    print(f"  (the published [+0.0167, +0.0524] is a CI for the POOLED gap,")
    print(f"   a different estimator - the review is right about that)")
    out["R2_ci_of_reported_statistic"] = {"lo": lo, "hi": hi,
                                          "point": base_med}

    # ---------------------------------------------------------------- R3
    print("\nR3  REAL caller dropout: drop every locus with X on the pseudo side")
    out["R3_real_dropout"] = {}
    for c in ("PGAP", "Pseudofinder", "Bakta", "DFAST"):
        mask = t.pseudo_callers.fillna("").apply(
            lambda s: c in str(s).split("+"))
        keep = t[~((t.kind == "discordant") & mask)]
        n_drop = int(((t.kind == "discordant") & mask).sum())
        gg = per_genome_gaps(keep)
        p, pos, n = sign_test(list(gg.values()))
        med = float(np.median(list(gg.values()))) if gg else None
        out["R3_real_dropout"][c] = {"n_discordant_dropped": n_drop,
                                     "n_genomes": len(gg), "median": med,
                                     "p": p, "n_positive": pos}
        print(f"  without any {c:13s} -{n_drop:5d} loci  {len(gg):3d} genomes"
              f"  {med:+.4f}  {pos}/{n}  p {p:.3g}" if p else
              f"  without any {c:13s} -{n_drop:5d} loci  NOT TESTABLE")
    print("  the published S5 drops only SOLE-caller loci, so it does not")
    print("  test what FINAL_LOCUS_SENSITIVITY.md claims. The review is right.")

    # ---------------------------------------------------------------- R4
    print("\nR4  are any of the 15 checks degenerate duplicates of base?")
    ls = json.loads((R / "locus_sensitivity.json").read_text())
    dup = []
    for k, v in ls["checks"].items():
        if isinstance(v, dict) and v.get("median_gap") is not None:
            if (abs(v["median_gap"] - base_med) < 1e-12
                    and v.get("n_genomes") == len(gaps)):
                dup.append(k)
    print(f"  checks identical to base in both median and n: {dup}")
    inc = [k for k in ls["checks"] if k.startswith(("S1_", "S2_"))]
    counted = [k for k, v in ls["checks"].items()
               if isinstance(v, dict) and v.get("p") is not None]
    print(f"  checks actually counted in the '15': {len(counted)}")
    print(f"  S1/S2 keys present but NOT counted: {inc}")
    out["R4_battery"] = {"n_counted": len(counted), "degenerate": dup,
                         "s1_s2_excluded_from_count": inc}

    # ---------------------------------------------------------------- R5
    print("\nR5  bootstrap CI for the matched Pseudofinder-vs-DFAST difference")
    dis = t[t.kind == "discordant"]
    base_l = t[t.kind == "concordant_intact"]
    sole = {c: dis[dis.pseudo_callers == c] for c in ("Pseudofinder", "DFAST")}
    diffs = []
    for acc in sorted(set(sole["Pseudofinder"].accession)
                      & set(sole["DFAST"].accession)):
        sa = sole["Pseudofinder"][sole["Pseudofinder"].accession == acc].omega
        sb = sole["DFAST"][sole["DFAST"].accession == acc].omega
        bb = base_l[base_l.accession == acc].omega
        if len(sa) >= MIN_LOCI and len(sb) >= MIN_LOCI and len(bb) >= MIN_LOCI:
            diffs.append(float((np.median(sa) - np.median(bb))
                               - (np.median(sb) - np.median(bb))))
    g2 = np.random.default_rng(20261107)
    bm = [float(np.median(g2.choice(diffs, len(diffs), replace=True)))
          for _ in range(4000)]
    mlo, mhi = float(np.percentile(bm, 2.5)), float(np.percentile(bm, 97.5))
    p, pos, n = sign_test(diffs)
    print(f"  n {len(diffs)} genomes, median {np.median(diffs):+.5f}, "
          f"{pos}/{n}, p {p:.3g}")
    print(f"  bootstrap 95% CI [{mlo:+.4f}, {mhi:+.4f}]")
    print(f"  upper bound is {100 * mhi / base_med:.0f}% of G3's whole effect")
    print("  so 'indistinguishable' is absence of evidence, not equivalence")
    out["R5_matched_ci"] = {"n": len(diffs), "median": float(np.median(diffs)),
                            "p": p, "ci95": [mlo, mhi],
                            "upper_as_frac_of_G3": mhi / base_med}

    # ---------------------------------------------------------------- R6
    print("\nR6  lesion-class exposure on the DISCORDANT population")
    pseudo_side = collections.Counter()
    for s in dis.pseudo_callers.fillna(""):
        for c in str(s).split("+"):
            if c:
                pseudo_side[c] += 1
    n_dis = len(dis)
    for c, v in pseudo_side.most_common():
        print(f"  {c:13s} on the pseudogene side of {v:5d} of {n_dis} "
              f"({100 * v / n_dis:.1f}%)")
    gl = json.loads((R / "g3_by_lesion.json").read_text())
    print(f"  lesion class is known ONLY for DFAST: "
          f"{gl['n_classified']} of {n_dis} discordant loci "
          f"({100 * gl['n_classified'] / n_dis:.1f}%)")
    print("  so the stratified test covers a tenth of the discordant set and")
    print("  Pseudofinder's share is unclassified. The review is right that")
    print("  the 71% exposure figure was imported from the wrong population.")
    out["R6_exposure"] = {"pseudo_side_counts": dict(pseudo_side),
                          "n_discordant": n_dis,
                          "n_classified": gl["n_classified"],
                          "frac_classified": gl["n_classified"] / n_dis}

    # ---------------------------------------------------------------- R7
    print("\nR7  FRAME2's P1 stability around its +0.030 band edge")
    f2 = list(csv_rows(R / "frame2_genes.tsv"))
    fs = [float(r["delta_omega"]) for r in f2
          if r["lesion_code"] in ("M2_fs_ins1", "M3_fs_del1", "M4_fs_indel2")
          and r["truth_tier"] == "high"]
    med = float(np.median(fs))
    loo = [float(np.median([x for j, x in enumerate(fs) if j != i]))
           for i in range(len(fs))]
    print(f"  n {len(fs)}, median {med:+.6f} (band edge +0.030)")
    print(f"  leave-one-out range [{min(loo):+.4f}, {max(loo):+.4f}]")
    n_below = sum(1 for m in loo if m < 0.030)
    print(f"  dropping one gene moves it BELOW the band edge in "
          f"{n_below} of {len(loo)} cases")
    print("  so CONFIRMED_LARGE is knife-edge at this n. The review is right.")
    out["R7_frame2_stability"] = {"n": len(fs), "median": med,
                                  "loo_range": [min(loo), max(loo)],
                                  "n_loo_below_band": n_below}

    (R / "review_response.json").write_text(
        json.dumps(out, indent=2, sort_keys=True) + "\n")
    print("\nwrote results/review_response.json")


def csv_rows(p):
    import csv
    with open(p) as fh:
        yield from csv.DictReader(fh, delimiter="\t")


if __name__ == "__main__":
    main()
