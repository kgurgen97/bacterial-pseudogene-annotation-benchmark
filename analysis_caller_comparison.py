#!/usr/bin/env python3
"""The corrected caller picture, compared matched rather than by p-values.

POST HOC. Not registered, not verdicts. Replaces the withdrawn DFAST_RESULTS.md,
whose every claim rested on a pseudogene set that was 100 % product-name matches.

The methodological point: it is wrong to rank callers by comparing the p-values
of separate tests. A caller with more loci gets a smaller p for the same effect.
So the headline comparison here is MATCHED - within each genome that carries
enough of both callers' unique calls, the difference between their gaps is
computed, and a sign test is run on those differences.

Five parts:
  M1  matched per-genome comparison of every caller pair
  M2  each caller's sole-call set against the baseline, with genus-clustered
      bootstrap intervals rather than bare p-values
  M3  what DFAST's corrected calls actually are: frameshift, internal stop, or
      both
  M4  the shape of each sole set - length, omega, omega > 1 fraction
  M5  distribution across reduction strata

    python3 analysis_caller_comparison.py
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
CALLERS = ["PGAP", "Pseudofinder", "Bakta", "DFAST"]
MIN_LOCI = 3


def sign_test(v, alternative="two-sided"):
    nz = [x for x in v if x != 0.0]
    if len(nz) < 5:
        return None, None, len(nz)
    pos = sum(1 for x in nz if x > 0)
    return (float(stats.binomtest(pos, len(nz), 0.5,
                                  alternative=alternative).pvalue),
            pos, len(nz))


def cluster_boot(sub, base, unit, B=2000, seed=20261104):
    """95 % CI on the pooled (sole-call minus baseline) gap, resampling clusters."""
    keys = sorted(set(sub[unit]) | set(base[unit]))
    if len(keys) < 3:
        return None, None
    a_by = {k: s.omega.values for k, s in sub.groupby(unit)}
    b_by = {k: s.omega.values for k, s in base.groupby(unit)}
    g = np.random.default_rng(seed)
    vals = []
    for _ in range(B):
        draw = g.choice(keys, len(keys), replace=True)
        a = np.concatenate([a_by[k] for k in draw if k in a_by]) if any(
            k in a_by for k in draw) else np.array([])
        b = np.concatenate([b_by[k] for k in draw if k in b_by]) if any(
            k in b_by for k in draw) else np.array([])
        if len(a) and len(b):
            vals.append(float(np.median(a) - np.median(b)))
    if not vals:
        return None, None
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def main() -> None:
    t = pd.read_csv(R / "loci.tsv", sep="\t")
    dis = t[t.kind == "discordant"]
    base = t[t.kind == "concordant_intact"]
    out = {"note": "POST HOC; matched comparisons, not p-value ranking"}
    print(f"{len(t)} loci: {len(dis)} discordant, {len(base)} concordant-intact\n")

    sole = {c: dis[dis.pseudo_callers == c] for c in CALLERS}

    # ---------------------------------------------------------------- M1
    print("M1  MATCHED per-genome comparison of caller pairs")
    print("    within genomes carrying enough of BOTH callers' sole calls,")
    print("    the difference of their gaps, then a sign test on it\n")
    m1 = {}
    for i, a in enumerate(CALLERS):
        for b in CALLERS[i + 1:]:
            diffs, accs = [], []
            for acc in sorted(set(sole[a].accession) & set(sole[b].accession)):
                sa = sole[a][sole[a].accession == acc].omega
                sb = sole[b][sole[b].accession == acc].omega
                bb = base[base.accession == acc].omega
                if len(sa) >= MIN_LOCI and len(sb) >= MIN_LOCI and len(bb) >= MIN_LOCI:
                    ga = np.median(sa) - np.median(bb)
                    gb = np.median(sb) - np.median(bb)
                    diffs.append(float(ga - gb))
                    accs.append(acc)
            p, pos, n = sign_test(diffs)
            m1[f"{a} vs {b}"] = {
                "n_genomes": len(diffs), "n_tested": n,
                "median_difference": float(np.median(diffs)) if diffs else None,
                "n_positive": pos, "p": p}
            if p is None:
                print(f"  {a:13s} vs {b:13s} {len(diffs):3d} shared genomes"
                      f"  NOT TESTABLE (needs 5, has {n})")
            else:
                who = a if np.median(diffs) > 0 else b
                print(f"  {a:13s} vs {b:13s} {len(diffs):3d} shared genomes"
                      f"  median diff {np.median(diffs):+.4f}"
                      f"  {pos}/{n}  p {p:.3g}   higher: {who}")
    out["M1_matched_pairs"] = m1

    # ---------------------------------------------------------------- M2
    print("\nM2  each caller's sole calls against the baseline, clustered")
    m2 = {}
    for c in CALLERS:
        s = sole[c]
        gaps = []
        for acc in sorted(set(s.accession)):
            sa = s[s.accession == acc].omega
            bb = base[base.accession == acc].omega
            if len(sa) >= MIN_LOCI and len(bb) >= MIN_LOCI:
                gaps.append(float(np.median(sa) - np.median(bb)))
        p, pos, n = sign_test(gaps)
        glo, ghi = cluster_boot(s, base, "genus")
        m2[c] = {"n_loci": int(len(s)), "n_genomes": len(gaps),
                 "median_gap": float(np.median(gaps)) if gaps else None,
                 "n_positive": pos, "p": p, "genus_cluster_ci95": [glo, ghi]}
        if p is None:
            print(f"  {c:13s} {len(s):5d} loci, {len(gaps):3d} genomes"
                  f"  NOT TESTABLE")
        else:
            ci = (f"[{glo:+.4f}, {ghi:+.4f}] "
                  f"{'excludes 0' if glo and glo > 0 else 'INCLUDES 0'}"
                  if glo is not None else "no CI")
            print(f"  {c:13s} {len(s):5d} loci, {len(gaps):3d} genomes"
                  f"  {np.median(gaps):+.4f}  {pos}/{n}  p {p:.3g}"
                  f"  genus CI {ci}")
    out["M2_sole_vs_baseline"] = m2

    # ---------------------------------------------------------------- M3
    print("\nM3  what DFAST's corrected pseudogene calls are, corpus-wide")
    cls = collections.Counter()
    for d in sorted((D / "n2-dfast").iterdir()):
        f = d / "pseudogene_summary.tsv"
        if not f.exists():
            continue
        for line in f.read_text(errors="replace").splitlines()[1:]:
            r = line.split("\t")
            if len(r) >= 4:
                cls[r[3]] += 1
    tot_dis = sum(v for k, v in cls.items() if k != "partial")
    for k, v in cls.most_common():
        mark = "  (NOT a pseudogene call)" if k == "partial" else ""
        print(f"    {k:34s} {v:7d}{mark}")
    print(f"    disrupted classes used as pseudogenes: {tot_dis}")
    out["M3_dfast_classes"] = dict(cls)

    # ---------------------------------------------------------------- M4
    print("\nM4  shape of each sole-call set")
    rows = []
    for c in CALLERS + ["(baseline)"]:
        s = base if c == "(baseline)" else sole[c]
        if not len(s):
            continue
        rows.append({"set": c, "n": int(len(s)),
                     "median_omega": float(np.median(s.omega)),
                     "median_codons": float(np.median(s.codons)),
                     "frac_omega_gt_1": float((s.omega > 1).mean())})
    print(pd.DataFrame(rows).to_string(index=False,
                                       float_format=lambda v: f"{v:.4f}"))
    out["M4_shape"] = rows

    # ---------------------------------------------------------------- M5
    print("\nM5  sole calls per 1 000 loci in each stratum")
    allloci = t.groupby("stratum").size()
    m5 = {}
    hdr = f"  {'stratum':28s}" + "".join(f"{c:>14s}" for c in CALLERS)
    print(hdr)
    for st in sorted(allloci.index):
        line = f"  {st:28s}"
        m5[st] = {}
        for c in CALLERS:
            n = int((sole[c].stratum == st).sum())
            rate = 1000 * n / allloci[st]
            m5[st][c] = {"n": n, "per_1000": float(rate)}
            line += f"{n:6d} ({rate:4.1f})"
        print(line)
    out["M5_strata"] = m5

    (R / "caller_comparison.json").write_text(
        json.dumps(out, indent=2, sort_keys=True) + "\n")
    print("\nwrote results/caller_comparison.json")


if __name__ == "__main__":
    main()
