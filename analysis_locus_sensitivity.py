#!/usr/bin/env python3
"""How fragile is the corrected G3? A post hoc sensitivity battery.

POST HOC THROUGHOUT. None of this is registered, none of it is a verdict, and
none of it may be reported as confirming LOCUS_PLAN.md. G3's registered verdict
stands or falls on the frozen analysis; this asks how much the number moves
when defensible choices are made differently.

The registered G3: per genome, median omega of discordant loci minus median of
concordant-intact loci, then a sign test over genomes with exact ties dropped.
Corrected value +0.031 at p = 6.2e-12 over 92 genomes.

Seven checks, each a choice someone could reasonably have made otherwise:

  S1  genus-clustered bootstrap on the pooled gap, since genomes within a genus
      are not independent
  S2  leave-one-genus-out, to find a genus carrying the result
  S3  a stricter baseline: concordant-intact requiring three callers to agree
      the gene is intact rather than two
  S4  a higher minimum of loci per genome than the registered three
  S5  caller dropout: remove the loci where one caller is the SOLE pseudogene
      caller, for each caller in turn, to see whether one tool drives G3
  S6  stronger discordance: restrict to loci where at least two callers call a
      pseudogene
  S7  stratum dropout, to check the reduced genomes are not carrying it

Not done here: varying the 50 % reciprocal-overlap threshold, which requires
rebuilding the locus set and is a separate job. Recorded as a gap.

    python3 analysis_locus_sensitivity.py
"""
import collections
import json
import pathlib

import numpy as np
import pandas as pd
from scipy import stats

P = pathlib.Path(__file__).resolve().parent
R = P / "results"
MIN_LOCI = 3            # the registered minimum per class per genome


def g3(t, min_loci=MIN_LOCI):
    """The registered G3 statistic on whatever subset is handed in."""
    gaps = {}
    for acc, g in t.groupby("accession"):
        d = g.loc[g.kind == "discordant", "omega"]
        c = g.loc[g.kind == "concordant_intact", "omega"]
        if len(d) >= min_loci and len(c) >= min_loci:
            gaps[acc] = float(np.median(d) - np.median(c))
    nz = [v for v in gaps.values() if v != 0.0]
    if len(nz) < 5:
        return {"n_genomes": len(gaps), "p": None, "median_gap": None,
                "note": f"only {len(nz)} non-tied genomes"}
    pos = sum(1 for v in nz if v > 0)
    return {"n_genomes": len(gaps), "n_tested": len(nz),
            "median_gap": float(np.median(list(gaps.values()))),
            "n_positive": pos,
            "p": float(stats.binomtest(pos, len(nz), 0.5).pvalue)}


def line(name, r, base=None):
    if r.get("p") is None:
        return f"  {name:44s} {r['n_genomes']:3d} genomes  NOT TESTABLE ({r['note']})"
    s = (f"  {name:44s} {r['n_genomes']:3d} genomes  {r['median_gap']:+.4f}  "
         f"{r['n_positive']}/{r['n_tested']}  p {r['p']:.3g}")
    if base and base.get("median_gap"):
        s += f"   ({100 * r['median_gap'] / base['median_gap']:.0f}% of base)"
    return s


def cluster_boot(t, unit, B=2000, seed=20261103):
    """Resample whole clusters, recomputing the POOLED gap each time."""
    keys = t[unit].unique()
    byk = {k: s for k, s in t.groupby(unit)}
    if len(keys) < 3:
        return None, None
    g = np.random.default_rng(seed)
    vals = []
    for _ in range(B):
        c = pd.concat([byk[k] for k in g.choice(keys, len(keys), replace=True)])
        d = c.loc[c.kind == "discordant", "omega"]
        i = c.loc[c.kind == "concordant_intact", "omega"]
        if len(d) and len(i):
            vals.append(float(np.median(d) - np.median(i)))
    if not vals:
        return None, None
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def callers_of(s):
    return set(str(s).split("+")) if isinstance(s, str) and s else set()


def main() -> None:
    t = pd.read_csv(R / "loci.tsv", sep="\t")
    t["pseudo_set"] = t.pseudo_callers.apply(callers_of)
    t["intact_set"] = t.intact_callers.apply(callers_of)
    base = g3(t)
    out = {"note": "POST HOC sensitivity battery; not registered, not verdicts",
           "base": base, "checks": {}}
    print(f"{len(t)} loci, {t.accession.nunique()} genomes, "
          f"{t.genus.nunique()} genera")
    print(line("registered G3, as recomputed", base))

    # ---- S1 cluster bootstraps on the pooled gap
    print("\nS1  cluster bootstraps on the POOLED gap (not the per-genome median)")
    for unit in ("accession", "genus"):
        lo, hi = cluster_boot(t, unit)
        out["checks"][f"S1_{unit}_cluster_ci95"] = [lo, hi]
        if lo is None:
            print(f"  {unit:12s} too few clusters")
        else:
            excl = "excludes 0" if lo > 0 else "INCLUDES 0"
            print(f"  {unit:12s} 95% CI [{lo:+.4f}, {hi:+.4f}]  {excl}")

    # ---- S2 leave one genus out
    print("\nS2  leave-one-genus-out (worst eight shown)")
    rows = []
    for gen in sorted(t.genus.dropna().unique()):
        r = g3(t[t.genus != gen])
        if r.get("p") is not None:
            rows.append((gen, r))
    rows.sort(key=lambda x: x[1]["median_gap"])
    out["checks"]["S2_leave_one_genus_out"] = {
        g: {k: v for k, v in r.items()} for g, r in rows}
    for gen, r in rows[:8]:
        print(line(f"without {gen}", r, base))
    worst = rows[0][1] if rows else None
    allsig = all(r["p"] < 0.05 and r["median_gap"] > 0 for _, r in rows)
    print(f"  -> every leave-one-out remains positive and significant: {allsig}")

    # ---- S3 a stricter baseline
    print("\nS3  baseline requiring three callers to agree the gene is intact")
    strict = t[(t.kind == "discordant") |
               ((t.kind == "concordant_intact") &
                (t.intact_set.apply(len) >= 3))]
    r = g3(strict)
    out["checks"]["S3_baseline_3_callers"] = r
    n_base = int((t.kind == "concordant_intact").sum())
    n_strict = int(((strict.kind == "concordant_intact")).sum())
    print(f"  baseline loci {n_base} -> {n_strict}")
    print(line("3-caller baseline", r, base))

    # ---- S4 minimum loci per genome
    print("\nS4  minimum loci per class per genome")
    for m in (3, 5, 10, 20, 50):
        r = g3(t, min_loci=m)
        out["checks"][f"S4_min_{m}"] = r
        print(line(f"min {m} per class", r, base))

    # ---- S5 caller dropout
    print("\nS5  dropping the loci where one caller is the SOLE pseudogene caller")
    for c in ("PGAP", "Pseudofinder", "Bakta", "DFAST"):
        keep = t[~((t.kind == "discordant") & (t.pseudo_callers == c))]
        n_drop = len(t) - len(keep)
        r = g3(keep)
        out["checks"][f"S5_without_sole_{c}"] = r
        print(line(f"without {c}-only loci (-{n_drop})", r, base))

    # ---- S6 stronger discordance
    print("\nS6  discordant loci where at least two callers say pseudogene")
    strong = t[(t.kind == "concordant_intact") |
               ((t.kind == "discordant") & (t.pseudo_set.apply(len) >= 2))]
    n_d = int((strong.kind == "discordant").sum())
    r = g3(strong)
    out["checks"]["S6_two_caller_discordance"] = r
    print(f"  discordant loci {int((t.kind == 'discordant').sum())} -> {n_d}")
    print(line("two-caller discordance", r, base))

    # ---- S7 stratum dropout
    print("\nS7  dropping each reduction stratum in turn")
    for s in sorted(t.stratum.dropna().unique()):
        r = g3(t[t.stratum != s])
        out["checks"][f"S7_without_{s}"] = r
        print(line(f"without {s}", r, base))

    # ---- summary
    testable = [(k, v) for k, v in out["checks"].items()
                if isinstance(v, dict) and v.get("p") is not None]
    pos_sig = [k for k, v in testable if v["p"] < 0.05 and v["median_gap"] > 0]
    print("\n" + "=" * 70)
    print(f"checks yielding a testable G3: {len(testable)}")
    print(f"  positive and significant at 0.05: {len(pos_sig)}")
    failed = [k for k, v in testable if not (v["p"] < 0.05 and v["median_gap"] > 0)]
    if failed:
        print("  NOT positive-and-significant under:")
        for k in failed:
            v = dict(testable)[k]
            print(f"    {k}: {v['median_gap']:+.4f} p {v['p']:.3g}")
    else:
        print("  G3 stays positive and significant under every check in this battery")
    out["summary"] = {"n_testable": len(testable),
                      "n_positive_significant": len(pos_sig),
                      "failed_checks": failed,
                      "leave_one_genus_out_all_significant": allsig,
                      "not_done": ["reciprocal overlap threshold variation: "
                                   "needs the locus set rebuilt"]}
    (R / "locus_sensitivity.json").write_text(
        json.dumps(out, indent=2, sort_keys=True, default=str) + "\n")
    print("\nwrote results/locus_sensitivity.json")


if __name__ == "__main__":
    main()
