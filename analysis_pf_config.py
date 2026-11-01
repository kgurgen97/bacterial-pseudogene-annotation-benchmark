#!/usr/bin/env python3
"""Is it the caller, or the configuration? As PF_CONFIG_PLAN.md specifies.

Config B differs from A in one documented parameter, --length_pseudo 0.50 against
the default 0.65. If that single change moves Pseudofinder as far as switching to
a different caller does, then every result in this project that is attributed to
caller choice is really about configuration, and the wording has to change.
"""
import hashlib
import json
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import stats

P = pathlib.Path.home() / "research/projects/n2-pseudogene-caller-bias"
R = P / "results"

PLAN = P / "PF_CONFIG_PLAN.md"
sha = hashlib.sha256(PLAN.read_bytes()).hexdigest()
print(f"plan sha256: {sha}")
print("  (quote this in PF_CONFIG_RESULTS.md; the project records plan hashes in "
      "its write-ups)")

b = pd.read_csv(R / "pf_configB_counts.tsv", sep="\t")
print(f"\nconfig B genomes: {len(b)}")
if len(b) < 523:
    print(f"REFUSING: config B has {len(b)} of 523 genomes. Wait for the run.")
    sys.exit(1)

d = pd.read_csv(R / "corpus523.tsv", sep="\t").merge(b, on="accession", how="left")
assert d.pfB_pseudo.notna().all(), "a genome in the corpus has no config B result"

d["frac_PGAP"] = d.pgap_pseudo / (d.pgap_pseudo + d.pgap_cds)
d["frac_PF_A"] = d.pf_pseudo / (d.pf_pseudo + d.pf_intact)
d["frac_PF_B"] = d.pfB_pseudo / (d.pfB_pseudo + d.pfB_intact)
d["frac_DFAST"] = d.dfast_disrupted / (d.dfast_disrupted + d.dfast_cds)
d["frac_Bakta"] = d.bakta_pseudo / (d.bakta_pseudo + d.bakta_cds)

CALLERS = ["frac_PGAP", "frac_PF_A", "frac_DFAST", "frac_Bakta"]
med = {c: float(np.nanmedian(d[c])) for c in CALLERS + ["frac_PF_B"]}
print("\nmedian pseudogene fraction over all 523 genomes")
for k, v in med.items():
    print(f"  {k:12s} {v:.4f}")

between = max(med[c] for c in CALLERS) / max(min(med[c] for c in CALLERS), 1e-12)
within = max(med['frac_PF_A'], med['frac_PF_B']) / max(min(med['frac_PF_A'], med['frac_PF_B']), 1e-12)
print(f"\nbetween-caller fold (4 callers): {between:.1f}x")
print(f"within-caller fold  (A vs B):    {within:.1f}x")

S = ["S4_free_living_breadth", "S1_extreme_obligate",
     "S2_obligate_intracellular", "S3_recent_restriction"]
exp = d[d.stratum.isin(S)]
print(f"\nstrata medians on the {len(exp)} additions")
print(f"  {'stratum':28s} {'PF_A':>9s} {'PF_B':>9s}")
for s in S:
    g = exp[exp.stratum == s]
    print(f"  {s:28s} {np.nanmedian(g.frac_PF_A):9.4f} {np.nanmedian(g.frac_PF_B):9.4f}")


def mw(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    if len(x) < 3 or len(y) < 3:
        return np.nan, np.nan
    u, p = stats.mannwhitneyu(x, y, alternative="two-sided")
    rbc = 2 * u / (len(x) * len(y)) - 1
    return float(p), float(rbc)


g = {s: exp[exp.stratum == s] for s in S}
p31b, r31b = mw(g["S3_recent_restriction"].frac_PF_B, g["S1_extreme_obligate"].frac_PF_B)
p34b, r34b = mw(g["S3_recent_restriction"].frac_PF_B, g["S4_free_living_breadth"].frac_PF_B)
p14b, r14b = mw(g["S1_extreme_obligate"].frac_PF_B, g["S4_free_living_breadth"].frac_PF_B)
p14a, r14a = mw(g["S1_extreme_obligate"].frac_PF_A, g["S4_free_living_breadth"].frac_PF_A)
print(f"\nconfig B: S3 vs S1 p {p31b:.4g} rbc {r31b:+.3f}; "
      f"S3 vs S4 p {p34b:.4g} rbc {r34b:+.3f}")
print(f"S1 vs S4: config A p {p14a:.4g} rbc {r14a:+.3f};  "
      f"config B p {p14b:.4g} rbc {r14b:+.3f}")

# C1 is tested on the original 363. corpus523.tsv carries size_bp only for the
# 160 additions - it came from the expansion manifest - so taking it from there
# gives an all-NaN column and a NaN correlation. The sizes for the original
# corpus live in literature_retest.tsv, which is the table C1 was always computed
# on, so take them from there.
lit = pd.read_csv(R / "literature_retest.tsv", sep="\t")[["accession", "size_bp"]]
orig = d[d.stratum == "original_corpus"].drop(columns=["size_bp"]).merge(
    lit, on="accession", how="left")
sz = pd.to_numeric(orig.size_bp, errors="coerce")
n_sz = int(np.isfinite(sz).sum())
print(f"\nC1 corpus: {len(orig)} genomes, {n_sz} with a size")
assert n_sz >= 300, f"only {n_sz} genomes have a size; C1 cannot be tested"
ok = np.isfinite(sz) & np.isfinite(orig.frac_PF_B)
rho_b, pr_b = stats.spearmanr(sz[ok], orig.frac_PF_B[ok])
ok_a = np.isfinite(sz) & np.isfinite(orig.frac_PF_A)
rho_a, pr_a = stats.spearmanr(sz[ok_a], orig.frac_PF_A[ok_a])
print(f"\nC1 on the original corpus: A rho {rho_a:+.4f} p {pr_a:.3g}; "
      f"B rho {rho_b:+.4f} p {pr_b:.3g}")

verdicts = {}


def verdict(k, text, confirmed, failed, value):
    v = "CONFIRMED" if confirmed else ("FAILED" if failed else "INCONCLUSIVE")
    verdicts[k] = dict(prediction=text, value=value, verdict=v)
    print(f"  {k}  {v:13s} {text}\n        {value}")


print("\n=== pre-committed predictions")
verdict("F1", "the more conservative setting calls fewer",
        med["frac_PF_B"] < med["frac_PF_A"], med["frac_PF_B"] >= med["frac_PF_A"],
        f"median fraction A {med['frac_PF_A']:.4f} vs B {med['frac_PF_B']:.4f}")
verdict("F2", "configuration moves it less than caller does",
        within <= between / 2, within >= between,
        f"within-caller {within:.1f}x against between-caller {between:.1f}x")
verdict("F3", "C3's ordering survives the change",
        (p31b < 0.05 and r31b > 0) and (p34b < 0.05 and r34b > 0),
        (p31b >= 0.05) or (p34b >= 0.05),
        f"S3>S1 p {p31b:.4g} rbc {r31b:+.3f}; S3>S4 p {p34b:.4g} rbc {r34b:+.3f}")
verdict("F4", "the S1-versus-S4 sign survives",
        r14b < 0, r14b > 0,
        f"rbc A {r14a:+.3f} (p {p14a:.4g}) -> B {r14b:+.3f} (p {p14b:.4g})")
verdict("F5", "C1 keeps its direction",
        rho_b > 0, rho_b < 0,
        f"Spearman A {rho_a:+.4f} -> B {rho_b:+.4f}")

d[["accession", "stratum", "frac_PGAP", "frac_PF_A", "frac_PF_B",
   "frac_DFAST", "frac_Bakta"]].to_csv(R / "pf_config.tsv", sep="\t", index=False)
(R / "pf_config.json").write_text(json.dumps(dict(
    plan_sha256=sha, n_genomes=int(len(d)), medians=med,
    between_caller_fold=between, within_caller_fold=within,
    S3_vs_S1_B=dict(p=p31b, rbc=r31b), S3_vs_S4_B=dict(p=p34b, rbc=r34b),
    S1_vs_S4_A=dict(p=p14a, rbc=r14a), S1_vs_S4_B=dict(p=p14b, rbc=r14b),
    C1_A=dict(rho=rho_a, p=pr_a), C1_B=dict(rho=rho_b, p=pr_b),
    verdicts=verdicts), indent=1))
nf = sum(1 for v in verdicts.values() if v["verdict"] == "FAILED")
print(f"\nwrote results/pf_config.{{tsv,json}} - {nf} of {len(verdicts)} FAILED")
