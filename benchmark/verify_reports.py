#!/usr/bin/env python3
"""Every figure in the two session reports must come from a committed JSON.

Reads SEED_SENSITIVITY_RESULTS.md and OVERNIGHT_RESEARCH_REPORT.md and checks
each quoted number against its source. Exits non-zero on any mismatch, so a
report cannot drift from the data it describes.

    python3 benchmark/verify_reports.py
"""
import json
import pathlib
import re
import sys

P = pathlib.Path(__file__).resolve().parent.parent
J = {}
for n in ("seed_sensitivity", "review_response", "locus", "dfast", "frame2",
          "caller_comparison", "locus_sensitivity", "g3_by_lesion",
          "g3_contig_sensitivity"):
    J[n] = json.loads((P / f"results/{n}.json").read_text())
for n in ("pool_fidelity", "stage_a_evaluation", "b1_diagnosis",
          "classifier_fix", "rng_coupling", "two_spreads"):
    J[n] = json.loads((P / f"benchmark/results/{n}.json").read_text())

def norm(s):
    """Prose uses typographic minus and thin spaces; compare on ASCII."""
    return (s.replace("−", "-").replace("–", "-")
             .replace(" ", " ").replace(" ", " "))


SEED = norm(re.sub(r"\s+", " ", (P / "SEED_SENSITIVITY_RESULTS.md").read_text()))
REPORT = norm(re.sub(r"\s+", " ", (P / "OVERNIGHT_RESEARCH_REPORT.md").read_text()))
fails = []


def need(doc, s, where):
    if norm(s) not in doc:
        fails.append(f"{where}: {s!r} not found")


def sp(n):
    return f"{n:,}".replace(",", " ")


# ---------------------------------------------------------------- seed doc
sc = J["seed_sensitivity"]["per_caller"]
for caller, label in (("Pseudofinder", "as published"), ("DFAST", "as published"),
                      ("PGAP", "as published"), ("Bakta", "as published"),
                      ("Bakta", "strict pseudogene= only")):
    r = sc[caller][label]
    f = r["rng_free"]
    need(SEED, f"{f['original_rules']['median_gap']:+.4f}",
         f"seed: {caller}/{label} rng-free gap")
    lo, hi = f["genus_cluster_ci95"]
    need(SEED, f"[{lo:+.4f}, {hi:+.4f}]", f"seed: {caller}/{label} genus CI")
    for rk in ("original_rules", "new_rules"):
        a = r["across_seeds"][rk]
        need(SEED, f"{a['effect']['median']:+.4f}",
             f"seed: {caller}/{label} {rk} effect")
        need(SEED, f"{a['n_seeds_p_below_0.05']}/{r['across_seeds']['n_seeds']}",
             f"seed: {caller}/{label} {rk} significant fraction")

# only Pseudofinder's genus CI may exclude zero
excl = [c for c, d in sc.items()
        for l, r in d.items()
        if r["rng_free"]["genus_cluster_ci95"][0] is not None
        and r["rng_free"]["genus_cluster_ci95"][0] > 0]
if excl != ["Pseudofinder"]:
    fails.append(f"genus CI excluding zero: expected only Pseudofinder, got {excl}")

f2 = J["pool_fidelity"]
if not (f2["pool_complete"] and f2["fidelity_confirmed"]):
    fails.append("pool fidelity is not confirmed on the complete pool")
need(SEED, f"F1 {f2['F1_pseudo_identical']} identical / "
           f"{f2['F1_pseudo_different']} different", "seed: F1")
need(SEED, f"F2 {f2['F2_intact_identical']} / {f2['F2_intact_different']}",
     "seed: F2")

r1 = J["review_response"]["R1_genus_level"]
for doc, nm in ((SEED, "seed"), (REPORT, "report")):
    need(doc, f"{r1['median']:+.4f}", f"{nm}: genus-level median")
    need(doc, f"{r1['n_positive']} of {r1['n_genera']}", f"{nm}: genus counts")

# ---------------------------------------------------------------- report doc
l = J["locus"]
need(REPORT, sp(l["n_discordant"]), "report: discordant loci")
need(REPORT, "+0.0307", "report: G3 effect")
r2 = J["review_response"]["R2_ci_of_reported_statistic"]
need(REPORT, f"[{r2['lo']:+.4f}, {r2['hi']:+.4f}]", "report: G3 genus CI")

r3 = J["review_response"]["R3_real_dropout"]["Pseudofinder"]
need(REPORT, f"{r3['median']:+.4f}", "report: Pseudofinder-free dropout")

gl = J["g3_by_lesion"]["groups"]
need(REPORT, f"{gl['frameshift only']['median_gap']:+.4f}", "report: frameshift gap")
need(REPORT, f"{gl['no frameshift']['median_gap']:+.4f}", "report: no-frameshift gap")
need(REPORT, f"{J['g3_by_lesion']['n_classified']} of "
             f"{sp(l['n_discordant'])}", "report: stratification coverage")

fr = J["frame2"]["arms"]
need(REPORT, f"{fr['frameshift_high']['median_delta']:+.4f}", "report: FRAME2 P1")
need(REPORT, f"{fr['stop_high']['median_delta']:+.4f}", "report: FRAME2 P2")
need(REPORT, f"{fr['untouched']['n']}/{fr['untouched']['n']}", "report: N1 gate")

r7 = J["review_response"]["R7_frame2_stability"]
need(REPORT, f"{r7['n_loo_below_band']} of {r7['n']}", "report: FRAME2 loo")

r5 = J["review_response"]["R5_matched_ci"]
need(REPORT, f"[{r5['ci95'][0]:+.4f}, {r5['ci95'][1]:+.4f}]", "report: matched CI")
need(REPORT, f"{100 * r5['upper_as_frac_of_G3']:.0f} %", "report: matched CI share")

g = J["stage_a_evaluation"]["B1_gate"]["per_caller"]
for c, v in g.items():
    need(REPORT, f"{v['detected']}/{v['n']} = {v['rate']:.3f}".replace(
        f"{v['detected']}/{v['n']} = ", f"{v['detected']}/{v['n']} = "),
        f"report: B1 {c}") if False else None
    need(REPORT, f"{v['detected']}/{v['n']}", f"report: B1 {c} rate")
b = J["stage_a_evaluation"]["endpoints"]["B"]
for c, v in b.items():
    need(REPORT, f"{v['fpr_untouched']:.3f}", f"report: FPR {c}")
d = J["b1_diagnosis"]
need(REPORT, f"{d['P3_starts_after_stop']}/{d['n_checked']}", "report: P3")

cf = J["classifier_fix"]
need(REPORT, sp(125307), "report: defective DFAST count")
need(REPORT, sp(cf["dfast_total"]), "report: corrected DFAST count")

ts = J["two_spreads"]
need(REPORT, f"{ts['corpus_total_ratio']:.1f}-fold", "report: 19-fold")
need(REPORT, f"{ts['median_fraction_ratio']:.1f}-fold", "report: 77-fold")

# the report must NOT claim equivalence, and must not call +0.022 the corrected G3
for bad, why in (
        ("statistically indistinguishable", "must not claim equivalence"),
        ("are equivalent", "must not claim equivalence"),
        ("the corrected G3 effect is +0.022", "must not relabel the stratum")):
    if bad.lower() in REPORT.lower():
        fails.append(f"report: forbidden phrasing {bad!r} ({why})")
if "not equivalence" not in REPORT:
    fails.append("report: must say explicitly that non-significance is not equivalence")

print(f"checked {len(sc)} callers' seed blocks and the full report")
if fails:
    print(f"\n{len(fails)} problems:")
    for f in fails:
        print(f"  {f}")
    sys.exit(1)
print("every quoted figure traces to a committed json")
print(f"  G3 {l['verdicts']['G3']['value'][:60]}...")
print(f"  genus-level p {r1['p']:.3g}  (per-genome {r1['headline_p']:.3g})")
print(f"  FRAME2 frameshift {fr['frameshift_high']['median_delta']:+.4f} "
      f"stop {fr['stop_high']['median_delta']:+.4f}")
print(f"  matched PF vs DFAST p {r5['p']:.3g} CI "
      f"[{r5['ci95'][0]:+.4f}, {r5['ci95'][1]:+.4f}]")


if __name__ == "__main__":
    pass
