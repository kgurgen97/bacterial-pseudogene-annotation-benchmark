#!/usr/bin/env python3
"""Every claim in the revised draft and the risk review must trace to a json.

Also enforces the two phrasings the brief forbids: a non-significant matched
comparison must never be called equivalence, and the +0.0220 internal-stop
figure must never be presented as the corrected whole-corpus G3 effect.

    python3 benchmark/verify_draft.py
"""
import json
import pathlib
import re
import sys

P = pathlib.Path(__file__).resolve().parent.parent


def norm(s):
    return (s.replace("−", "-").replace("–", "-")
             .replace(" ", " ").replace(" ", " "))


def load(rel):
    return norm(re.sub(r"\s+", " ", (P / rel).read_text()))


DRAFT = load("MANUSCRIPT_REVISED_DRAFT.md")
RISK = load("FINAL_RISK_REVIEW.md")
SUB = load("SUBMISSION_READINESS.md")
J = {}
for n in ("locus", "dfast", "dnds", "frame", "frame2", "caller_comparison",
          "locus_sensitivity", "g3_by_lesion", "g3_contig_sensitivity",
          "review_response", "seed_sensitivity", "overlap_sensitivity",
          "pf_residual_dependence"):
    J[n] = json.loads((P / f"results/{n}.json").read_text())
for n in ("two_spreads", "stage_a_evaluation", "b1_diagnosis", "rng_coupling",
          "frame2_confound"):
    J[n] = json.loads((P / f"benchmark/results/{n}.json").read_text())

fails = []


def need(doc, s, where):
    if norm(str(s)) not in doc:
        fails.append(f"{where}: {s!r} not found")


def sp(n):
    return f"{n:,}".replace(",", " ")


# ---- the four callers' RNG-free figures, which the draft's 3.7 table uses
for c in ("Pseudofinder", "DFAST", "PGAP", "Bakta"):
    lbl = "as published"
    f = J["seed_sensitivity"]["per_caller"][c][lbl]["rng_free"]
    o = f["original_rules"]
    need(DRAFT, f"{o['median_gap']:+.4f}", f"draft 3.7 gap {c}")
    need(DRAFT, f"{o['n_positive']}/{o['n_genomes']}", f"draft 3.7 positive {c}")
    lo, hi = f["genus_cluster_ci95"]
    need(DRAFT, f"[{lo:+.4f}, {hi:+.4f}]", f"draft 3.7 genus CI {c}")

# ---- G3 and its sensitivity
l = J["locus"]
need(DRAFT, sp(l["n_loci"]), "draft: scored loci")
need(DRAFT, sp(l["n_discordant"]), "draft: discordant loci")
need(DRAFT, sp(l["n_concordant_intact"]), "draft: concordant loci")
r1 = J["review_response"]["R1_genus_level"]
need(DRAFT, f"{r1['median']:+.4f}", "draft: genus-level effect")
need(DRAFT, f"{r1['n_positive']} of {r1['n_genera']}", "draft: genus counts")
r2 = J["review_response"]["R2_ci_of_reported_statistic"]
need(DRAFT, f"[{r2['lo']:+.4f}, {r2['hi']:+.4f}]", "draft: genus bootstrap CI")
r3 = J["review_response"]["R3_real_dropout"]["Pseudofinder"]
need(DRAFT, f"{r3['median']:+.4f}", "draft: Pseudofinder dropout")
ov = J["overlap_sensitivity"]["thresholds"]
lo_g3 = min(v["g3_on_matched"]["median_gap"] for v in ov.values())
need(DRAFT, f"{lo_g3:+.4f}", "draft: overlap sensitivity lower bound")
cs = J["g3_contig_sensitivity"]["arms"]["single-contig only"]
need(DRAFT, f"{cs['median_gap']:+.4f}", "draft: single-contig G3")

# ---- FRAME2
fr = J["frame2"]["arms"]
for arm in ("frameshift_high", "stop_high", "synonymous"):
    need(DRAFT, f"{fr[arm]['median_delta']:+.4f}", f"draft: frame2 {arm}")
need(DRAFT, str(fr["untouched"]["n"]), "draft: frame2 untouched n")
need(DRAFT, f"{J['frame2']['n_paired_genes']} of "
            f"{J['frame2']['n_paired_genes']}", "draft: same-orthologue count")
fc = J["frame2_confound"]["spearman_kept_vs_delta"]
need(DRAFT, f"{fc['rho']:+.3f}", "draft: frame2 confound rho")

# ---- stratification
gl = J["g3_by_lesion"]["groups"]
need(DRAFT, f"{gl['frameshift only']['median_gap']:+.4f}", "draft: frameshift stratum")
need(DRAFT, f"{gl['no frameshift']['median_gap']:+.4f}", "draft: stop stratum")
rr6 = J["review_response"]["R6_exposure"]
need(DRAFT, f"{100 * rr6['frac_classified']:.1f} %", "draft: stratification coverage")
pf_share = rr6["pseudo_side_counts"]["Pseudofinder"] / rr6["n_discordant"]
need(DRAFT, f"{100 * pf_share:.1f} %", "draft: Pseudofinder discordant share")
for doc, nm in ((DRAFT, "draft"), (RISK, "risk")):
    need(doc, f"{100 * pf_share:.1f} %", f"{nm}: PF share")

# ---- matched comparison, and the forbidden reading
r5 = J["review_response"]["R5_matched_ci"]
need(DRAFT, f"[{r5['ci95'][0]:+.4f}, {r5['ci95'][1]:+.4f}]", "draft: matched CI")
need(DRAFT, f"{100 * r5['upper_as_frac_of_G3']:.0f} %", "draft: matched CI share")

# ---- residual Pseudofinder dependence
d = J["pf_residual_dependence"]
need(RISK, sp(d["D1_reference"]["withheld_as_pseudogene"]), "risk: withheld genes")
need(RISK, f"{100 * d['D1_reference']['frac_withheld']:.1f} %", "risk: withheld frac")
need(RISK, f"{100 * d['D1_reference']['per_partner_max_frac']:.1f} %", "risk: max frac")
need(RISK, sp(d["D2_baseline"]["n_two_caller_dependent_on_pf"]), "risk: PF-dependent baseline")
need(RISK, f"{d['D2_baseline']['g3_baseline_without_pf_at_all']['median_gap']:+.4f}",
     "risk: G3 without PF baseline")
need(DRAFT, f"{d['D2_baseline']['g3_baseline_without_pf_at_all']['median_gap']:+.4f}",
     "draft: G3 without PF baseline")
if d["D2_verdict"] != "CLEARED":
    fails.append("risk: D2 should be CLEARED by the json")

# ---- two spreads
ts = J["two_spreads"]
need(DRAFT, f"{ts['corpus_total_ratio']:.1f}-fold", "draft: corpus spread")
need(DRAFT, f"{ts['median_fraction_ratio']:.1f}-fold", "draft: fraction spread")
bk = next(r for r in ts["per_caller"] if r["caller"] == "Bakta")
need(DRAFT, f"{bk['genomes_with_zero']} of {bk['genomes']}", "draft: Bakta zeros")

# ---- Stage A
b4 = J["stage_a_evaluation"]["endpoints"]["B"]
for c, v in b4.items():
    need(DRAFT, f"{v['fpr_untouched']:.3f}", f"draft: FPR {c}")
bd = J["b1_diagnosis"]
need(DRAFT, f"{bd['P3_starts_after_stop']} of {bd['n_checked']}", "draft: P3")
rngc = J["rng_coupling"]["PGAP"]
need(DRAFT, f"{rngc['pseudo_identical']} of {rngc['pseudo_compared']}",
     "draft: PGAP pseudo identical")

# ---- forbidden phrasings, in all three documents
for doc, nm in ((DRAFT, "draft"), (RISK, "risk"), (SUB, "submission")):
    low = doc.lower()
    for bad in ("statistically indistinguishable", "are equivalent",
                "shown to be equivalent", "proves equivalence"):
        if bad in low:
            fails.append(f"{nm}: forbidden equivalence phrasing {bad!r}")
    if re.search(r"corrected\s+g3\s+(effect\s+)?(is|of)\s*\+?0\.022", low):
        fails.append(f"{nm}: presents +0.022 as the corrected G3 effect")
if "absence of evidence, not equivalence" not in DRAFT:
    fails.append("draft: must state that non-significance is not equivalence")
if "not** a corrected whole-corpus G3 effect" not in DRAFT and \
        "not a corrected whole-corpus G3 effect" not in DRAFT:
    fails.append("draft: must disclaim +0.0220 as a whole-corpus effect")

print(f"checked the revised draft, the risk review and the submission decision")
if fails:
    print(f"\n{len(fails)} problems:")
    for f in fails:
        print(f"  {f}")
    sys.exit(1)
print("every claim traces to a validated output; both forbidden readings absent")
print(f"  G3 {l['verdicts']['G3']['verdict']}  genus p {r1['p']:.3g}")
print(f"  frameshift {fr['frameshift_high']['median_delta']:+.4f} vs stop "
      f"{fr['stop_high']['median_delta']:+.4f}")
print(f"  stratification covers {100 * rr6['frac_classified']:.1f}% of discordant loci")
print(f"  reference withholds {100 * d['D1_reference']['frac_withheld']:.1f}% "
      f"of partner genes (max {100 * d['D1_reference']['per_partner_max_frac']:.1f}%)")
