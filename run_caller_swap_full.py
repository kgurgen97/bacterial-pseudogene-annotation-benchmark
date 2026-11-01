#!/usr/bin/env python3
"""Rerun each study's COMPLETE pre-specified machinery under the caller swap.

The earlier pass refit only the primary models. This one runs the studies' own
`bioai.studies.association.analyze` end to end, so the stratified permutation test,
the Holm correction and the leave-one-stratum-out estimates are all recomputed, and
the verdict is produced by their rule rather than read off a p-value.

Two runs per study, on identical genome sets:
  repro  untouched features  -> must match the published statistics.json (gate)
  pf     pseudogene_fraction replaced by Pseudofinder's  -> the swap

Nothing in the studies' code, plans or seeds is modified.
"""
import csv, json, shutil, subprocess, sys
from pathlib import Path

LAB = Path.home() / "BioAI_Lab"
P = Path.home() / "research/projects/n2-pseudogene-caller-bias"
RUN = Path.home() / "research/datasets/n2-swap"
COUNTS = P / "results/pseudogene_counts.tsv"

STUDIES = {
    "h10": "study-h10v5-pseudogene-fraction-crispr-cas-complete",
    "h11": "study-h11v5-mge-per-mb-pseudogene-fraction",
    "h12": "study-h12v4-n-16s-pseudogene-fraction",
}

pf = {}
for r in csv.DictReader(COUNTS.open(), delimiter="\t"):
    try:
        ps, it = int(r["pf_pseudo"]), int(r["pf_intact"])
    except (ValueError, KeyError):
        continue
    if ps + it:
        pf[r["accession"]] = ps / (ps + it)
print(f"Pseudofinder fractions available for {len(pf)} accessions\n")


def stage(src: Path, dst: Path, swap: bool):
    """Minimal copy: everything analyze() reads, nothing else."""
    if dst.exists():
        shutil.rmtree(dst)
    (dst / "results").mkdir(parents=True)
    (dst / "metadata").mkdir()
    (dst / "figures").mkdir()
    shutil.copy2(src / "plan.json", dst / "plan.json")
    shutil.copy2(src / "metadata/manifest.tsv", dst / "metadata/manifest.tsv")
    shutil.copy2(src / "results/qc.json", dst / "results/qc.json")
    n = swapped = 0
    for d in sorted((src / "work").iterdir()):
        f = d / "features.json"
        if not f.exists():
            continue
        j = json.loads(f.read_text())
        if swap:
            v = pf.get(d.name)
            if v is None:
                continue            # keep the two runs on the same genomes
            j["pseudogene_fraction"] = round(v, 6)
            swapped += 1
        elif d.name not in pf:
            continue                # same restriction for the baseline
        o = dst / "work" / d.name
        o.mkdir(parents=True)
        (o / "features.json").write_text(json.dumps(j))
        n += 1
    return n, swapped


def run(dst: Path):
    r = subprocess.run([sys.executable, "-m", "bioai.studies.association",
                        "analyze", "--project", str(dst)],
                       cwd=LAB, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-2500:]); print(r.stderr[-2500:])
        raise SystemExit(f"analyze failed for {dst}")
    return json.loads((dst / "results/statistics.json").read_text())


def est(s):
    """The study's own primary estimate, whichever form it takes."""
    a = s["adjusted"]
    return ("OR/SD", a["or_per_sd"], a["p_predictor"]) if "or_per_sd" in a \
        else ("partial rho", a["rho"], a["p"])


def loso_range(s):
    for k in ("leave_one_stratum_out", "leave_one_phylum_out"):
        if k in s:
            v = s[k]
            e = v.get("estimates", v) if isinstance(v, dict) else v
            if isinstance(e, dict) and e:
                vals = [x for x in e.values() if isinstance(x, (int, float))]
                if vals:
                    return min(vals), max(vals), len(vals)
    return None


out = {}
for key, name in STUDIES.items():
    src = LAB / "projects" / name
    pub = json.loads((src / "results/statistics.json").read_text())
    print("=" * 78)
    print(f"{key}  {name}")
    plan = json.loads((src / "plan.json").read_text())["design"]
    print(f"  design: {plan['predictor']} -> {plan['outcome']}  "
          f"predicted {plan['direction']}  covariates {plan['covariates']}")

    rows = {}
    for tag, swap in (("repro", False), ("pf", True)):
        dst = RUN / f"{key}-{tag}"
        n, sw = stage(src, dst, swap)
        s = run(dst)
        rows[tag] = s
        label, e, p = est(s)
        lo = loso_range(s)
        print(f"\n  [{tag}] genomes staged {n}" + (f", fraction replaced in {sw}" if swap else ""))
        print(f"        n analysed      {s['n_genomes']}")
        print(f"        unadjusted p    {s['unadjusted']['p']:.4g}")
        print(f"        {label:<14s}  {e:.4f}   p {p:.4g}")
        print(f"        permutation p   {s['permutation']['p']:.4g}  "
              f"(B={s['permutation']['B']}, within {s['permutation']['within']})")
        print(f"        Holm adjusted   {s['holm']['adjusted']:.4g}")
        if lo:
            print(f"        leave-one-out   {lo[0]:+.4f} .. {lo[1]:+.4f}  over {lo[2]} strata")
        print(f"        VERDICT         {s['verdict']}")

    # gate
    pl, pe, pp = est(pub)
    rl, re_, rp = est(rows["repro"])
    gate = (abs(pe - re_) < 1e-6 and abs(pp - rp) < 1e-9
            and pub["verdict"] == rows["repro"]["verdict"]
            and abs(pub["permutation"]["p"] - rows["repro"]["permutation"]["p"]) < 1e-9)
    print(f"\n  GATE vs published: {'PASS' if gate else 'DIFFERS'}"
          f"   published {pl} {pe:.4f} p {pp:.4g} perm {pub['permutation']['p']:.4g} "
          f"{pub['verdict']}  (n {pub['n_genomes']})")
    if not gate:
        print("    -> the staged rerun is not bit-identical to the published run; "
              "differences below are read against the repro column, not the paper")

    a, b = rows["repro"], rows["pf"]
    al, ae, ap = est(a); bl, be, bp = est(b)
    print(f"\n  SWAP EFFECT   {al}: {ae:+.4f} -> {be:+.4f}   "
          f"({(be-ae)/abs(ae)*100:+.1f}% )" if ae else "")
    print(f"                permutation p: {a['permutation']['p']:.4g} -> {b['permutation']['p']:.4g}")
    print(f"                Holm adjusted: {a['holm']['adjusted']:.4g} -> {b['holm']['adjusted']:.4g}")
    print(f"                VERDICT      : {a['verdict']} -> {b['verdict']}"
          f"   {'UNCHANGED' if a['verdict'] == b['verdict'] else '*** CHANGED ***'}")

    out[key] = {
        "design": plan, "gate_pass": bool(gate),
        "published": {"n": pub["n_genomes"], "estimate": pe, "p": pp,
                      "permutation_p": pub["permutation"]["p"],
                      "holm_adjusted": pub["holm"]["adjusted"], "verdict": pub["verdict"]},
        "repro": {"n": a["n_genomes"], "estimate": ae, "p": ap,
                  "permutation_p": a["permutation"]["p"],
                  "holm_adjusted": a["holm"]["adjusted"], "verdict": a["verdict"],
                  "loso": loso_range(a)},
        "pseudofinder": {"n": b["n_genomes"], "estimate": be, "p": bp,
                         "permutation_p": b["permutation"]["p"],
                         "holm_adjusted": b["holm"]["adjusted"], "verdict": b["verdict"],
                         "loso": loso_range(b)},
        "verdict_changed": a["verdict"] != b["verdict"],
        "estimate_shift_pct": (be - ae) / abs(ae) * 100 if ae else None,
    }
    print()

(P / "results").mkdir(exist_ok=True)
json.dump(out, open(P / "results/caller_swap_full.json", "w"), indent=1)
print("=" * 78)
print(f"{'study':6s} {'gate':6s} {'verdict repro':>16s} {'verdict pf':>16s} {'perm p repro':>13s} {'perm p pf':>11s}")
for k, v in out.items():
    print(f"{k:6s} {'PASS' if v['gate_pass'] else 'DIFF':6s} "
          f"{v['repro']['verdict']:>16s} {v['pseudofinder']['verdict']:>16s} "
          f"{v['repro']['permutation_p']:>13.4g} {v['pseudofinder']['permutation_p']:>11.4g}")
print(f"\nverdicts changed by the caller swap: "
      f"{sum(v['verdict_changed'] for v in out.values())} of {len(out)}")
print(f"wrote {P/'results/caller_swap_full.json'}")
