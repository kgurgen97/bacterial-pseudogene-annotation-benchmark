#!/usr/bin/env python3
"""Is FRAME2's inflation systematic, or noise from a shortened alignment?

POST HOC, not in FRAME2_PLAN.md, and it does not change any verdict there.

P1 found a single simulated frameshift inflates omega by a median +0.0315.
Before that is used to threaten section 3.7, one alternative has to be ruled
out: a frameshift truncates the region DIAMOND can align, so the modified gene
may be scored on far fewer codons. omega estimated from few codons is noisy,
and noise is symmetric - it inflates and deflates equally - so it cannot by
itself produce a one-sided median shift. But it can inflate the MEAN and widen
the CI, and if the surviving alignment is tiny the estimate is not meaningful
at all.

Three things to look at:
  C1  how many codons each arm loses between the two assemblies
  C2  whether delta omega tracks the codon loss
  C3  the frameshift arm restricted to genes that keep a usable alignment

    python3 benchmark/check_frame2_confound.py
"""
import csv
import json
import pathlib

import numpy as np
from scipy import stats

P = pathlib.Path(__file__).resolve().parent.parent
FRAMESHIFT = ("M2_fs_ins1", "M3_fs_del1", "M4_fs_indel2")


def main() -> None:
    rows = list(csv.DictReader((P / "results/frame2_genes.tsv").open(),
                               delimiter="\t"))
    for r in rows:
        r["delta_omega"] = float(r["delta_omega"])
        r["cp"] = int(r["codons_parental"])
        r["cm"] = int(r["codons_modified"])
        r["lost"] = r["cp"] - r["cm"]
        r["kept"] = r["cm"] / r["cp"] if r["cp"] else 0.0

    arms = {
        "frameshift_high": [r for r in rows if r["lesion_code"] in FRAMESHIFT
                            and r["truth_tier"] == "high"],
        "stop_high": [r for r in rows if r["lesion_code"] == "M1_stop_sub"
                      and r["truth_tier"] == "high"],
        "untouched": [r for r in rows if r["lesion_code"] == "untouched"],
    }

    print("C1  aligned codons, parental against modified")
    out = {"arms": {}}
    for nm, v in arms.items():
        if not v:
            continue
        print(f"  {nm:18s} n {len(v):4d}  median codons "
              f"{np.median([r['cp'] for r in v]):6.0f} -> "
              f"{np.median([r['cm'] for r in v]):6.0f}"
              f"   median kept {np.median([r['kept'] for r in v]):.3f}"
              f"   min kept {min(r['kept'] for r in v):.3f}")
        out["arms"][nm] = {
            "n": len(v),
            "median_codons_parental": float(np.median([r["cp"] for r in v])),
            "median_codons_modified": float(np.median([r["cm"] for r in v])),
            "median_fraction_kept": float(np.median([r["kept"] for r in v])),
            "min_fraction_kept": float(min(r["kept"] for r in v))}

    fs = arms["frameshift_high"]
    print("\nC2  does delta omega track the codon loss?")
    if len(fs) >= 5:
        x = np.array([r["kept"] for r in fs])
        y = np.array([r["delta_omega"] for r in fs])
        rho, p = stats.spearmanr(x, y)
        print(f"  Spearman(fraction of codons kept, delta omega) = "
              f"{rho:+.3f}, p = {p:.3g}, n = {len(fs)}")
        print("  a strong NEGATIVE rho would mean the inflation is largest")
        print("  where the alignment shrank most, i.e. driven by the loss")
        out["spearman_kept_vs_delta"] = {"rho": float(rho), "p": float(p),
                                         "n": len(fs)}
    print("  per gene:")
    for r in sorted(fs, key=lambda r: -r["delta_omega"]):
        print(f"    {r['locus_tag']:16s} {r['lesion_code']:12s} "
              f"{r['position']:7s} codons {r['cp']:4d}->{r['cm']:4d} "
              f"({r['kept']:.2f})  omega {float(r['omega_parental']):.4f}"
              f" -> {float(r['omega_modified']):.4f}"
              f"  delta {r['delta_omega']:+.4f}")

    print("\nC3  frameshift arm restricted to genes keeping >= 50 % of codons")
    keep = [r for r in fs if r["kept"] >= 0.5]
    for nm, v in (("all frameshift", fs), (">=50% codons kept", keep)):
        if len(v) < 3:
            print(f"  {nm:22s} n {len(v)} too few")
            continue
        d = [r["delta_omega"] for r in v]
        nz = [x for x in d if x != 0]
        pos = sum(1 for x in nz if x > 0)
        pv = (stats.binomtest(pos, len(nz), 0.5).pvalue
              if len(nz) >= 5 else None)
        print(f"  {nm:22s} n {len(v):3d}  median {np.median(d):+.4f}"
              f"  mean {np.mean(d):+.4f}  {pos}/{len(nz)} positive"
              f"  p {'n/a' if pv is None else f'{pv:.3g}'}")
        out.setdefault("restricted", {})[nm] = {
            "n": len(v), "median": float(np.median(d)),
            "mean": float(np.mean(d)), "n_positive": pos,
            "p": None if pv is None else float(pv)}

    print("\nverdict")
    med_all = np.median([r["delta_omega"] for r in fs]) if fs else None
    med_keep = np.median([r["delta_omega"] for r in keep]) if len(keep) >= 3 else None
    if med_keep is not None and med_all is not None:
        if med_keep >= 0.010:
            print("  The inflation SURVIVES restricting to genes that keep most")
            print("  of their alignment, so it is not an artefact of scoring a")
            print("  short remnant. It is a property of how a frameshifted gene")
            print("  aligns, which is what section 3.7 is exposed to.")
        else:
            print("  The inflation does NOT survive the restriction: it is")
            print("  concentrated in genes whose alignment collapsed, so it is")
            print("  better described as noise from a short remnant than as a")
            print("  systematic artefact.")
    (P / "benchmark/results/frame2_confound.json").write_text(
        json.dumps(out, indent=2, sort_keys=True) + "\n")
    print("\nwrote benchmark/results/frame2_confound.json")


if __name__ == "__main__":
    main()
