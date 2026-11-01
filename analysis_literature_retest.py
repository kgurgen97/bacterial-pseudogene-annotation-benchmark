#!/usr/bin/env python3
"""The canonical claims under test, exactly as LITERATURE_RETEST_PLAN.md specifies.

Plan sha256: a576cda95ff60a368d6f1726d8117e9909cff795d5890c67efa75ce118e1c5a9

C1 (primary)   pseudogene fraction rises as genome size falls. Predicted: NEGATIVE
               Spearman. Label-free, all genomes.
C2 (secondary) host-associated exceed free-living. Predicted: host > free.
               74 genomes with tax_id/species-level labels.

Each claim is tested independently for every caller whose counts are complete, on the
identical genome set. Runs on whatever is available and is meant to be rerun when the
remaining caller finishes.
"""
import json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

P = Path.home() / "research/projects/n2-pseudogene-caller-bias"
R = P / "results"
ALPHA = 0.05

# ---------------------------------------------------------------- assemble
pg = pd.read_csv(R / "pseudogene_counts.tsv", sep="\t")
cor = pd.read_csv(R / "corpus_traits.tsv", sep="\t", dtype=str)
size = pd.to_numeric(cor.genome_size, errors="coerce")
mb = pd.to_numeric(cor.get("genome_size_mb"), errors="coerce") * 1e6
cor["size_bp"] = size.fillna(mb)
cor["gc"] = pd.to_numeric(cor.gc_percent, errors="coerce")

d = pg[["accession", "pgap_pseudo", "pgap_cds", "pf_pseudo", "pf_intact"]].copy()
d = d.merge(cor[["accession", "size_bp", "gc", "phylum", "iso_madin", "join_level"]],
            on="accession", how="left")

callers = {}
d["frac_pgap"] = d.pgap_pseudo / (d.pgap_pseudo + d.pgap_cds)
d["cnt_pgap"] = d.pgap_pseudo
callers["PGAP"] = ("frac_pgap", "cnt_pgap")
d["frac_pfinder"] = d.pf_pseudo / (d.pf_pseudo + d.pf_intact)
d["cnt_pfinder"] = d.pf_pseudo
callers["Pseudofinder"] = ("frac_pfinder", "cnt_pfinder")

# The plan, section 5: "Every test is run independently for each of the four
# callers on the identical genome set." A caller part-way through its run breaks
# that, so admission requires FULL coverage of this corpus - not a loose floor.
# The earlier floor of 300 of 363 would have admitted Bakta 63 genomes short and
# computed its associations on a different, non-random subset from the others'.
NEEDED = set(d.accession)
for nm, f, pcol, ccol in (("Bakta", "bakta_counts.tsv", "bakta_pseudo", "bakta_cds"),
                          ("DFAST", "dfast_counts.tsv", "dfast_disrupted", "dfast_cds")):
    p = R / f
    if not p.exists():
        continue
    t = pd.read_csv(p, sep="\t")
    t = t[[c for c in ("accession", pcol, ccol) if c in t.columns]].dropna()
    have = NEEDED & set(t.accession)
    if len(have) < len(NEEDED):
        print(f"skipping {nm}: covers {len(have)} of {len(NEEDED)} corpus genomes "
              f"({len(NEEDED) - len(have)} missing); rerun when its run completes")
        continue
    d = d.merge(t, on="accession", how="left")
    d[f"frac_{nm.lower()}"] = d[pcol] / (d[pcol] + d[ccol])
    d[f"cnt_{nm.lower()}"] = d[pcol]
    callers[nm] = (f"frac_{nm.lower()}", f"cnt_{nm.lower()}")

# Having admitted only fully-covering callers, every fraction column must be
# complete. If this ever fires, some genome produced a count but not a usable
# fraction, and the "identical genome set" guarantee is gone.
for nm, (fcol, _) in callers.items():
    missing = int(d[fcol].isna().sum())
    assert missing == 0, (f"{nm} has {missing} genomes without a fraction; the "
                          f"callers would not be compared on one genome set")

d["per_mb_pgap"] = d.pgap_pseudo / (d.size_bp / 1e6)
d["per_mb_pfinder"] = d.pf_pseudo / (d.size_bp / 1e6)
for nm in ("bakta", "dfast"):
    if f"cnt_{nm}" in d:
        d[f"per_mb_{nm}"] = d[f"cnt_{nm}"] / (d.size_bp / 1e6)

print(f"genomes: {len(d)}   callers complete: {list(callers)}")
print(f"phyla: {dict(d.phylum.fillna('unknown').value_counts())}")

out = {"callers": list(callers), "n_genomes": int(len(d))}


def spear(x, y):
    ok = x.notna() & y.notna()
    if ok.sum() < 10:
        return np.nan, np.nan, int(ok.sum())
    r, p_ = stats.spearmanr(x[ok], y[ok])
    return float(r), float(p_), int(ok.sum())


def partial_spear(x, y, z):
    ok = x.notna() & y.notna() & z.notna()
    if ok.sum() < 12:
        return np.nan, np.nan
    rk = lambda v: stats.rankdata(v[ok])
    Z = np.column_stack([np.ones(ok.sum()), rk(z)])
    rx = rk(x) - Z @ np.linalg.lstsq(Z, rk(x), rcond=None)[0]
    ry = rk(y) - Z @ np.linalg.lstsq(Z, rk(y), rcond=None)[0]
    r, p_ = stats.pearsonr(rx, ry)
    return float(r), float(p_)


# ---------------------------------------------------------------- C1
print("\n" + "=" * 86)
print("C1  pseudogene fraction vs genome size.  PREDICTED: negative")
print("=" * 86)
c1 = {}
for scope, sub in [("pooled", d)] + [(ph, g) for ph, g in d.groupby(d.phylum.fillna("unknown"))]:
    print(f"\n--- {scope}  (n={len(sub)})")
    print(f"  {'caller':14s} {'measure':10s} {'rho':>8s} {'p':>11s} "
          f"{'partial rho|GC':>15s} {'p':>11s} {'n':>5s}  dir")
    for nm, (fcol, ccol) in callers.items():
        for measure, col in (("fraction", fcol), ("count", ccol),
                             ("per Mb", f"per_mb_{nm.lower().replace('pseudofinder','pfinder')}")):
            if col not in sub:
                continue
            r, p_, n = spear(sub[col], sub.size_bp)
            pr, pp = partial_spear(sub[col], sub.size_bp, sub.gc)
            if r != r:
                continue
            ok = "as predicted" if (r < 0 and p_ < ALPHA) else \
                 ("OPPOSITE" if (r > 0 and p_ < ALPHA) else "n.s.")
            print(f"  {nm:14s} {measure:10s} {r:+8.3f} {p_:11.3g} "
                  f"{pr:+15.3f} {pp:11.3g} {n:5d}  {ok}")
            c1[f"{scope}|{nm}|{measure}"] = dict(rho=r, p=p_, partial_rho=pr,
                                                 partial_p=pp, n=n, verdict=ok)

print("\n--- C1 PRIMARY ENDPOINT (fraction, pooled)")
hits = [(nm, c1[f"pooled|{nm}|fraction"]) for nm in callers
        if f"pooled|{nm}|fraction" in c1]
nsig = sum(1 for _, v in hits if v["rho"] < 0 and v["p"] < ALPHA)
nopp = sum(1 for _, v in hits if v["rho"] > 0 and v["p"] < ALPHA)
mags = [abs(v["rho"]) for _, v in hits]
print(f"  callers tested                     : {len(hits)}")
print(f"  significant in predicted direction : {nsig}")
print(f"  significant in the OPPOSITE direction: {nopp}")
print(f"  effect-size spread (max/min |rho|) : "
      f"{max(mags)/min(mags):.2f}  ({min(mags):.3f} to {max(mags):.3f})")
for nm, v in hits:
    print(f"    {nm:14s} rho {v['rho']:+.3f}  p {v['p']:.3g}  -> {v['verdict']}")
out["C1"] = {"per_test": c1, "n_callers": len(hits), "n_as_predicted": nsig,
             "n_opposite": nopp,
             "effect_spread": float(max(mags) / min(mags)) if mags else None}

# ---------------------------------------------------------------- C2
print("\n" + "=" * 86)
print("C2  host-associated vs free-living.  PREDICTED: host > free")
print("=" * 86)
HOST = lambda s: isinstance(s, str) and s.startswith("host")
FREE = lambda s: isinstance(s, str) and s.split("_")[0] in ("soil", "water", "sediment", "rock")
d["lifestyle"] = np.where(d.iso_madin.map(HOST), "host",
                   np.where(d.iso_madin.map(FREE), "free", None))

c2 = {}
for setname, levels in (("primary (tax_id/species)", {"tax_id", "species"}),
                        ("exploratory (+genus)", {"tax_id", "species", "genus"})):
    sub = d[d.join_level.isin(levels) & d.lifestyle.notna()]
    nh = int((sub.lifestyle == "host").sum()); nf = int((sub.lifestyle == "free").sum())
    print(f"\n--- {setname}: {len(sub)} genomes ({nh} host, {nf} free-living)")
    if nh < 8 or nf < 8:
        print("  too few, skipped"); continue
    print(f"  {'caller':14s} {'host med':>10s} {'free med':>10s} {'MWU p':>11s} "
          f"{'rank-biserial':>14s}  dir")
    for nm, (fcol, _) in callers.items():
        a = sub.loc[sub.lifestyle == "host", fcol].dropna()
        b = sub.loc[sub.lifestyle == "free", fcol].dropna()
        if len(a) < 8 or len(b) < 8:
            continue
        u = stats.mannwhitneyu(a, b, alternative="two-sided")
        rbc = 2 * u.statistic / (len(a) * len(b)) - 1
        ok = "as predicted" if (a.median() > b.median() and u.pvalue < ALPHA) else \
             ("OPPOSITE" if (a.median() < b.median() and u.pvalue < ALPHA) else "n.s.")
        print(f"  {nm:14s} {a.median():10.4f} {b.median():10.4f} {u.pvalue:11.3g} "
              f"{rbc:+14.3f}  {ok}")
        c2[f"{setname}|{nm}"] = dict(host_median=float(a.median()),
                                     free_median=float(b.median()),
                                     p=float(u.pvalue), rbc=float(rbc),
                                     n_host=len(a), n_free=len(b), verdict=ok)
    if setname.startswith("primary"):
        hs = [(nm, v) for k, v in c2.items() if k.startswith(setname)
              for nm in [k.split("|")[1]]]
        ns = sum(1 for _, v in hs if v["host_median"] > v["free_median"] and v["p"] < ALPHA)
        mg = [abs(v["rbc"]) for _, v in hs if v["rbc"] == v["rbc"] and v["rbc"] != 0]
        print(f"\n  C2 PRIMARY ENDPOINT: {ns} of {len(hs)} callers significant in the "
              f"predicted direction")
        if mg:
            print(f"  effect-size spread (max/min |rank-biserial|): "
                  f"{max(mg)/min(mg):.2f}  ({min(mg):.3f} to {max(mg):.3f})")
        out["C2"] = {"n_callers": len(hs), "n_as_predicted": ns,
                     "effect_spread": float(max(mg) / min(mg)) if mg else None}
out["C2_detail"] = c2

# ---------------------------------------------------------------- caller spread
print("\n" + "=" * 86)
print("For context: how far apart are the callers on the same genomes?")
print(f"  {'caller':14s} {'median fraction':>16s} {'median count':>13s}  "
      f"{'IQR of fraction':>22s}")
for nm, (fcol, ccol) in callers.items():
    f_ = d[fcol].dropna(); c_ = d[ccol].dropna()
    print(f"  {nm:14s} {f_.median():16.4f} {c_.median():13.0f}  "
          f"{f_.quantile(.25):10.4f}-{f_.quantile(.75):<10.4f}")
fr = [d[f].median() for f, _ in callers.values()]
print(f"  ratio of the largest to the smallest median fraction: {max(fr)/min(fr):.1f}x")
out["median_fraction_spread"] = float(max(fr) / min(fr))

d.to_csv(R / "literature_retest.tsv", sep="\t", index=False)
json.dump(out, open(R / "literature_retest.json", "w"), indent=1, default=str)
print(f"\nwrote results/literature_retest.{{tsv,json}}")
