#!/usr/bin/env python3
"""How much of the dN/dS gap is the seed? A sensitivity analysis, not a rewrite.

WHY THIS EXISTS
---------------
analysis_dnds.py draws its intact controls from ONE module-level generator:

    RNG = np.random.default_rng(20261028)                       # line 39
    ...
    idx = RNG.choice(len(it), len(ps), replace=False)            # line 537

shared by every caller and every genome, consumed in pair order. The classifier
fix changed len(ps) and len(it) for Bakta and DFAST, which moved the stream, so
PGAP's controls changed in 88 of 89 genomes -- identical only in pair position 1
-- while its pseudogene omegas stayed identical in all 85. Its sign test moved
from p = 0.009 to p = 0.134 without its parser changing a line.
check_rng_coupling.py demonstrates this.

That means the published per-caller sign tests are conditional on an arbitrary
stream position. This script measures how much that matters.

TWO DOCUMENTED DEVIATIONS from the frozen method in DNDS_PLAN.md
----------------------------------------------------------------
D1  No subsampling when building the omega pool. The frozen method computes
    omega only for a matched subsample of intact genes. To resample controls
    across seeds without recomputing omega 40 times, omega is computed once for
    EVERY pseudogene and EVERY intact gene. The pool is therefore deterministic
    and RNG-free; the sampling happens afterwards, on numbers already fixed.

D2  Independent streams. Controls are drawn per (accession, caller) from a
    generator seeded by blake2b(f"{accession}|{caller}|{seed}"), so no genome's
    draw can move another's, and a caller's draw cannot move a different
    caller's. This is the proposed repair, reported ALONGSIDE the original
    single-stream method rather than replacing it.

Neither DNDS_PLAN.md nor analysis_dnds.py is modified. The original method's
results stay exactly as published and recomputed.

    python3 analysis_seed_sensitivity.py --pool     # ~70 min, once
    python3 analysis_seed_sensitivity.py --sweep    # seconds
"""
from __future__ import annotations

import argparse
import ast
import collections
import csv
import hashlib
import json
import pathlib
import tempfile
import types

import numpy as np
from scipy import stats

P = pathlib.Path(__file__).resolve().parent
R = P / "results"

# A fixed, documented seed set. The first two are the seeds the project itself
# used, so the published configuration is inside the sweep rather than compared
# against from outside. The rest are 1..38, chosen before any result was seen.
SEEDS = [20261028, 20261030] + list(range(1, 39))


def borrow(path: pathlib.Path):
    """analysis_dnds.py's machinery, without running its top level."""
    tree = ast.parse(path.read_text())
    keep = [n for n in tree.body
            if isinstance(n, (ast.Import, ast.ImportFrom, ast.FunctionDef,
                              ast.ClassDef, ast.Assign, ast.AnnAssign))]
    mod = types.ModuleType("dnds_borrowed")
    mod.__dict__["__file__"] = str(path)
    exec(compile(ast.Module(body=keep, type_ignores=[]), str(path), "exec"),
         mod.__dict__)
    return mod


X = borrow(P / "analysis_dnds.py")
POOL = R / "omega_pool.tsv"


def stream(acc: str, caller: str, seed: int) -> np.random.Generator:
    """D2: a generator that depends only on this accession, caller and seed."""
    h = hashlib.blake2b(f"{acc}|{caller}|{seed}".encode(), digest_size=8)
    return np.random.default_rng(int.from_bytes(h.digest(), "big"))


# ---------------------------------------------------------------- phase 1
def build_pool() -> None:
    # analysis_dnds.py builds its pair list inside main(), so it cannot be
    # borrowed. The recomputed table records exactly which pairs it used and in
    # what order, so the pool covers the same pairs rather than a reselection.
    rows = list(csv.DictReader((R / "dnds_pairs.tsv").open(), delimiter="\t"))
    seen: dict[str, tuple] = {}
    for r in rows:
        seen.setdefault(r["accession"], (r["partner"], r["stratum"], r["genus"]))
    pairs = [(a, p, s, g) for a, (p, s, g) in seen.items()]
    print(f"omega pool over {len(pairs)} pairs, no subsampling (deviation D1)")

    out = POOL.open("w", newline="")
    w = csv.writer(out, delimiter="\t")
    w.writerow(["accession", "partner", "genus", "stratum", "caller", "call",
                "idx", "dN", "dS", "omega", "codons"])
    n_rows = 0
    for k, (acc, partner, stratum, gen) in enumerate(pairs, 1):
        contigs = X.genome_seqs(acc)
        pcds = X.read_fasta(X.D / f"n2-pseudofinder/{partner}/pf_intact.ffn")
        pfaa = X.D / f"n2-pseudofinder/{partner}/pf_intact.faa"
        if not contigs or not pcds or not pfaa.exists():
            continue
        resolve = X.contig_key(contigs)
        reglen = {
            "Bakta": X.gff_region_lengths(X.D / "n2-bakta" / acc / f"{acc}.gff3"),
            "DFAST": X.gff_region_lengths(X.D / "n2-dfast" / acc / "genome.gff"),
        }
        with tempfile.TemporaryDirectory() as td:
            qf = pathlib.Path(td) / "q.fna"
            meta = {}
            with qf.open("w") as fh:
                for caller in X.CALLERS:
                    ps = X.calls_for(acc, caller, True)
                    it = X.calls_for(acc, caller, False)
                    if not ps:
                        continue
                    # D1: every gene, both classes, no RNG anywhere
                    seen_c = {}
                    for tag, lst in (("pseudo", ps), ("intact", it)):
                        for i, (c, s, e, st) in enumerate(lst):
                            if c not in seen_c:
                                seen_c[c] = len(seen_c)
                            real, _ = resolve(
                                c, reglen.get(caller, {}).get(c), seen_c[c])
                            if real is None:
                                continue
                            if e - s < X.MIN_CODONS * 3 or e > len(contigs[real]):
                                continue
                            name = f"{caller}|{tag}|{i}"
                            seq = X.sub(contigs[real], s, e, st)
                            fh.write(f">{name}\n{seq}\n")
                            meta[name] = seq
            if not meta:
                continue
            hits = X.best_hits(qf, pfaa)
        for name, (sid, qs, qe) in hits.items():
            pnt = pcds.get(sid)
            if pnt is None:
                continue
            cp = X.codon_pairs(meta[name], pnt, qs, qe)
            if cp is None:
                continue
            dn, ds, used = X.ng86(*cp)
            if dn is None or ds is None or ds <= 0:
                continue
            caller, tag, i = name.split("|")
            w.writerow([acc, partner, gen, stratum, caller, tag, i,
                        f"{dn:.8f}", f"{ds:.8f}", f"{dn / ds:.8f}", used])
            n_rows += 1
        if k % 5 == 0 or k == len(pairs):
            out.flush()
            print(f"  {k}/{len(pairs)} pairs, {n_rows} omega values", flush=True)
    out.close()
    print(f"wrote {POOL.relative_to(P)}, {n_rows} rows")


# ---------------------------------------------------------------- phase 2
FULLCOUNTS = R / "pool_fullcounts.tsv"


def build_fullcounts(accs):
    """How many intervals each caller OFFERED, before omega attrition.

    Fidelity depends on this. The original procedure samples from the full
    intact list and only then loses genes that produce no omega -- no
    orthologue hit, no codon pairs, dS <= 0. Sampling from the survivors
    instead would draw a different number of controls from a different
    population, so seed-to-seed variance would not be comparable to the
    original's. The sweep therefore draws indices in the FULL index space and
    keeps whichever drawn indices produced an omega, which is exactly what the
    original does.
    """
    rows = []
    for acc in accs:
        for caller in X.CALLERS:
            rows.append({"accession": acc, "caller": caller,
                         "n_pseudo_full": len(X.calls_for(acc, caller, True)),
                         "n_intact_full": len(X.calls_for(acc, caller, False))})
    with FULLCOUNTS.open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t",
                           fieldnames=["accession", "caller", "n_pseudo_full",
                                       "n_intact_full"])
        w.writeheader()
        w.writerows(rows)
    return rows


def load_fullcounts():
    out = {}
    with FULLCOUNTS.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[(r["accession"], r["caller"])] = (int(r["n_pseudo_full"]),
                                                  int(r["n_intact_full"]))
    return out


def load_pool():
    """-> {(acc, caller): {"pseudo": {idx: omega}, "intact": {idx: omega}}}

    Keyed by the interval's index in calls_for's list, so a draw made in the
    full index space can be intersected with what actually produced an omega.
    """
    pool = collections.defaultdict(lambda: {"pseudo": {}, "intact": {}})
    with POOL.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pool[(r["accession"], r["caller"])][r["call"]][int(r["idx"])] = \
                float(r["omega"])
    return pool


def sign_test(vals):
    nz = [v for v in vals if v != 0.0]
    if len(nz) < 5:
        return None, None, len(vals), len(vals) - len(nz)
    pos = sum(1 for v in nz if v > 0)
    return (float(stats.binomtest(pos, len(nz), 0.5).pvalue), pos,
            len(vals), len(vals) - len(nz))


# The original clustering analysis's rules, replicated exactly so p-values are
# comparable rather than merely similar. analysis_clustering.py uses
# MIN_PER_SIDE = 5, a ONE-sided binomtest whose n is every genome (ties counted
# as non-positive, not dropped), and a cluster bootstrap that resamples whole
# clusters and recomputes the POOLED median difference.
MIN_PER_SIDE_ORIGINAL = 5
MIN_PER_SIDE_NEW = 3


def draw_for_seed(pool, full, caller, seed, subsample=True,
                  definition=None, bakta_def=None):
    """-> {accession: (pseudogene omegas, intact omegas)} for one seed.

    Returns the omega lists rather than a gap so both rule sets and the cluster
    bootstraps can be computed from the same draw.

    The original procedure samples len(ps_full) indices from range(len(it_full))
    without replacement, then computes omega and loses whatever fails.
    Reproduced by drawing in the full index space and keeping the drawn indices
    that produced an omega. `subsample=False` is the RNG-free estimator.

    `definition="strict"` restricts Bakta to the pseudogene= rows. The control
    count must then come from the STRICT pseudogene total, not the inclusive
    one: an earlier version took nps from the full count including pseudo=True,
    so the strict arm drew as many controls as the inclusive arm and its
    sampling distribution was not the strict design's.
    """
    out = {}
    for (acc, c), d in pool.items():
        if c != caller:
            continue
        ps_map, it_map = d["pseudo"], d["intact"]
        nps, nit = full.get((acc, caller), (0, 0))
        if definition == "strict" and bakta_def is not None:
            keep = bakta_def.get(acc, {})
            ps_map = {i: v for i, v in ps_map.items() if keep.get(i) == "strict"}
            # the strict design offers fewer pseudogenes, so it matches fewer
            # controls; count them over the FULL index space, as nps is
            nps = sum(1 for v in keep.values() if v == "strict")
        ps = list(ps_map.values())
        if subsample and nit > nps and nps > 0:
            g = stream(acc, caller, seed)
            drawn = g.choice(nit, nps, replace=False)
            it = [it_map[i] for i in drawn if i in it_map]
        else:
            it = list(it_map.values())
        out[acc] = (ps, it)
    return out


def genome_gaps(draw, min_per_side):
    return {acc: float(np.median(ps) - np.median(it))
            for acc, (ps, it) in draw.items()
            if len(ps) >= min_per_side and len(it) >= min_per_side}


def sign_original(gaps):
    """analysis_clustering.py's rule: one-sided, ties counted as non-positive."""
    v = list(gaps.values())
    if len(v) < 5:
        return None, None, len(v)
    pos = sum(1 for x in v if x > 0)
    p = stats.binomtest(pos, len(v), 0.5, alternative="greater").pvalue
    return float(p), pos, len(v)


def wilcoxon_greater(gaps):
    v = [x for x in gaps.values()]
    if len(v) < 5 or all(x == 0 for x in v):
        return None
    try:
        return float(stats.wilcoxon(v, alternative="greater").pvalue)
    except ValueError:
        return None


def cluster_boot(draw, genus, unit, seed=20261102, B=2000):
    """Resample whole clusters, recomputing the POOLED median difference.

    This is analysis_clustering.py's cluster_boot, not a median of per-cluster
    medians. The two are different statistics and the distinction is kept.
    """
    buckets = collections.defaultdict(lambda: ([], []))
    for acc, (ps, it) in draw.items():
        k = acc if unit == "accession" else genus.get(acc, "?")
        buckets[k][0].extend(ps)
        buckets[k][1].extend(it)
    keys = list(buckets)
    if len(keys) < 3:
        return None, None
    g = np.random.default_rng(seed)
    vals = []
    for _ in range(B):
        a, b = [], []
        for k in g.choice(keys, len(keys), replace=True):
            a.extend(buckets[k][0])
            b.extend(buckets[k][1])
        if a and b:
            vals.append(float(np.median(a) - np.median(b)))
    if not vals:
        return None, None
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def boot_ci(vals, seed=20261031, n=2000):
    if len(vals) < 5:
        return None, None
    g = np.random.default_rng(seed)
    a = np.asarray(vals, float)
    meds = [float(np.median(g.choice(a, len(a), replace=True))) for _ in range(n)]
    return float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5))


def load_genus():
    g = {}
    with (R / "dnds_pairs.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            g[r["accession"]] = r["genus"]
    return g


def load_bakta_def():
    f = R / "bakta_definition.tsv"
    if not f.exists():
        return {}
    out = collections.defaultdict(dict)
    with f.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["accession"]][int(r["idx"])] = r["definition"]
    return out


def by_genus(gaps, genus):
    """One value per genus: the median of its genomes' gaps."""
    buck = collections.defaultdict(list)
    for acc, v in gaps.items():
        buck[genus.get(acc, "?")].append(v)
    return {k: float(np.median(v)) for k, v in buck.items()}


def one_draw(draw, genus):
    """Every statistic, both rule sets, from a single draw."""
    g_orig = genome_gaps(draw, MIN_PER_SIDE_ORIGINAL)
    g_new = genome_gaps(draw, MIN_PER_SIDE_NEW)
    po, poso, no = sign_original(g_orig)
    pn, posn, nn, tied = sign_test(list(g_new.values()))
    gen = by_genus(g_orig, genus)
    return {
        # analysis_clustering.py's rules, for a comparable p-value
        "original_rules": {
            "min_per_side": MIN_PER_SIDE_ORIGINAL, "n_genomes": no,
            "median_gap": float(np.median(list(g_orig.values()))) if g_orig else None,
            "n_positive": poso, "sign_p_one_sided": po,
            "wilcoxon_p_one_sided": wilcoxon_greater(g_orig)},
        # this analysis's rules, stated separately and not mixed in
        "new_rules": {
            "min_per_side": MIN_PER_SIDE_NEW, "n_genomes": nn,
            "median_gap": float(np.median(list(g_new.values()))) if g_new else None,
            "n_positive": posn, "sign_p_two_sided": pn, "n_tied_dropped": tied},
        # The median of per-genus medians. It is not a cluster bootstrap and is
        # no substitute for one; the genus-CLUSTER bootstrap is computed
        # separately by cluster_boot and reported under its own key.
        "genus_medians_summary": {
            "n_genera": len(gen),
            "median_of_genus_medians": float(np.median(list(gen.values()))) if gen else None,
            "n_positive": sum(1 for v in gen.values() if v > 0)},
        "_gaps_original": g_orig,
    }


def summarise(label, caller, pool, full, genus, seeds, definition=None,
              bakta_def=None):
    """The RNG-free estimate with cluster bootstraps, then the seed sweep."""
    free_draw = draw_for_seed(pool, full, caller, 0, subsample=False,
                              definition=definition, bakta_def=bakta_def)
    free = one_draw(free_draw, genus)
    free.pop("_gaps_original")
    # the cluster bootstraps are the expensive part, so only on the RNG-free
    # estimator and on the project's own seed
    free["genome_cluster_ci95"] = cluster_boot(free_draw, genus, "accession")
    free["genus_cluster_ci95"] = cluster_boot(free_draw, genus, "genus")

    rows = []
    acc_o, acc_n, acc_g = [], [], []
    for s in seeds:
        d = draw_for_seed(pool, full, caller, s, subsample=True,
                          definition=definition, bakta_def=bakta_def)
        r = one_draw(d, genus)
        r.pop("_gaps_original")
        o, n, gm = r["original_rules"], r["new_rules"], r["genus_medians_summary"]
        if o["median_gap"] is not None:
            acc_o.append((o["median_gap"], o["sign_p_one_sided"],
                          o["wilcoxon_p_one_sided"]))
        if n["median_gap"] is not None:
            acc_n.append((n["median_gap"], n["sign_p_two_sided"]))
        if gm["median_of_genus_medians"] is not None:
            acc_g.append(gm["median_of_genus_medians"])
        rows.append({"caller": caller, "definition": label, "seed": s,
                     "n_genomes_original": o["n_genomes"],
                     "median_gap_original": o["median_gap"],
                     "sign_p_one_sided": o["sign_p_one_sided"],
                     "wilcoxon_p_one_sided": o["wilcoxon_p_one_sided"],
                     "n_genomes_new": n["n_genomes"],
                     "median_gap_new": n["median_gap"],
                     "sign_p_two_sided": n["sign_p_two_sided"],
                     "n_genera": gm["n_genera"],
                     "median_of_genus_medians": gm["median_of_genus_medians"]})
    if not acc_o:
        return None, rows

    def dist(vals, label_):
        v = [x for x in vals if x is not None]
        if not v:
            return None
        return {"n": len(v), "median": float(np.median(v)),
                "range": [float(min(v)), float(max(v))],
                "all_positive": bool(min(v) > 0) if label_ == "effect" else None}

    sig_o = sum(1 for _, p, _ in acc_o if p is not None and p < 0.05)
    sig_n = sum(1 for _, p in acc_n if p is not None and p < 0.05)
    return {
        "rng_free": free,
        "across_seeds": {
            "n_seeds": len(acc_o),
            "original_rules": {
                "effect": dist([m for m, _, _ in acc_o], "effect"),
                "sign_p_one_sided": dist([p for _, p, _ in acc_o], "p"),
                "wilcoxon_p_one_sided": dist([w for _, _, w in acc_o], "p"),
                "n_seeds_p_below_0.05": sig_o,
                "frac_seeds_p_below_0.05": sig_o / len(acc_o)},
            "new_rules": {
                "effect": dist([m for m, _ in acc_n], "effect"),
                "sign_p_two_sided": dist([p for _, p in acc_n], "p"),
                "n_seeds_p_below_0.05": sig_n,
                "frac_seeds_p_below_0.05": sig_n / max(len(acc_n), 1)},
            "genus_medians_summary": {"effect": dist(acc_g, "effect")},
        }}, rows


def show(caller, label, r):
    f = r["rng_free"]
    o, n = f["original_rules"], f["new_rules"]
    print(f"{caller}  [{label}]")
    print("  RNG-free: every intact gene used, nothing sampled")
    print(f"    original rules (min {o['min_per_side']}/side, one-sided): "
          f"{o['median_gap']:+.4f}  {o['n_positive']}/{o['n_genomes']}  "
          f"sign p {o['sign_p_one_sided']:.4g}  "
          f"wilcoxon p {o['wilcoxon_p_one_sided']:.4g}"
          if o["median_gap"] is not None else "    original rules: not testable")
    if n["median_gap"] is not None:
        print(f"    new rules (min {n['min_per_side']}/side, two-sided, ties "
              f"dropped): {n['median_gap']:+.4f}  "
              f"{n['n_positive']}/{n['n_genomes']}  p {n['sign_p_two_sided']:.4g}")
    for key, name in (("genome_cluster_ci95", "genome-cluster bootstrap"),
                      ("genus_cluster_ci95", "genus-cluster bootstrap")):
        ci = f.get(key)
        if ci and ci[0] is not None:
            excl = "excludes 0" if ci[0] > 0 or ci[1] < 0 else "INCLUDES 0"
            print(f"    {name} 95% CI on the pooled gap "
                  f"[{ci[0]:+.4f}, {ci[1]:+.4f}]  {excl}")
    g = f["genus_medians_summary"]
    print(f"    genus medians summary (not a bootstrap): "
          f"{g['median_of_genus_medians']:+.4f} over {g['n_genera']} genera, "
          f"{g['n_positive']} positive")

    a = r["across_seeds"]
    print(f"  across {a['n_seeds']} seeds")
    for rules, pk in (("original_rules", "sign_p_one_sided"),
                      ("new_rules", "sign_p_two_sided")):
        s = a[rules]
        e, p = s["effect"], s[pk]
        if not e:
            continue
        print(f"    {rules.replace('_', ' ')}: effect {e['median']:+.4f} "
              f"range [{e['range'][0]:+.4f}, {e['range'][1]:+.4f}] "
              f"{'all positive' if e['all_positive'] else 'SIGN VARIES'}")
        if p:
            print(f"      p median {p['median']:.4g} range "
                  f"[{p['range'][0]:.4g}, {p['range'][1]:.4g}];  "
                  f"below 0.05 in {s['n_seeds_p_below_0.05']}/{a['n_seeds']} "
                  f"seeds ({100 * s['frac_seeds_p_below_0.05']:.0f}%)")
    print("      (the fraction of significant seeds describes the stability of")
    print("       the computation, and is not itself a p-value)")
    print()


def sweep() -> None:
    pool = load_pool()
    accs = sorted({a for a, _ in pool})
    if not FULLCOUNTS.exists():
        print(f"building {FULLCOUNTS.name} for {len(accs)} accessions")
        build_fullcounts(accs)
    full = load_fullcounts()
    genus = load_genus()
    bakta_def = load_bakta_def()
    callers = sorted({c for _, c in pool})
    print(f"pool: {len(pool)} (accession, caller) cells, callers {callers}")
    print(f"seeds: {len(SEEDS)} fixed values, the project's own two first")
    print(f"draws are made in the FULL interval index space and intersected")
    print(f"with the omegas that exist, which is what the original does\n")

    report = {"seeds": SEEDS,
              "deviations": ["D1 pool built with no subsampling",
                             "D2 per-accession per-caller blake2b streams"],
              "per_caller": {}}
    rows = []
    for caller in callers:
        r, rr = summarise("as published", caller, pool, full, genus, SEEDS)
        rows += rr
        if r:
            report["per_caller"][caller] = {"as published": r}
            show(caller, "as published", r)
        if caller == "Bakta" and bakta_def:
            rs, rrs = summarise("strict pseudogene= only", caller, pool, full,
                                genus, SEEDS, definition="strict",
                                bakta_def=bakta_def)
            rows += rrs
            if rs:
                report["per_caller"][caller]["strict pseudogene= only"] = rs
                show(caller, "strict pseudogene= only", rs)

    with (R / "seed_sensitivity.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t",
                           fieldnames=["caller", "definition", "seed",
                                       "n_genomes_original",
                                       "median_gap_original",
                                       "sign_p_one_sided",
                                       "wilcoxon_p_one_sided",
                                       "n_genomes_new", "median_gap_new",
                                       "sign_p_two_sided", "n_genera",
                                       "median_of_genus_medians"])
        w.writeheader()
        w.writerows(rows)
    (R / "seed_sensitivity.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("wrote results/seed_sensitivity.json and results/seed_sensitivity.tsv")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", action="store_true")
    ap.add_argument("--sweep", action="store_true")
    a = ap.parse_args()
    if a.pool:
        build_pool()
    if a.sweep or not a.pool:
        if not POOL.exists():
            raise SystemExit("run --pool first")
        sweep()


if __name__ == "__main__":
    main()
