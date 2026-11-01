#!/usr/bin/env python3
"""POST HOC, 2026-10-30. A caller-neutral locus set, and G3 run as registered.

Fixed in LOCUS_PLAN.md before this file existed; the plan's SHA-256 is asserted
below. DNDS_PLAN.md and FRAME_PLAN.md are untouched.

G3 is quoted from DNDS_PLAN.md unchanged. It could not be run because the four
annotations share no gene identity. Once contig names are reconciled the
intervals live in one coordinate system, and overlap supplies one.

The pipeline is borrowed from analysis_dnds.py through the AST so the two
cannot drift; that module's top level never runs.
"""
import ast
import collections
import hashlib
import json
import pathlib
import sys
import tempfile
import types

import numpy as np
import pandas as pd
from scipy import stats

P = pathlib.Path(__file__).resolve().parent
D = pathlib.Path.home() / "research/datasets"
R = P / "results"
PLAN = P / "LOCUS_PLAN.md"
FROZEN = "7bd089381e4fb8f6355b4ef2aa429fb305e3521c816306e57043a9e64d3c4a80"

h = hashlib.sha256(PLAN.read_bytes()).hexdigest()
if h != FROZEN:
    sys.exit(f"LOCUS_PLAN.md changed since it was frozen\n  expected {FROZEN}\n"
             f"  found    {h}")
print(f"plan sha256 {h[:16]}… matches the frozen plan")

RECIPROCAL = 0.50       # plan section 2
MIN_DISCORDANT = 500    # L1 confirm
FAIL_DISCORDANT = 100   # L1 fail
MAX_UNJOINED = 0.40     # L1 fail
OK_UNJOINED = 0.20      # L1 confirm


def borrow(path):
    tree = ast.parse(path.read_text())
    simple = (ast.Constant, ast.List, ast.Dict, ast.Set, ast.Tuple,
              ast.Call, ast.BinOp)
    keep = [n for n in tree.body
            if isinstance(n, (ast.FunctionDef, ast.Import, ast.ImportFrom,
                              ast.ClassDef))
            or (isinstance(n, (ast.Assign, ast.AnnAssign))
                and n.value is not None and isinstance(n.value, simple))]
    mod = types.ModuleType("dnds_borrowed")
    exec(compile(ast.Module(body=keep, type_ignores=[]), str(path), "exec"),
         mod.__dict__)
    return mod


X = borrow(P / "analysis_dnds.py")
MIN_CODONS, CALLERS = X.MIN_CODONS, X.CALLERS
RNG = np.random.default_rng(20261030)
X.selftest()


def components(ivs):
    """Union-find over intervals linked by RECIPROCAL overlap of each.

    ivs: list of (contig, start, end, strand, caller, klass). Sorted by start
    within each contig+strand, so only forward neighbours need checking until
    one starts beyond the current end.
    """
    parent = list(range(len(ivs)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        a, b = find(i), find(j)
        if a != b:
            parent[b] = a

    bykey = collections.defaultdict(list)
    for i, (c, s, e, st, _, _) in enumerate(ivs):
        bykey[(c, st)].append(i)
    for key, idx in bykey.items():
        idx.sort(key=lambda i: ivs[i][1])
        for a in range(len(idx)):
            ia = idx[a]
            _, sa, ea, _, _, _ = ivs[ia]
            for b in range(a + 1, len(idx)):
                ib = idx[b]
                _, sb, eb, _, _, _ = ivs[ib]
                if sb >= ea:
                    break
                ov = min(ea, eb) - max(sa, sb)
                if ov <= 0:
                    continue
                if ov >= RECIPROCAL * (ea - sa) and ov >= RECIPROCAL * (eb - sb):
                    union(ia, ib)
    out = collections.defaultdict(list)
    for i in range(len(ivs)):
        out[find(i)].append(i)
    return list(out.values())


def main():
    probe = int(sys.argv[1]) if len(sys.argv) > 1 else 0

    d = pd.read_csv(R / "corpus523.tsv", sep="\t")
    d["genus"] = d.organism.map(X.genus_of)
    d["species"] = d.organism.map(X.species_of)
    have = d[d.genus.notna() & d.species.notna()]
    pairs = []
    for g, grp in have.groupby("genus"):
        if grp.species.nunique() < 2:
            continue
        for row in grp.itertuples():
            c = grp[grp.species != row.species].accession.sort_values()
            if len(c):
                pairs.append((row.accession, c.iloc[0], row.stratum, g))
    pairs.sort()
    if probe:
        pairs = pairs[:probe]
    print(f"\ntarget/partner pairs: {len(pairs)}")

    rows = []
    n_iv = n_unjoined = n_ambig = n_multi = n_by_order = 0
    for k, (acc, partner, stratum, gen) in enumerate(pairs, 1):
        contigs = X.genome_seqs(acc)
        pcds = X.read_fasta(D / f"n2-pseudofinder/{partner}/pf_intact.ffn")
        pfaa = D / f"n2-pseudofinder/{partner}/pf_intact.faa"
        if not contigs or not pcds or not pfaa.exists():
            continue
        resolve = X.contig_key(contigs)
        reglen = {
            "Bakta": X.gff_region_lengths(D / "n2-bakta" / acc / f"{acc}.gff3"),
            "DFAST": X.gff_region_lengths(D / "n2-dfast" / acc / "genome.gff"),
        }

        ivs = []
        n_by_order_here = 0
        for caller in CALLERS:
            # `seen` must span BOTH classes of one caller. The order-based
            # fall-back in contig_key keys off first-appearance index, so a
            # per-class dict would give the pseudogene and intact intervals of
            # the same caller different indices for the same contig and map
            # them to different contigs. analysis_dnds.py scopes it per caller
            # and this must match, or the whole locus premise breaks.
            seen = {}
            for klass, want in (("pseudo", True), ("intact", False)):
                for (c, s, e, st) in X.calls_for(acc, caller, want):
                    if c not in seen:
                        seen[c] = len(seen)
                    real, by_order = resolve(
                        c, reglen.get(caller, {}).get(c), seen[c])
                    if by_order:
                        n_by_order_here += 1
                    if real is None or e > len(contigs[real]) or e - s < 30:
                        continue
                    ivs.append((real, s, e, st, caller, klass))
        n_by_order += n_by_order_here
        if not ivs:
            continue
        n_iv += len(ivs)

        comps = components(ivs)
        # a component of one never had a partner to agree or disagree with
        n_unjoined += sum(1 for c in comps if len(c) == 1)

        keep = []
        for comp in comps:
            if len(comp) < 2:
                continue
            # count intervals per caller first, so "two intervals from one
            # caller" is counted even when the pair is pseudo+intact and the
            # locus is then dropped as ambiguous. Counting it only on the
            # surviving loci would make L3, the merging check, look better
            # than it is by hiding the worst cases.
            per_caller = collections.Counter(ivs[i][4] for i in comp)
            multi = any(v > 1 for v in per_caller.values())
            if multi:
                n_multi += 1
            verdict = {}
            for i in comp:
                caller, klass = ivs[i][4], ivs[i][5]
                if caller in verdict and verdict[caller] != klass:
                    verdict[caller] = "ambiguous"
                elif caller not in verdict:
                    verdict[caller] = klass
            if any(v == "ambiguous" for v in verdict.values()):
                n_ambig += 1
                continue
            ps = [c for c, v in verdict.items() if v == "pseudo"]
            it = [c for c, v in verdict.items() if v == "intact"]
            if ps and it:
                kind = "discordant"
            elif it and len(it) >= 2 and not ps:
                kind = "concordant_intact"
            else:
                continue
            longest = max(comp, key=lambda i: ivs[i][2] - ivs[i][1])
            keep.append((kind, ivs[longest], tuple(sorted(ps)), tuple(sorted(it))))

        if not keep:
            continue

        with tempfile.TemporaryDirectory() as td:
            qf = pathlib.Path(td) / "q.fna"
            meta = {}
            with qf.open("w") as fh:
                for n, (kind, iv, ps, it) in enumerate(keep):
                    c, s, e, st, _, _ = iv
                    seq = X.sub(contigs[c], s, e, st)
                    name = f"L{n}"
                    fh.write(f">{name}\n{seq}\n")
                    meta[name] = (seq, kind, ps, it, c, s, e)
            hits = X.best_hits(qf, pfaa)

        for name, (sid, qs, qe) in hits.items():
            pnt = pcds.get(sid)
            if pnt is None:
                continue
            cp = X.codon_pairs(meta[name][0], pnt, qs, qe)
            if cp is None:
                continue
            dn, ds, used = X.ng86(*cp)
            if dn is None or ds is None or ds <= 0:
                continue
            _, kind, ps, it, c, s, e = meta[name]
            rows.append(dict(accession=acc, genus=gen, stratum=stratum,
                             contig=c, start=s, end=e, kind=kind,
                             pseudo_callers="+".join(ps) or "-",
                             intact_callers="+".join(it) or "-",
                             omega=dn / ds, codons=used))
        if k % 10 == 0 or k == len(pairs):
            print(f"  {k}/{len(pairs)} genomes, {len(rows)} scored loci")

    t = pd.DataFrame(rows)
    if t.empty:
        sys.exit("no loci scored; nothing to report")
    t.to_csv(R / "loci.tsv", sep="\t", index=False)
    nd = int((t.kind == "discordant").sum())
    nc = int((t.kind == "concordant_intact").sum())
    unj = n_unjoined / n_iv if n_iv else 1.0
    multi_frac = n_multi / max(1, nd + nc)
    print(f"\nwrote results/loci.tsv, {len(t)} scored loci "
          f"({nd} discordant, {nc} concordant-intact)")
    print(f"  intervals {n_iv}, unjoined {n_unjoined} ({100 * unj:.0f}%), "
          f"ambiguous loci {n_ambig}, loci with a repeated caller {n_multi} "
          f"({100 * multi_frac:.1f}%)")

    out = {"plan_sha256": h, "n_loci": len(t), "n_discordant": nd,
           "n_concordant_intact": nc, "n_intervals": n_iv,
           "frac_unjoined": unj, "n_ambiguous": n_ambig,
           "frac_multi_caller": multi_frac,
           "n_contigs_resolved_by_order": n_by_order, "verdicts": {}}

    # invariants: cheap, and each one has a way of being silently false
    print("\ninvariants")
    inv = []

    def check(name, ok, detail):
        inv.append((name, bool(ok), detail))
        print(f"  {'ok ' if ok else 'XX '} {name}: {detail}")

    check("kinds are only the two defined",
          set(t.kind) <= {"discordant", "concordant_intact"},
          f"{sorted(set(t.kind))}")
    check("every discordant locus names both sides",
          not len(t[(t.kind == "discordant")
                    & ((t.pseudo_callers == "-") | (t.intact_callers == "-"))]),
          "no discordant row has an empty pseudo or intact list")
    check("no concordant-intact locus names a pseudogene caller",
          not len(t[(t.kind == "concordant_intact") & (t.pseudo_callers != "-")]),
          "none")
    check("no caller appears on both sides of one locus",
          not any(set(r.pseudo_callers.split("+")) & set(r.intact_callers.split("+"))
                  for r in t.itertuples() if r.pseudo_callers != "-"),
          "ambiguous loci were excluded, so this must hold")
    check("omega is finite and positive",
          bool(np.isfinite(t.omega).all() and (t.omega >= 0).all()),
          f"min {t.omega.min():.4f}, max {t.omega.max():.3f}")
    check("every locus meets the codon floor",
          bool((t.codons >= MIN_CODONS).all()),
          f"min codons {int(t.codons.min())}, floor {MIN_CODONS}")
    check("spans are positive and ordered",
          bool((t.end > t.start).all()), "end > start everywhere")
    check("loci are unique per genome and position",
          not t.duplicated(["accession", "contig", "start", "end"]).any(),
          "no duplicated genome/contig/span")
    out["invariants"] = [{"name": n, "ok": o, "detail": d} for n, o, d in inv]
    if not all(o for _, o, _ in inv):
        (R / "locus.json").write_text(json.dumps(out, indent=1))
        sys.exit("an invariant failed; refusing to report verdicts")

    # a per-genome diagnostic table, so the aggregate can be audited
    diag = (t.groupby("accession")
             .agg(genus=("genus", "first"), stratum=("stratum", "first"),
                  n_loci=("omega", "size"),
                  n_discordant=("kind", lambda s: int((s == "discordant").sum())),
                  n_concordant=("kind", lambda s: int((s == "concordant_intact").sum())),
                  median_omega_discordant=("omega", lambda s: np.nan),
                  median_codons=("codons", "median"))
             .reset_index())
    for i, r in diag.iterrows():
        s = t[(t.accession == r.accession) & (t.kind == "discordant")].omega
        diag.loc[i, "median_omega_discordant"] = float(np.median(s)) if len(s) else np.nan
        c = t[(t.accession == r.accession) & (t.kind == "concordant_intact")].omega
        diag.loc[i, "median_omega_concordant"] = float(np.median(c)) if len(c) else np.nan
    diag.to_csv(R / "locus_diagnostics.tsv", sep="\t", index=False)
    print(f"\nwrote results/locus_diagnostics.tsv, {len(diag)} genomes")

    def verdict(k, pred, ok, bad, val):
        v = "CONFIRMED" if ok else ("FAILED" if bad else "INCONCLUSIVE")
        out["verdicts"][k] = dict(prediction=pred, value=val, verdict=v)
        print(f"  {k}  {v:13s} {val}")

    def genome_test(a, b):
        """One difference of medians per genome, then a sign test over genomes.

        Exact zeros are DROPPED, which is what a sign test requires. A bare
        `v > 0` lumps ties in with the negatives and biases the p-value. Ties
        are plausible here because each genome's value is a median of few loci.
        """
        g = []
        for acc in set(a.accession) & set(b.accession):
            x, y = a[a.accession == acc].omega, b[b.accession == acc].omega
            if len(x) >= 3 and len(y) >= 3:
                g.append(float(np.median(x) - np.median(y)))
        nz = [v for v in g if v != 0.0]
        n_tied = len(g) - len(nz)
        if len(nz) < 5:
            return None, None, len(g), n_tied
        pos = sum(1 for v in nz if v > 0)
        p = stats.binomtest(pos, len(nz), 0.5).pvalue
        return float(np.median(g)), float(p), len(g), n_tied

    print("\nverdicts, as fixed in LOCUS_PLAN.md")
    verdict("L1", "the locus set is usable",
            nd >= MIN_DISCORDANT and unj < OK_UNJOINED,
            nd < FAIL_DISCORDANT or unj > MAX_UNJOINED,
            f"{nd} discordant loci, {100 * unj:.0f}% of intervals unjoined")

    if out["verdicts"]["L1"]["verdict"] == "FAILED":
        print("  L1 is a gate and it failed; nothing below may be read")
    else:
        dis = t[t.kind == "discordant"]
        con = t[t.kind == "concordant_intact"]
        med, p, ng, tied = genome_test(dis, con)
        if med is None:
            verdict("G3", "where callers disagree, selection takes a side",
                    False, False,
                    f"only {ng} genomes had both kinds, {tied} of them tied")
        else:
            verdict("G3", "where callers disagree, selection takes a side",
                    p < 0.05, p >= 0.05,
                    f"discordant minus concordant-intact, median over {ng} "
                    f"genomes {med:+.3f}, sign test p {p:.3g} "
                    f"({tied} tied genomes dropped)")

        pf_only = dis[(dis.pseudo_callers == "Pseudofinder")
                      & dis.intact_callers.str.contains("PGAP")]
        pg_only = dis[(dis.pseudo_callers == "PGAP")
                      & dis.intact_callers.str.contains("Pseudofinder")]
        med2, p2, ng2, tied2 = genome_test(pf_only, pg_only)
        if med2 is None:
            verdict("L2", "where Pseudofinder and PGAP disagree, one is closer",
                    False, False,
                    f"NOT TESTABLE: only {ng2} genomes carry both directions "
                    f"({tied2} tied); {len(pf_only)} loci are Pseudofinder-only "
                    f"pseudogene against PGAP intact, {len(pg_only)} the "
                    f"reverse. The asymmetry is itself the observation.")
        else:
            verdict("L2", "where Pseudofinder and PGAP disagree, one is closer",
                    p2 < 0.05, p2 >= 0.05,
                    f"PF-pseudo/PGAP-intact minus the reverse, median over "
                    f"{ng2} genomes {med2:+.3f}, sign test p {p2:.3g} "
                    f"({len(pf_only)} vs {len(pg_only)} loci)")

        verdict("L3", "overlap is a crude gene identity",
                multi_frac > 0.05, multi_frac <= 0.05,
                f"{100 * multi_frac:.1f}% of loci contain two or more "
                f"intervals from one caller")

    (R / "locus.json").write_text(json.dumps(out, indent=1))
    nbad = sum(1 for v in out["verdicts"].values() if v["verdict"] != "CONFIRMED")
    print(f"\nwrote results/locus.json - {nbad} of {len(out['verdicts'])} "
          f"not confirmed")


if __name__ == "__main__":
    main()
