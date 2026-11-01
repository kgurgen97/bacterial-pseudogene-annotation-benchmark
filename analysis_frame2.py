#!/usr/bin/env python3
"""Does aligning a frameshifted gene mechanically inflate omega? Paired test.

Analysis as fixed in FRAME2_PLAN.md (sha256 asserted below). FRAME_PLAN.md is
untouched and the original frame analysis stays void under its own F4 gate.

The design in one line: take the same gene with and without a KNOWN simulated
frameshift, run both through the identical pipeline against the identical
orthologue, and look at the difference. Everything except the lesion is held
constant, so the difference is the inflation attributable to the lesion and to
how the aligner renders it.

    python3 analysis_frame2.py [--probe]
"""
from __future__ import annotations

import argparse
import ast
import bisect
import collections
import csv
import hashlib
import json
import pathlib
import sys
import tempfile
import types

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from scipy import stats

P = pathlib.Path(__file__).resolve().parent
R = P / "results"
B = P / "benchmark"
ASM = pathlib.Path.home() / "research/datasets/n2-benchmark/assemblies"

PLAN = P / "FRAME2_PLAN.md"
PLAN_SHA = "a926ec3c7e579959940d2ad1adf1bb6a6f86e95cc7f0710d85bceb7227755e83"

FRAMESHIFT = ("M2_fs_ins1", "M3_fs_del1", "M4_fs_indel2")
STOP = ("M1_stop_sub",)
SYNONYMOUS = ("C1_synonymous",)
BAR_N2 = 0.010
BAR_LARGE = 0.030
BAR_SMALL = 0.010


def assert_plan():
    h = hashlib.sha256(PLAN.read_bytes()).hexdigest()
    if h != PLAN_SHA:
        sys.exit(f"FRAME2_PLAN.md changed.\n  expected {PLAN_SHA}\n"
                 f"  found    {h}\nA frozen plan is never edited because of a "
                 f"result.")


def borrow(path):
    tree = ast.parse(path.read_text())
    keep = [n for n in tree.body
            if isinstance(n, (ast.Import, ast.ImportFrom, ast.FunctionDef,
                              ast.ClassDef, ast.Assign, ast.AnnAssign))]
    m = types.ModuleType("b")
    m.__dict__["__file__"] = str(path)
    exec(compile(ast.Module(body=keep, type_ignores=[]), str(path), "exec"),
         m.__dict__)
    return m


X = borrow(P / "analysis_dnds.py")


class Lift:
    """Rebuilt from the saved liftover table; altlen recovered exactly."""

    def __init__(self, rows):
        rows = sorted(rows)
        self.starts, self.ends, self.altlen, self.cum = [], [], [], []
        prev = 0
        for s, e, c in rows:
            self.starts.append(s)
            self.ends.append(e)
            self.altlen.append((e - s) + (c - prev))
            self.cum.append(c)
            prev = c

    def delta_before(self, pos):
        i = bisect.bisect_right(self.ends, pos) - 1
        return self.cum[i] if i >= 0 else 0

    def to_modified(self, pos):
        i = bisect.bisect_right(self.starts, pos) - 1
        if i >= 0 and self.starts[i] <= pos < self.ends[i]:
            k = pos - self.starts[i]
            if k < self.altlen[i]:
                return self.starts[i] + self.delta_before(self.starts[i]) + k
            return None
        return pos + self.delta_before(pos)


def load_lift(rid):
    p = B / f"results/liftover/{rid}.tsv"
    by = collections.defaultdict(list)
    with p.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            by[r["contig"]].append((int(r["parental_pos0"]),
                                    int(r["parental_end0"]),
                                    int(r["cum_delta"])))
    return {c: Lift(v) for c, v in by.items()}


def load_fa(p):
    return {r.id: str(r.seq).upper() for r in SeqIO.parse(str(p), "fasta")}


def gene_seq(contigs, cg, s1, e1, strand, delta=0):
    """The gene's nucleotides, in its own orientation."""
    s = contigs[cg][s1 - 1:e1 + delta]
    return str(Seq(s).reverse_complement()) if strand < 0 else s


def omegas(named: dict[str, str], pfaa, pcds):
    """-> {name: omega} through exactly the project's pipeline."""
    out = {}
    with tempfile.TemporaryDirectory() as td:
        qf = pathlib.Path(td) / "q.fna"
        with qf.open("w") as fh:
            for k, v in named.items():
                fh.write(f">{k}\n{v}\n")
        hits = X.best_hits(qf, pfaa)
    for name, (sid, qs, qe) in hits.items():
        pnt = pcds.get(sid)
        if pnt is None:
            continue
        cp = X.codon_pairs(named[name], pnt, qs, qe)
        if cp is None:
            continue
        dn, ds, used = X.ng86(*cp)
        if dn is None or ds is None or ds <= 0:
            continue
        out[name] = (dn / ds, used, sid)
    return out


def sign_test(v):
    nz = [x for x in v if x != 0.0]
    if len(nz) < 5:
        return None, None, len(nz)
    pos = sum(1 for x in nz if x > 0)
    return float(stats.binomtest(pos, len(nz), 0.5).pvalue), pos, len(nz)


def boot_ci(v, seed=20261105, n=4000):
    if len(v) < 5:
        return None, None
    g = np.random.default_rng(seed)
    a = np.asarray(v, float)
    m = [float(np.median(g.choice(a, len(a), replace=True))) for _ in range(n)]
    return float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def main() -> None:
    assert_plan()
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", action="store_true",
                    help="one replicate only, for validation")
    a = ap.parse_args()
    print(f"plan {PLAN_SHA[:16]} ok")

    # the only simulated genome with a dN/dS partner
    partner = {}
    with (R / "dnds_pairs.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            partner.setdefault(r["accession"], r["partner"])
    truth = list(csv.DictReader((B / "results/ground_truth.tsv").open(),
                                delimiter="\t"))
    elig = list(csv.DictReader((B / "results/eligible_genes.tsv").open(),
                               delimiter="\t"))
    accs = sorted({r["accession"] for r in truth} & set(partner))
    print(f"simulated genomes with a partner: {accs}")
    if not accs:
        sys.exit("no simulated genome has a dN/dS partner; NOT TESTABLE")

    reps = sorted({r["replicate"] for r in truth
                   if r["accession"] in accs})
    if a.probe:
        reps = reps[:1]
    print(f"replicates: {reps}")

    rows = []
    for rid in reps:
        acc = rid.split("__")[0]
        pt = partner[acc]
        pfaa = X.D / f"n2-pseudofinder/{pt}/pf_intact.faa"
        pcds = X.read_fasta(X.D / f"n2-pseudofinder/{pt}/pf_intact.ffn")
        par = load_fa(ASM / f"{acc}__C0_parental.fna")
        mod = load_fa(ASM / f"{rid}.fna")
        lift = load_lift(rid)
        les = {r["locus_tag"]: r for r in truth if r["replicate"] == rid}
        genes = [g for g in elig if g["accession"] == acc]
        # N1 needs untouched genes, but all of them would be thousands of
        # alignments; a fixed deterministic slice is enough for an exactness
        # gate and is documented as such
        untouched = [g for g in genes if g["locus_tag"] not in les][:150]
        use = [g for g in genes if g["locus_tag"] in les] + untouched
        print(f"\n{rid}: {len(use)} genes "
              f"({len(use) - len(untouched)} lesioned, {len(untouched)} untouched)")

        q_par, q_mod, meta = {}, {}, {}
        for g in use:
            tag, cg = g["locus_tag"], g["contig"]
            s1, e1, st = int(g["start1"]), int(g["end1"]), int(g["strand"])
            l = les.get(tag)
            delta = int(l["delta_bp"]) if l else 0
            lf = lift.get(cg)
            if lf is None:
                continue
            try:
                if st >= 0:
                    m0 = lf.to_modified(s1 - 1)
                    mseq = mod[cg][m0:m0 + (e1 - s1 + 1) + delta]
                else:
                    m1 = lf.to_modified(e1 - 1)
                    n = (e1 - s1 + 1) + delta
                    mseq = str(Seq(mod[cg][m1 - n + 1:m1 + 1]).reverse_complement())
            except (KeyError, TypeError):
                continue
            pseq = gene_seq(par, cg, s1, e1, st)
            if not pseq or not mseq:
                continue
            q_par[tag] = pseq
            q_mod[tag] = mseq
            meta[tag] = {"lesion_code": l["lesion_code"] if l else "untouched",
                         "truth_tier": l["truth_tier"] if l else "untouched",
                         "position": l["position"] if l else "",
                         "n_codons": g["n_codons"]}

        op = omegas(q_par, pfaa, pcds)
        om = omegas(q_mod, pfaa, pcds)
        both = sorted(set(op) & set(om))
        print(f"  omega computed for {len(op)} parental, {len(om)} modified, "
              f"{len(both)} in both")
        for tag in both:
            wp, cp_, sp = op[tag]
            wm, cm_, sm = om[tag]
            rows.append({"replicate": rid, "accession": acc, "locus_tag": tag,
                         **meta[tag],
                         "omega_parental": wp, "omega_modified": wm,
                         "delta_omega": wm - wp,
                         "same_orthologue": int(sp == sm),
                         "codons_parental": cp_, "codons_modified": cm_})

    with (R / "frame2_genes.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote results/frame2_genes.tsv with {len(rows)} paired genes")

    def arm(codes, tier="high"):
        return [r for r in rows if r["lesion_code"] in codes
                and (tier is None or r["truth_tier"] == tier)]

    out = {"plan_sha256": PLAN_SHA, "n_paired_genes": len(rows),
           "replicates": reps, "verdicts": {}}

    def record(key, name, vals, verdict, value):
        out["verdicts"][key] = {"name": name, "verdict": verdict,
                                "value": value, "n": len(vals)}
        print(f"  {key:4s} {verdict:14s} {value}")

    print("\nverdicts, as fixed in FRAME2_PLAN.md")

    # ---- N1 gate: untouched genes must be identical
    un = [r for r in rows if r["lesion_code"] == "untouched"]
    nonzero = [r for r in un if r["delta_omega"] != 0.0]
    n1_ok = len(un) > 0 and not nonzero
    record("N1", "untouched genes give delta exactly 0",
           un, "CONFIRMED" if n1_ok else "FAILED",
           f"{len(un) - len(nonzero)}/{len(un)} exactly zero"
           + ("" if n1_ok else f"; worst |delta| "
              f"{max(abs(r['delta_omega']) for r in nonzero):.4g}"))
    if not n1_ok:
        print("\n  N1 is a gate and it failed: the pipeline is not deterministic")
        print("  across assemblies, so nothing else here may be interpreted.")
        (R / "frame2.json").write_text(json.dumps(out, indent=2,
                                                  sort_keys=True) + "\n")
        sys.exit(1)

    # ---- N2 synonymous
    syn = arm(SYNONYMOUS, tier=None)
    sv = [r["delta_omega"] for r in syn]
    med_s = float(np.median(sv)) if sv else None
    record("N2", "a synonymous change does not inflate omega", syn,
           "CONFIRMED" if med_s is not None and abs(med_s) < BAR_N2
           else ("FAILED" if med_s is not None else "NOT_TESTABLE"),
           f"median delta {med_s:+.4f} against a bar of {BAR_N2}"
           if med_s is not None else "no synonymous genes paired")

    # ---- P1 primary
    fs = arm(FRAMESHIFT)
    fv = [r["delta_omega"] for r in fs]
    p, pos, n = sign_test(fv)
    lo, hi = boot_ci(fv)
    if len(fv) < 5 or p is None:
        record("P1", "a frameshift inflates omega mechanically", fs,
               "NOT_TESTABLE", f"only {len(fv)} paired frameshift genes")
        med_f = None
    else:
        med_f = float(np.median(fv))
        if med_f >= BAR_LARGE:
            v = "CONFIRMED_LARGE"
        elif med_f < BAR_SMALL:
            v = "CONFIRMED_SMALL"
        else:
            v = "INCONCLUSIVE"
        record("P1", "a frameshift inflates omega mechanically", fs, v,
               f"median delta {med_f:+.4f}, {pos}/{n} positive, "
               f"sign test p {p:.3g}, CI95 [{lo:+.4f}, {hi:+.4f}]")

    # ---- P2 specificity to frame change
    st = arm(STOP)
    stv = [r["delta_omega"] for r in st]
    med_t = float(np.median(stv)) if stv else None
    if med_f is None or med_t is None:
        record("P2", "the inflation is specific to frame change", st,
               "NOT_TESTABLE", "one arm has no paired genes")
    else:
        record("P2", "the inflation is specific to frame change", st,
               "CONFIRMED" if med_f > med_t else "FAILED",
               f"frameshift {med_f:+.4f} against premature stop "
               f"{med_t:+.4f} (n {len(stv)})")

    out["arms"] = {}
    for nm, sel in (("frameshift_high", fs), ("stop_high", st),
                    ("synonymous", syn), ("untouched", un)):
        v = [r["delta_omega"] for r in sel]
        out["arms"][nm] = {
            "n": len(v),
            "median_delta": float(np.median(v)) if v else None,
            "mean_delta": float(np.mean(v)) if v else None,
            "n_positive": sum(1 for x in v if x > 0),
            "n_zero": sum(1 for x in v if x == 0.0),
            "max_abs": float(max(abs(x) for x in v)) if v else None}
    print("\narms")
    for nm, s in out["arms"].items():
        if s["n"]:
            print(f"  {nm:18s} n {s['n']:4d}  median {s['median_delta']:+.4f}"
                  f"  mean {s['mean_delta']:+.4f}  positive {s['n_positive']}"
                  f"  exactly zero {s['n_zero']}")

    same = sum(r["same_orthologue"] for r in rows)
    out["same_orthologue_fraction"] = same / len(rows)
    print(f"\nthe same orthologue was chosen in both assemblies for "
          f"{same}/{len(rows)} genes ({100 * same / len(rows):.1f}%)")

    (R / "frame2.json").write_text(json.dumps(out, indent=2,
                                              sort_keys=True) + "\n")
    print("wrote results/frame2.json")


if __name__ == "__main__":
    main()
