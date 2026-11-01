#!/usr/bin/env python3
"""POST HOC, 2026-10-29. Is the selection signal in section 3.7 a frame artefact?

Not covered by DNDS_PLAN.md, which stays frozen. Fixed in FRAME_PLAN.md before
this file existed, and this script asserts that plan's SHA-256 before computing.

Candidates are aligned to their orthologue in the ORTHOLOGUE's frame, so past a
frameshift the codons are near random with respect to the partner and dN rises
with no selection involved. This splits each protein alignment at its first gap
column and recomputes dN/dS on the codons before it.

The heavy machinery is borrowed from analysis_dnds.py rather than copied, so the
two cannot drift. Its top level is never executed: the module is parsed, the
definitions and simple constants are taken, and nothing else runs.
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

P = pathlib.Path.home() / "research/projects/n2-pseudogene-caller-bias"
D = pathlib.Path.home() / "research/datasets"
R = P / "results"
PLAN = P / "FRAME_PLAN.md"
FROZEN = "05c154af6e1dfb037458ed76436e00200b5932bc1e03549a47141bc2da31c2d6"

h = hashlib.sha256(PLAN.read_bytes()).hexdigest()
if h != FROZEN:
    sys.exit(f"FRAME_PLAN.md has changed since it was frozen\n"
             f"  expected {FROZEN}\n  found    {h}")
print(f"plan sha256 {h[:16]}… matches the frozen plan")


def borrow(path):
    """Import a sibling's definitions without running its top level."""
    tree = ast.parse(path.read_text())
    # Definitions, imports, and plain data assignments. AnnAssign is included
    # because an annotated constant such as "SYN_CACHE: dict = {}" is not an
    # Assign node and was silently dropped, which the self-test caught.
    simple = (ast.Constant, ast.List, ast.Dict, ast.Set, ast.Tuple,
              ast.Call, ast.BinOp)
    keep = []
    for n in tree.body:
        if isinstance(n, (ast.FunctionDef, ast.Import, ast.ImportFrom,
                          ast.ClassDef)):
            keep.append(n)
        elif isinstance(n, (ast.Assign, ast.AnnAssign)) and \
                n.value is not None and isinstance(n.value, simple):
            keep.append(n)
    mod = types.ModuleType("dnds_borrowed")
    mod.__dict__["__file__"] = str(path)
    exec(compile(ast.Module(body=keep, type_ignores=[]), str(path), "exec"),
         mod.__dict__)
    return mod


X = borrow(P / "analysis_dnds.py")
MIN_CODONS = X.MIN_CODONS
CALLERS = X.CALLERS
# the same seed as the published run, so the same intact controls are drawn
RNG = np.random.default_rng(20261028)
print(f"borrowed the pipeline: MIN_CODONS {MIN_CODONS}, callers {CALLERS}")

X.selftest()


def split_codon_pairs(cand_nt, partner_nt, qstart, qend):
    """As analysis_dnds.codon_pairs, but split at the first gap column.

    Returns (pre, post, first_gap_col, n_total) where pre and post are each a
    pair of codon lists. A sequence whose alignment has no gap at all puts
    everything in pre, which is correct: there is no damage to sit behind.
    """
    Seq = X.Seq
    span = cand_nt[min(qstart, qend) - 1:max(qstart, qend)]
    if len(span) < MIN_CODONS * 3:
        return None
    if qstart > qend:
        span = str(Seq(span).reverse_complement())
    best = None
    for frame in (0, 1, 2):
        trimmed = span[frame:]
        trimmed = trimmed[:len(trimmed) - len(trimmed) % 3]
        if len(trimmed) < MIN_CODONS * 3:
            continue
        prot = str(Seq(trimmed).translate())
        ref = str(Seq(partner_nt[:len(partner_nt) - len(partner_nt) % 3]).translate())
        try:
            aln = X.ALIGNER.align(prot, ref)[0]
        except Exception:                                      # noqa: BLE001
            continue
        if best is None or aln.score > best[0]:
            best = (aln.score, aln, trimmed, partner_nt)
    if best is None:
        return None
    _, aln, qnt, snt = best
    qa, sa = aln[0], aln[1]

    first_gap = None
    for col, (a, b) in enumerate(zip(qa, sa)):
        if a == "-" or b == "-":
            first_gap = col
            break

    pre1, pre2, post1, post2 = [], [], [], []
    all1, all2 = [], []
    qi = si = 0
    for col, (a, b) in enumerate(zip(qa, sa)):
        if a != "-" and b != "-":
            x, y = qnt[qi * 3:qi * 3 + 3], snt[si * 3:si * 3 + 3]
            if len(x) == 3 and len(y) == 3:
                all1.append(x); all2.append(y)
                if first_gap is None or col < first_gap:
                    pre1.append(x); pre2.append(y)
                else:
                    post1.append(x); post2.append(y)
        if a != "-":
            qi += 1
        if b != "-":
            si += 1
    return ((pre1, pre2), (post1, post2), first_gap, len(pre1) + len(post1),
            (all1, all2))


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
    dropped = collections.Counter()
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
        with tempfile.TemporaryDirectory() as td:
            qf = pathlib.Path(td) / "q.fna"
            meta = {}
            with qf.open("w") as fh:
                for caller in CALLERS:
                    ps = X.calls_for(acc, caller, True)
                    it = X.calls_for(acc, caller, False)
                    if not ps:
                        continue
                    if len(it) > len(ps):
                        idx = RNG.choice(len(it), len(ps), replace=False)
                        it = [it[i] for i in idx]
                    seen = {}
                    for tag, lst in (("pseudo", ps), ("intact", it)):
                        for i, (c, s, e, st) in enumerate(lst):
                            if c not in seen:
                                seen[c] = len(seen)
                            real, _ = resolve(
                                c, reglen.get(caller, {}).get(c), seen[c])
                            if real is None:
                                continue
                            if e - s < MIN_CODONS * 3 or e > len(contigs[real]):
                                continue
                            name = f"{caller}|{tag}|{i}"
                            fh.write(f">{name}\n{X.sub(contigs[real], s, e, st)}\n")
                            meta[name] = X.sub(contigs[real], s, e, st)
            if not meta:
                continue
            hits = X.best_hits(qf, pfaa)

        for name, (sid, qs, qe) in hits.items():
            pnt = pcds.get(sid)
            if pnt is None:
                continue
            got = split_codon_pairs(meta[name], pnt, qs, qe)
            if got is None:
                continue
            pre, post, fgap, ntot, whole = got
            caller, tag, _ = name.split("|")

            def omega(cp, floor=MIN_CODONS):
                if len(cp[0]) < floor:
                    return None, len(cp[0])
                dn, ds, used = X.ng86(cp[0], cp[1], min_codons=floor)
                if dn is None or ds is None or ds <= 0:
                    return None, used or len(cp[0])
                return dn / ds, used

            w_pre, n_pre = omega(pre)
            w_post, n_post = omega(post)

            # F2 as FRAME_PLAN.md registers it: the first 15 gap-free codon
            # pairs against the last 15, within the same alignment. The first
            # implementation used the pre/post-gap split instead, which is F1's
            # own split, so it was neither the registered test nor independent
            # of gap detection. Both are reported; see FRAME_RESULTS.md.
            w_first = w_last = None
            a1, a2 = whole
            if len(a1) >= 30:                      # 15 + 15, non-overlapping
                w_first, _ = omega((a1[:15], a2[:15]), floor=15)
                w_last, _ = omega((a1[-15:], a2[-15:]), floor=15)

            if w_pre is None:
                dropped[f"{caller}|{tag}"] += 1
            rows.append(dict(accession=acc, partner=partner, genus=gen,
                             stratum=stratum, caller=caller, call=tag,
                             omega_pre=w_pre, n_pre=n_pre,
                             omega_post=w_post, n_post=n_post,
                             omega_first15=w_first, omega_last15=w_last,
                             first_gap=fgap if fgap is not None else -1,
                             n_total=ntot))
        if k % 10 == 0 or k == len(pairs):
            print(f"  {k}/{len(pairs)} genomes, {len(rows)} sequences")

    t = pd.DataFrame(rows)
    if t.empty:
        sys.exit("no sequences; nothing to report")
    t.to_csv(R / "frame_pairs.tsv", sep="\t", index=False)
    print(f"\nwrote results/frame_pairs.tsv, {len(t)} rows")

    full = json.loads((R / "dnds.json").read_text())["per_caller"]
    out = {"plan_sha256": h, "n_rows": len(t), "per_caller": {}}

    print(f"\n{'caller':14s} {'n pre ps':>9s} {'n pre in':>9s} "
          f"{'pre ps':>8s} {'pre in':>8s} {'pre gap':>8s} {'full gap':>9s} "
          f"{'kept':>6s} {'p':>10s}")
    for c in CALLERS:
        a = t[(t.caller == c) & (t.call == "pseudo")].omega_pre.dropna()
        b = t[(t.caller == c) & (t.call == "intact")].omega_pre.dropna()
        if len(a) < 5 or len(b) < 5:
            print(f"{c:14s}  too few sequences ({len(a)} vs {len(b)})")
            continue
        gap = float(np.median(a) - np.median(b))
        p = stats.mannwhitneyu(a, b, alternative="greater").pvalue
        fullgap = full[c]["gap"]
        kept = gap / fullgap if fullgap else float("nan")
        print(f"{c:14s} {len(a):9d} {len(b):9d} {np.median(a):8.3f} "
              f"{np.median(b):8.3f} {gap:+8.3f} {fullgap:+9.3f} "
              f"{kept:6.2f} {p:10.3g}")
        out["per_caller"][c] = dict(n_pseudo=len(a), n_intact=len(b),
                                    median_pseudo=float(np.median(a)),
                                    median_intact=float(np.median(b)),
                                    gap=gap, full_gap=fullgap,
                                    fraction_kept=kept, p=float(p))

    verdicts = {}

    def verdict(k, pred, ok, bad, val):
        v = "CONFIRMED" if ok else ("FAILED" if bad else "INCONCLUSIVE")
        verdicts[k] = dict(prediction=pred, value=val, verdict=v)
        print(f"  {k}  {v:13s} {val}")

    print("\nverdicts, as fixed in FRAME_PLAN.md")

    # F4 first: it is the gate
    dg = out["per_caller"].get("DFAST", {}).get("gap")
    if dg is None:
        verdict("F4", "the split must not manufacture signal", False, False,
                "DFAST yielded too few sequences to check the gate")
    else:
        verdict("F4", "the split must not manufacture signal",
                dg < 0.010, dg >= 0.010,
                f"DFAST pre gap {dg:+.3f} against its full {full['DFAST']['gap']:+.3f}")

    kept = {c: out["per_caller"][c]["fraction_kept"]
            for c in ("PGAP", "Pseudofinder") if c in out["per_caller"]}
    if len(kept) == 2:
        ok = all(v >= 0.5 for v in kept.values())
        bad = any(v < 0.25 for v in kept.values())
        verdict("F1", "the gap survives on pre-damage codons", ok, bad,
                "fraction of the full gap kept: "
                + ", ".join(f"{c} {v:.2f}" for c, v in kept.items()))
    else:
        verdict("F1", "the gap survives on pre-damage codons", False, False,
                "PGAP or Pseudofinder missing from the pre set")

    # F2 as registered: the last 15 gap-free codon pairs against the first 15,
    # paired within each alignment, so it needs no gap detection.
    rises, npair = {}, {}
    for tag in ("pseudo", "intact"):
        s = t[(t.call == tag)].dropna(subset=["omega_first15", "omega_last15"])
        npair[tag] = len(s)
        if len(s) >= 5:
            rises[tag] = float(np.median(s.omega_last15 - s.omega_first15))
    if len(rises) == 2:
        diff = rises["pseudo"] - rises["intact"]
        verdict("F2", "omega rises along the alignment faster in calls",
                diff < 0.02, diff >= 0.02,
                f"median within-sequence rise, last 15 codon pairs minus "
                f"first 15: calls {rises['pseudo']:+.3f} (n {npair['pseudo']}), "
                f"controls {rises['intact']:+.3f} (n {npair['intact']}), "
                f"difference {diff:+.3f}")
    else:
        verdict("F2", "omega rises along the alignment faster in calls",
                False, False, "too few sequences with 30 gap-free codon pairs")

    # the first implementation, kept and labelled, because it was reported
    old = {}
    for tag in ("pseudo", "intact"):
        s = t[t.call == tag]
        a, b = s.omega_pre.dropna(), s.omega_post.dropna()
        if len(a) >= 5 and len(b) >= 5:
            old[tag] = float(np.median(b) - np.median(a))
    if len(old) == 2:
        out["F2_first_implementation"] = dict(
            rise_calls=old["pseudo"], rise_controls=old["intact"],
            difference=old["pseudo"] - old["intact"],
            note="pre/post-gap split, not the registered first15/last15 "
                 "contrast, and not independent of gap detection")
        print(f"  F2 as first implemented (a deviation, not the plan): "
              f"difference {old['pseudo'] - old['intact']:+.3f}")

    gp = t[(t.call == "pseudo") & (t.first_gap >= 0)].first_gap
    gi = t[(t.call == "intact") & (t.first_gap >= 0)].first_gap
    if len(gp) >= 5 and len(gi) >= 5:
        verdict("F3", "the split captures damage",
                np.median(gi) > np.median(gp), np.median(gi) <= np.median(gp),
                f"median first-gap column: controls {np.median(gi):.0f}, "
                f"calls {np.median(gp):.0f}")
    else:
        verdict("F3", "the split captures damage", False, False,
                "too few gapped alignments")

    nogap = int((t.first_gap < 0).sum())
    print(f"\nalignments with no gap at all: {nogap} of {len(t)}")
    print("sequences whose pre set fell below the 30-codon floor:")
    for k2, v in sorted(dropped.items()):
        print(f"  {k2}: {v}")

    out["verdicts"] = verdicts
    out["no_gap_alignments"] = nogap
    out["dropped_below_floor"] = dict(dropped)
    (R / "frame.json").write_text(json.dumps(out, indent=1))
    nfail = sum(1 for v in verdicts.values() if v["verdict"] != "CONFIRMED")
    print(f"\nwrote results/frame.json - {nfail} of {len(verdicts)} not confirmed")


if __name__ == "__main__":
    main()
