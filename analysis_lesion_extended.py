#!/usr/bin/env python3
"""Caller-independent lesion evidence on EXTENDED intervals.

Executes LESION_HARMONISATION_PLAN.md, whose SHA-256 is asserted below.

WHY. analysis_lesion_from_sequence.py recovered only 18 of 413 known
frameshifts (4.4%). The mechanism is known: a caller that detects a frameshift
reports the TRUNCATED feature, so the post-lesion fragment lies outside the
annotated interval and the two-frame signature cannot be observed. The interval
was the limiting factor, not the detector. Each span is therefore extended by
max(300 bp, half the orthologue length) before the evidence is collected.

WHAT IS NEW BESIDES THE EXTENSION. A biological NEGATIVE control. The earlier
analysis measured recall and never measured specificity, so it could not say
whether its classifier labels intact genes as lesioned. Here every
concordant-intact locus -- genes all four callers call intact -- is classified
through the identical path, and the false-lesion rate among them is a gate.

Nothing here recomputes omega, reassigns an orthologue for the dN/dS analysis,
or re-annotates anything. Outputs are isolated under results/lesion_extended/.

    python3 analysis_lesion_extended.py [--targets N]
"""
from __future__ import annotations

import argparse
import ast
import collections
import csv
import hashlib
import json
import pathlib
import subprocess
import sys
import tempfile
import types

import numpy as np
import pandas as pd
from Bio.Seq import Seq
from scipy import stats

P = pathlib.Path(__file__).resolve().parent
R = P / "results"
OUT = R / "lesion_extended"
D = pathlib.Path.home() / "research/datasets"
REFDIR = D / "n2-refseq-reference"

PLAN = P / "LESION_HARMONISATION_PLAN.md"
PLAN_SHA = "27552b9dc974d8d4407a2a4ea66e9bbc0e37417e8374c8c771ece129211ae11c"

MIN_HSP_AA = 20
MIN_MARGIN = 300
MARGIN_FRAC = 0.5
SEG_OVERLAP_MAX = 20        # aa of subject overlap allowed between frames
FRAG_GAP_MIN = 30           # aa gap in the subject for `fragmentation`
STOP_COVER_MAX = 0.90
TRUNC_COVER_MAX = 0.50
INTACT_COVER_MIN = 0.90
N_SYNTH = 100
SYNTH_POS_FRAC = 0.40
THREADS = 4
SEED = 20261031
FLOOR_START_GB, FLOOR_ABORT_GB = 10, 8

GATES = {"G1_frameshift_recall": 0.40, "G2_false_frameshift": 0.10,
         "G3_false_lesion": 0.10, "G4_decided": 0.60,
         "G5_internal_stop_recall": 0.70, "G6_synthetic_positive": 0.80}
LESION_CLASSES = ("frameshift", "internal_stop", "truncation", "fragmentation")


def borrow(path):
    tree = ast.parse(path.read_text())
    keep = [n for n in tree.body
            if isinstance(n, (ast.Import, ast.ImportFrom, ast.FunctionDef,
                              ast.ClassDef, ast.Assign, ast.AnnAssign))]
    m = types.ModuleType(path.stem + "_b")
    m.__dict__["__file__"] = str(path)
    exec(compile(ast.Module(body=keep, type_ignores=[]), str(path), "exec"),
         m.__dict__)
    return m


S = borrow(P / "analysis_reference_sensitivity.py")
L, X = S.L, S.X
LFS = borrow(P / "analysis_lesion_from_sequence.py")      # for dfast_labels


def host_free_gb():
    out = subprocess.run(["df", "--output=avail", "-BG", "/mnt/c"],
                         capture_output=True, text=True, check=True).stdout
    return int("".join(c for c in out.splitlines()[-1] if c.isdigit()))


def guard(stage):
    free = host_free_gb()
    if free < FLOOR_ABORT_GB:
        sys.exit(f"EMERGENCY STOP at {stage}: {free} GB < {FLOOR_ABORT_GB}")
    return free


def partner_faa(partner):
    faa, sent = REFDIR / f"{partner}.faa", REFDIR / f"{partner}.json"
    if not (faa.exists() and sent.exists()):
        return None
    try:
        if json.loads(sent.read_text()).get("version") != 2:
            return None
    except Exception:
        return None
    return faa


FMT = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart",
       "send", "slen", "bitscore"]


def blastx(query_fa, db_faa, workdir, all_hsps=True):
    """-> {query: [hsp dicts]}, best subject only, optionally every HSP."""
    workdir.mkdir(parents=True, exist_ok=True)
    dmnd = workdir / "db.dmnd"
    if not dmnd.exists():
        subprocess.run(["diamond", "makedb", "--in", str(db_faa), "-d",
                        str(dmnd), "--quiet", "--threads", str(THREADS)],
                       check=True)
    out = workdir / "hits.tsv"
    cmd = ["diamond", "blastx", "-q", str(query_fa), "-d", str(dmnd),
           "-o", str(out), "--quiet", "--threads", str(THREADS),
           "--max-target-seqs", "1", "--evalue", "1e-5",
           "--outfmt", "6", *FMT]
    if all_hsps:
        cmd += ["--max-hsps", "0"]
    subprocess.run(cmd, check=True)
    per = collections.defaultdict(list)
    for line in out.read_text().splitlines():
        f = line.split("\t")
        if len(f) < len(FMT):
            continue
        per[f[0]].append({
            "subject": f[1], "pident": float(f[2]), "aln": int(f[3]),
            "qstart": int(f[4]), "qend": int(f[5]), "sstart": int(f[6]),
            "send": int(f[7]), "slen": int(f[8]), "bits": float(f[9])})
    # keep only HSPs against the single best subject (highest bitscore)
    outd = {}
    for q, hs in per.items():
        best = max(hs, key=lambda h: h["bits"])["subject"]
        outd[q] = [h for h in hs if h["subject"] == best]
    return outd


def frame_of(h):
    return (min(h["qstart"], h["qend"]) - 1) % 3


def reverse_hsp(h):
    return h["qstart"] > h["qend"]


def seg_overlap(a, b):
    lo = max(min(a["sstart"], a["send"]), min(b["sstart"], b["send"]))
    hi = min(max(a["sstart"], a["send"]), max(b["sstart"], b["send"]))
    return max(0, hi - lo + 1)


def premature_stop(span, hsps):
    """-> (has_stop, fraction of the subject covered before the first stop).

    Translated in the frame and orientation of the strongest HSP, starting at
    that HSP's query start and running no further than the projected end of the
    gene, so that stops in imported neighbouring sequence cannot be counted.
    """
    h = max(hsps, key=lambda x: x["bits"])
    nt = span
    if reverse_hsp(h):
        nt = str(Seq(span).reverse_complement())
        qs = len(span) - max(h["qstart"], h["qend"]) + 1
    else:
        qs = min(h["qstart"], h["qend"])
    start = qs - 1
    # project the gene's nucleotide end from the subject length and where the
    # alignment sits within the subject
    s_begin = min(h["sstart"], h["send"])
    proj_end = start + 3 * (h["slen"] - s_begin + 1)
    end = min(len(nt), max(proj_end, start + 3))
    seg = nt[start:end]
    seg = seg[:len(seg) - len(seg) % 3]
    if len(seg) < 3:
        return False, None
    prot = str(Seq(seg).translate())
    i = prot.find("*")
    if i < 0:
        return False, None
    covered = (s_begin - 1 + i) / h["slen"] if h["slen"] else None
    return True, covered


def classify(span, hsps):
    """-> (class, detail dict). Rules applied in the plan's declared order."""
    if not hsps:
        return "no_alignment", {}
    strong = [h for h in hsps if h["aln"] >= MIN_HSP_AA]
    if not strong:
        return "ambiguous", {"reason": f"no HSP >= {MIN_HSP_AA} aa"}
    slen = strong[0]["slen"]
    # "total subject coverage" is the UNION of covered subject residues. A
    # first version used the max-minus-min SPAN, which counts an unaligned gap
    # between two HSPs as covered and so understates truncation and
    # overstates intactness.
    iv = sorted((min(h["sstart"], h["send"]), max(h["sstart"], h["send"]))
                for h in strong)
    merged, covered = [], 0
    for lo, hi in iv:
        if merged and lo <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    for lo, hi in merged:
        covered += hi - lo + 1
    cov = covered / slen if slen else 0.0
    det = {"n_hsp": len(hsps), "n_strong": len(strong), "slen": slen,
           "subject_coverage": float(cov),
           "frames": sorted({frame_of(h) for h in strong})}

    # 1 frameshift: two frames covering non-overlapping subject segments
    byf = collections.defaultdict(list)
    for h in strong:
        byf[frame_of(h)].append(h)
    if len(byf) >= 2:
        reps = [max(v, key=lambda h: h["bits"]) for v in byf.values()]
        reps.sort(key=lambda h: -h["bits"])
        for i in range(len(reps)):
            for j in range(i + 1, len(reps)):
                if seg_overlap(reps[i], reps[j]) <= SEG_OVERLAP_MAX:
                    det["frameshift_pair_overlap_aa"] = seg_overlap(reps[i],
                                                                    reps[j])
                    return "frameshift", det
        det["reason"] = "multiple frames but the segments overlap"

    # 2 internal stop
    has_stop, at = premature_stop(span, strong)
    det["stop_at_subject_frac"] = at
    if has_stop and at is not None and at < STOP_COVER_MAX:
        return "internal_stop", det

    # 3 truncation
    if cov < TRUNC_COVER_MAX and not has_stop:
        return "truncation", det

    # 4 fragmentation: same frame, non-overlapping, separated in the subject
    f0 = sorted(byf, key=lambda k: -sum(h["bits"] for h in byf[k]))[0]
    same = sorted(byf[f0], key=lambda h: min(h["sstart"], h["send"]))
    for a, b in zip(same, same[1:]):
        gap = min(b["sstart"], b["send"]) - max(a["sstart"], a["send"])
        if gap >= FRAG_GAP_MIN:
            det["fragment_gap_aa"] = int(gap)
            return "fragmentation", det

    # 5 intact
    if cov >= INTACT_COVER_MIN and not has_stop:
        return "intact", det
    return "ambiguous", det


def write_fa(path, seqs):
    with open(path, "w") as fh:
        for k, v in seqs.items():
            fh.write(f">{k}\n{v}\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-locus",
                    default="results/reference_sensitivity/per_locus.tsv")
    ap.add_argument("--out-json", default="results/lesion_extended.json")
    ap.add_argument("--targets", type=int, default=0)
    # re-score the gates from already-written per-locus evidence, without
    # repeating any search. Used when the gate LOGIC is corrected but the
    # classifications are unchanged.
    ap.add_argument("--from-existing", action="store_true")
    a = ap.parse_args()

    got = hashlib.sha256(PLAN.read_bytes()).hexdigest()
    print(f"plan sha256: {got[:16]}…")
    if got != PLAN_SHA:
        sys.exit(f"FROZEN PLAN CHANGED\n  expected {PLAN_SHA}\n  found    {got}")
    free0 = guard("start")
    print(f"host C: free {free0} GB (from /mnt/c)")
    if free0 < FLOOR_START_GB:
        sys.exit(f"refusing to start: {free0} < {FLOOR_START_GB} GB")
    OUT.mkdir(parents=True, exist_ok=True)

    src = P / a.per_locus
    in_sha = hashlib.sha256(src.read_bytes()).hexdigest()
    print(f"input {src.relative_to(P)}  sha256 {in_sha[:16]}…")
    df = pd.read_csv(src, sep="\t")
    df = df[df.control.notna()]
    pairs = df[["accession", "partner"]].drop_duplicates().values.tolist()
    if a.targets:
        pairs = pairs[:a.targets]
    print(f"{len(df)} loci over {len(pairs)} (query, partner) pairs")

    print("\nrecovering locus strands from a 0.50 partition rebuild")
    strands = S.locus_strands([p[0] for p in pairs])
    print(f"  {len(strands)} representative intervals recovered")

    rng = np.random.default_rng(SEED)
    rows, synth_rows, audit = [], [], collections.Counter()

    if a.from_existing:
        out = pd.read_csv(OUT / "per_locus_extended.tsv", sep="\t")
        syn = pd.read_csv(OUT / "controls.tsv", sep="\t")
        print(f"\nre-scoring gates from existing evidence: {len(out)} loci, "
              f"{len(syn)} control rows (no search re-run)")
        pairs = out[["accession"]].drop_duplicates().values.tolist()
        return report(out, syn, pairs, got, in_sha, audit, free0, a)

    for i, (acc, partner) in enumerate(pairs, 1):
        faa = partner_faa(partner)
        contigs = X.genome_seqs(acc)
        if faa is None or not contigs:
            print(f"  [{i}/{len(pairs)}] {acc}: skipped (reference or contigs)")
            continue
        sub = df[df.accession == acc]
        base, meta = {}, {}
        for n, r in enumerate(sub.itertuples()):
            st = strands.get((r.accession, r.contig, r.start, r.end))
            if st is None or r.contig not in contigs:
                continue
            if r.end > len(contigs[r.contig]):
                continue
            seq = X.sub(contigs[r.contig], r.start, r.end, st)
            if seq:
                k = f"L{n}"
                base[k], meta[k] = seq, (r, st)
        if not base:
            continue

        # DFAST's own class per feature: (contig, start, end, class) 4-tuples.
        # A first version unpacked five values and matched on exact
        # coordinates; the arity error was swallowed by a broad except and the
        # exact match would have failed anyway, because a locus representative
        # interval is not a DFAST feature's coordinates. Labels are therefore
        # assigned by >= 50 % reciprocal overlap, as elsewhere in this project,
        # and a genuine failure is now visible rather than counted.
        dfl = LFS.dfast_labels(acc)
        audit["dfast_features"] += len(dfl)

        def dfast_for(contig, s, e):
            best, bov = None, 0
            qlen = e - s + 1
            for cg, s2, e2, cl in dfl:
                if cg != contig:
                    continue
                ov = min(e, e2) - max(s, s2) + 1
                if ov <= 0:
                    continue
                if ov >= 0.5 * qlen and ov >= 0.5 * (e2 - s2 + 1) and ov > bov:
                    best, bov = cl, ov
            return best

        with tempfile.TemporaryDirectory() as td:
            td = pathlib.Path(td)
            # pass 1: unextended, one HSP, to learn the subject length
            write_fa(td / "p1.fna", base)
            p1 = blastx(td / "p1.fna", faa, td / "db", all_hsps=False)

            # pass 2: extended spans
            ext, extmeta = {}, {}
            for k, (r, st) in meta.items():
                hs = p1.get(k)
                if not hs:
                    rows.append({"accession": acc, "genus": r.genus,
                                 "kind": r.kind, "contig": r.contig,
                                 "start": r.start, "end": r.end,
                                 "lesion_class": "no_alignment",
                                 "margin": None,
                                 "dfast_class": dfast_for(r.contig, r.start,
                                                          r.end)})
                    continue
                slen = hs[0]["slen"]
                m = max(MIN_MARGIN, int(np.ceil(MARGIN_FRAC * 3 * slen)))
                s2 = max(1, int(r.start) - m)
                e2 = min(len(contigs[r.contig]), int(r.end) + m)
                seq = X.sub(contigs[r.contig], s2, e2, st)
                if not seq:
                    continue
                ext[k], extmeta[k] = seq, (r, st, m, s2, e2)
            if ext:
                write_fa(td / "p2.fna", ext)
                p2 = blastx(td / "p2.fna", faa, td / "db", all_hsps=True)
                for k, (r, st, m, s2, e2) in extmeta.items():
                    cls, det = classify(ext[k], p2.get(k, []))
                    rows.append({
                        "accession": acc, "genus": r.genus, "kind": r.kind,
                        "contig": r.contig, "start": r.start, "end": r.end,
                        "span": int(r.end) - int(r.start) + 1,
                        "margin": m, "ext_start": s2, "ext_end": e2,
                        "lesion_class": cls,
                        "subject_coverage": det.get("subject_coverage"),
                        "n_frames": len(det.get("frames", [])),
                        "stop_at": det.get("stop_at_subject_frac"),
                        "dfast_class": dfast_for(r.contig, r.start, r.end)})

            # controls: synthetic positive and negative on intact loci
            cand = [k for k, (r, st) in meta.items()
                    if r.kind == "concordant_intact" and len(base[k]) >= 300]
            if cand:
                take = list(rng.permutation(cand))[:max(1, N_SYNTH // len(pairs))]
                pos, neg = {}, {}
                for k in take:
                    s = base[k]
                    cut = int(len(s) * SYNTH_POS_FRAC)
                    pos[k] = s[:cut] + s[cut + 1:]        # 1 bp deletion
                    neg[k] = s
                if pos:
                    write_fa(td / "sp.fna", pos)
                    write_fa(td / "sn.fna", neg)
                    hp = blastx(td / "sp.fna", faa, td / "db", all_hsps=True)
                    hn = blastx(td / "sn.fna", faa, td / "db", all_hsps=True)
                    for k in pos:
                        cp, _ = classify(pos[k], hp.get(k, []))
                        cn, _ = classify(neg[k], hn.get(k, []))
                        synth_rows.append({"accession": acc, "locus": k,
                                           "arm": "deleted_1bp",
                                           "lesion_class": cp})
                        synth_rows.append({"accession": acc, "locus": k,
                                           "arm": "unmodified",
                                           "lesion_class": cn})
        print(f"  [{i}/{len(pairs)}] {acc}: {len(ext)} extended spans "
              f"classified", flush=True)
        guard(acc)

    out = pd.DataFrame(rows)
    syn = pd.DataFrame(synth_rows)
    out.to_csv(OUT / "per_locus_extended.tsv", sep="\t", index=False)
    syn.to_csv(OUT / "controls.tsv", sep="\t", index=False)
    print(f"\n{len(out)} loci classified, {len(syn)} control rows")
    return report(out, syn, pairs, got, in_sha, audit, free0, a)


def report(out, syn, pairs, got, in_sha, audit, free0, a):
    """Scoring and gates, separated so they can be re-run without searching."""
    dis = out[out.kind == "discordant"]
    con = out[out.kind == "concordant_intact"]

    def wilson(k, n):
        if not n:
            return None, None
        lo, hi = stats.binomtest(int(k), int(n)).proportion_ci(0.95,
                                                              method="wilson")
        return float(lo), float(hi)

    print("\nCLASS DISTRIBUTION")
    tab = {}
    for kind, g in out.groupby("kind"):
        vc = g.lesion_class.value_counts()
        tab[kind] = {str(k): int(v) for k, v in vc.items()}
        print(f"  {kind}  (n={len(g)})")
        for k, v in vc.items():
            print(f"    {str(k):16s} {v:6d}  ({100 * v / len(g):5.1f} %)")

    # ---- G1, G2, G5: against DFAST's published class
    lab = dis[dis.dfast_class.notna()]
    fs = lab[lab.dfast_class.str.contains("frameshift", na=False)]
    ist = lab[lab.dfast_class == "internal_stop_codon"]
    g1 = float((fs.lesion_class == "frameshift").mean()) if len(fs) else None
    g2 = float((ist.lesion_class == "frameshift").mean()) if len(ist) else None
    g5 = float((ist.lesion_class == "internal_stop").mean()) if len(ist) else None
    # ---- G3: false lesion on the biological negative control
    g3 = float(con.lesion_class.isin(LESION_CLASSES).mean()) if len(con) else None
    # ---- G4: the plan defines "decided" as NOT ambiguous, which includes
    # `intact` -- finding a locus intact is a decision, and an informative one.
    # A first version counted only the four lesion classes, which is stricter
    # than the frozen criterion; the criterion is not changed to match the
    # code, the code is changed to match the criterion. The lesion-only
    # fraction is still reported, as a separate diagnostic.
    undecided = ("ambiguous", "no_alignment")
    g4 = float((~dis.lesion_class.isin(undecided)).mean()) if len(dis) else None
    g4_lesion_only = (float(dis.lesion_class.isin(LESION_CLASSES).mean())
                      if len(dis) else None)
    # ---- G6: synthetic positive control
    sp = syn[syn.arm == "deleted_1bp"]
    sn = syn[syn.arm == "unmodified"]
    g6 = float((sp.lesion_class == "frameshift").mean()) if len(sp) else None

    print("\nEXTERNAL VALIDATION against DFAST's published class")
    print(f"  loci with a DFAST class            {len(lab)}")
    print(f"  DFAST frameshift loci              {len(fs)}")
    print(f"  DFAST internal-stop loci           {len(ist)}")
    if len(lab):
        ct = pd.crosstab(lab.dfast_class, lab.lesion_class)
        ct.to_csv(OUT / "crosstab_vs_dfast.tsv", sep="\t")
        print(ct.to_string())

    print("\nCONTROLS")
    if len(sp):
        print(f"  synthetic, 1 bp deleted   n {len(sp):4d}  "
              f"frameshift {100 * g6:5.1f} %")
        for k, v in sp.lesion_class.value_counts().items():
            print(f"      {str(k):16s} {v}")
    if len(sn):
        fl = float(sn.lesion_class.isin(LESION_CLASSES).mean())
        print(f"  synthetic, unmodified     n {len(sn):4d}  "
              f"any lesion {100 * fl:5.1f} %")
    if len(con):
        lo, hi = wilson(con.lesion_class.isin(LESION_CLASSES).sum(), len(con))
        print(f"  biological negative       n {len(con):4d}  "
              f"any lesion {100 * g3:5.1f} %  "
              f"Wilson [{100 * lo:.1f}, {100 * hi:.1f}]")

    # ---- the comparison the classification is for, per the plan's section 8:
    # summarised PER GENOME and PER GENUS, never pooled over loci.
    per_genome, gkeys = [], []
    for acc, g in out.groupby("accession"):
        d = g[g.kind == "discordant"]
        c = g[g.kind == "concordant_intact"]
        if len(d) < 3 or len(c) < 3:
            continue
        fd = float((~d.lesion_class.isin(["intact"] + list(undecided))).mean())
        fc = float((~c.lesion_class.isin(["intact"] + list(undecided))).mean())
        per_genome.append({"accession": acc, "genus": g.genus.iloc[0],
                           "n_discordant": len(d), "n_intact_class": len(c),
                           "lesion_rate_discordant": fd,
                           "lesion_rate_concordant": fc,
                           "difference": fd - fc})
    pg = pd.DataFrame(per_genome)
    contrast = {}
    if len(pg) >= 5:
        pg.to_csv(OUT / "lesion_rate_by_genome.tsv", sep="\t", index=False)
        d = pg.difference.to_numpy(float)
        nz = d[d != 0]
        gb = collections.defaultdict(list)
        for _, r in pg.iterrows():
            gb[r.genus].append(r.difference)
        keys = list(gb)
        g = np.random.default_rng(20261108)
        meds = [float(np.median([x for k in g.choice(keys, len(keys),
                                                     replace=True)
                                 for x in gb[k]])) for _ in range(2000)]
        contrast = {
            "n_genomes": int(len(pg)), "n_genera": int(len(keys)),
            "median_lesion_rate_discordant":
                float(pg.lesion_rate_discordant.median()),
            "median_lesion_rate_concordant":
                float(pg.lesion_rate_concordant.median()),
            "median_difference": float(np.median(d)),
            "n_positive": int((d > 0).sum()),
            "sign_test_p_two_sided": float(stats.binomtest(
                int((nz > 0).sum()), len(nz), 0.5).pvalue),
            "genus_cluster_ci95": [float(np.percentile(meds, 2.5)),
                                   float(np.percentile(meds, 97.5))]}
        print("\nLESION RATE, discordant vs concordant-intact "
              "(per genome, genus-clustered interval)")
        print(f"  discordant  median {100 * contrast['median_lesion_rate_discordant']:.1f} %"
              f"   concordant-intact  median "
              f"{100 * contrast['median_lesion_rate_concordant']:.1f} %")
        print(f"  difference  median {100 * contrast['median_difference']:+.1f} pp"
              f"   {contrast['n_positive']}/{contrast['n_genomes']} genomes"
              f"   sign p {contrast['sign_test_p_two_sided']:.3g}"
              f"   genus CI "
              f"[{100 * contrast['genus_cluster_ci95'][0]:+.1f}, "
              f"{100 * contrast['genus_cluster_ci95'][1]:+.1f}] pp")

    gres = {}
    for name, thr in GATES.items():
        val = {"G1_frameshift_recall": g1, "G2_false_frameshift": g2,
               "G3_false_lesion": g3, "G4_decided": g4,
               "G5_internal_stop_recall": g5,
               "G6_synthetic_positive": g6}[name]
        lower_is_better = name in ("G2_false_frameshift", "G3_false_lesion")
        ok = (val is not None and
              (val <= thr if lower_is_better else val >= thr))
        gres[name] = {"value": val, "threshold": thr,
                      "direction": "<=" if lower_is_better else ">=",
                      "pass": bool(ok)}
    print("\nGATES (declared in LESION_HARMONISATION_PLAN.md section 10)")
    for k, v in gres.items():
        sv = "n/a" if v["value"] is None else f"{v['value']:.4f}"
        print(f"  {k:26s} {sv:>7} {v['direction']} {v['threshold']}  "
              f"-> {'PASS' if v['pass'] else 'FAIL'}")
    # The plan draws two DIFFERENT lines and a first version of this code
    # merged them, which would have thrown away a usable method. Section 11
    # makes the method unusable only if G1 or G3 fails. Section 12 requires
    # G1, G3, G5 and G6 for EXPANSION. G5 can fail while the method remains
    # usable at pilot scale, with the regression stated.
    usable_keys = ("G1_frameshift_recall", "G3_false_lesion")
    expand_keys = ("G1_frameshift_recall", "G3_false_lesion",
                   "G5_internal_stop_recall", "G6_synthetic_positive")
    usable = all(gres[k]["pass"] for k in usable_keys)
    expand = all(gres[k]["pass"] for k in expand_keys)
    print(f"\n  section 11, method usable (needs G1 and G3): "
          f"{'YES' if usable else 'NO'}")
    print(f"  section 12, expansion to 92 genomes "
          f"(needs G1, G3, G5, G6): "
          f"{'AUTHORISED' if expand else 'NOT authorised'}")
    if not usable:
        print("  per section 11, caller disagreement must NOT be decomposed "
              "by lesion class on these results")
    elif not expand:
        failed = [k for k in expand_keys if not gres[k]["pass"]]
        print(f"  usable at pilot scale; expansion blocked by {failed}")
    nec_ok = usable

    res = {"note": "executes LESION_HARMONISATION_PLAN.md; pilot",
           "plan_sha256": got, "input_sha256": in_sha,
           "n_pairs": len(pairs), "n_loci": int(len(out)),
           "params": {"min_hsp_aa": MIN_HSP_AA, "min_margin": MIN_MARGIN,
                      "margin_frac": MARGIN_FRAC,
                      "seg_overlap_max": SEG_OVERLAP_MAX,
                      "frag_gap_min": FRAG_GAP_MIN,
                      "stop_cover_max": STOP_COVER_MAX,
                      "trunc_cover_max": TRUNC_COVER_MAX,
                      "intact_cover_min": INTACT_COVER_MIN, "seed": SEED},
           "class_distribution": tab,
           "decided_lesion_only_fraction": g4_lesion_only,
           "lesion_rate_contrast": contrast,
           "n_with_dfast_class": int(len(lab)),
           "n_dfast_frameshift": int(len(fs)),
           "n_dfast_internal_stop": int(len(ist)),
           "gates": gres,
           "method_usable_section_11": bool(usable),
           "expansion_authorised_section_12": bool(expand),
           "necessary_gates_pass": bool(nec_ok),
           "expansion_authorised": bool(expand),
           "audit": dict(audit),
           "host_free_gb_start": free0, "host_free_gb_end": host_free_gb()}
    (P / a.out_json).write_text(json.dumps(res, indent=1, default=str,
                                           sort_keys=True))
    print(f"\nwrote {a.out_json} and {OUT}/*.tsv")


if __name__ == "__main__":
    main()
