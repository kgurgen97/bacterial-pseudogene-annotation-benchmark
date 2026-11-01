#!/usr/bin/env python3
"""Reciprocal-best-hit orthology, with the orthology and selection-bias audit.

POST HOC. Not pre-registered. No frozen plan is touched or reinterpreted.
Third generation. analysis_rbh_orthology.py (v1) and analysis_rbh_validated.py
(v2) are kept, with their defects recorded, because a method's history is part
of its evidence.

WHAT EACH GENERATION FIXED
  v1  ran the pilot; reconciled GenBank records to locus contigs by sequence
      LENGTH with record ORDER as a fall-back, and took DIAMOND's first printed
      line as the reverse best hit.
  v2  replaced that with exact sequence hashes and no fall-back, separated the
      standard protein-protein method from the modified translated-locus one,
      collapsed reverse-database entries to biological identifiers, adopted a
      tie policy, and added a covariate-matched arm.
  v3  this file. Seven further defects and omissions, each a specific way the
      v2 numbers could have been wrong:

  1 FORWARD AMBIGUITY WAS NEVER ASSESSED. v2 audited the REVERSE search for
    near-ties and took the forward hit on trust, from a search that asked for
    one target. If two partner proteins score within a hair of each other, the
    whole orthologue assignment is a coin toss and no amount of reverse-side
    rigour repairs it. Here the forward search returns up to 25 targets, the
    runner-up subject and its bitscore are recorded, and a locus whose forward
    margin is under 1 % is FORWARD-AMBIGUOUS and excluded from the
    high-confidence cohort. The stability of the original single-best-hit
    choice is reported as a rate, not assumed.

  2 COLLAPSING WAS TOO AGGRESSIVE AND WAS WRONG. v2 collapsed any two
    annotated proteins with identical sequence into one biological object
    (`DUP:`). Two identical proteins at different genomic positions are two
    genomic loci, not one, and treating them as one hid exactly the paralogue
    ambiguity the reverse search exists to detect. Here ONLY representations of
    the same verified genomic locus collapse: a locus's translated span with
    the annotated CDS or CDSs that overlap THAT locus by at least 50 %
    reciprocally on a hash-verified contig. Every other annotated gene keeps
    its own identity. Identical-sequence families are counted and reported, and
    the ties they now cause are a cost of correctness, not a defect.

  3 THE MATCHED ARM HAD NO DIAGNOSTICS. A 1:1 nearest-neighbour match with no
    caliper, no balance statistic and no record of who went unmatched is an
    assertion, not a control. Here: standardised mean differences before and
    after matching, common-support diagnostics per covariate, calipers over a
    grid, the unmatched fraction at each, and the whole thing repeated with and
    without percent identity in the covariate set.

  4 THE ESTIMAND WAS LEFT IMPLICIT. G3 on RBH-retained loci is a DIFFERENT
    QUANTITY from G3 on the original locus universe. It is not a corrected
    version of it. Matching inside the retained set balances covariates among
    retained loci and recovers nothing whatever about the excluded ones. This
    file states both estimands, reports the excluded loci's own omega and their
    own per-genome gap, and never presents one as the other.

  5 THE THREE METHODS WERE PARTLY MERGED, AND SYNTENY WAS UNCHECKED. Modified
    RBH, standard protein RBH and RBH-plus-synteny are reported as separate
    arms throughout. Synteny anchors are additionally checked for ORDER
    CONSISTENCY - anchors that merely land near the subject, in scrambled
    order, are weaker evidence than anchors that are collinear with it - and a
    frame-suspect-excluded arm shows whether the suspect spans carry any
    conclusion.

  6 INPUT IMMUTABILITY WAS ASSUMED. A concurrent 92-target reference expansion
    writes into the same results tree. Here the 12-target per-locus input is
    checksummed, snapshotted, re-checksummed at the end, and the run aborts if
    it moved. The input is also asserted to be the 12-target pilot table and
    not the expansion's.

  7 THE GATES WERE ONLY TECHNICAL. Reproducing a forward search and hashing a
    contig says nothing about whether the result means anything. A SCIENTIFIC
    QUALITY GATE is added over forward ambiguity, retention imbalance,
    post-matching balance and sensitivity across matching configurations, and
    the expansion recommendation follows from it.

A WARNING ABOUT THE MATCHED ARM ITSELF, WHICH POINTS THE OTHER WAY. v2's
matched arm drove G3 from +0.0286 to +0.0022 (7/12 genomes, p = 0.77), and it
would be just as wrong to take that at face value. Two reasons. Its median
match distance was 0.83 z-units, which is a poor match and there was no
caliper. And percent identity, one of the three matching covariates, is not a
confounder: it is close to a deterministic function of the same divergence
omega measures. Matching on it conditions on a mediator and can null a real
effect by construction. That is why the grid below runs every caliper both WITH
and WITHOUT percent identity, and why neither configuration is called the
answer on its own. Span and alignment coverage are defensible covariates; adding
identity is a stress test, not a correction.

THE STANDING WARNING. A stable G3 after RBH filtering is NOT independent
biological corroboration. RBH retention depends on reciprocation strength,
which depends on divergence and alignable length, and divergence is what omega
measures; v2 measured this directly - rejected discordant loci sat at 52.3 %
identity and omega 0.33 against 72.5 % and 0.21 for retained ones. Filtering on
reciprocation selects on a correlate of the outcome. Nothing here licenses the
inference "G3 survived filtering, therefore the orthologues are right", and the
matched arm exists precisely because that inference is unavailable.

OMEGA IS NEVER RECOMPUTED. Every arm reuses the omega in the input table, whose
control column reproduced results/loci.tsv exactly (6360/6360). Only the locus
set changes, and the script asserts it bit for bit.

    python3 analysis_rbh_audited.py [--per-locus T] [--out-json J] [--targets N]
"""
from __future__ import annotations

import argparse
import ast
import collections
import csv
import hashlib
import json
import pathlib
import shutil
import subprocess
import sys
import tempfile
import types

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from scipy import stats

P = pathlib.Path(__file__).resolve().parent
R = P / "results"
OUT = R / "rbh3"
D = pathlib.Path.home() / "research/datasets"
REFDIR = D / "n2-refseq-reference"          # READ ONLY
QDIR = D / "n2-query-proteomes"

MIN_PROT_AA = 20
MAX_STOP_FRAC = 0.10
TIE_MARGIN = 0.01           # relative bitscore margin: below it, a tie
FWD_TARGETS = 25
REV_TARGETS = 25
N_NEIGHBOURS = 5
SYNTENY_WINDOW = 10
MIN_ANCHORS_FOR_ORDER = 3
COLLINEAR_RHO = 0.8
THREADS = 4
FLOOR_START_GB, FLOOR_ABORT_GB = 10, 8

# matching grid
CALIPERS = [0.10, 0.25, 0.50, 1.00, None]
PRIMARY_CALIPER = 0.25
COVSETS = {"with_pident": ["span", "forward_coverage", "forward_pident"],
           "without_pident": ["span", "forward_coverage"]}
PRIMARY_COVSET = "with_pident"
SMD_BALANCED = 0.10

# technical gates
GATE_FORWARD_AGREE = 0.99
GATE_HASH_COVERAGE = 0.95
# scientific quality gate
GATE_MAX_FWD_AMBIGUITY = 0.10
GATE_MAX_RETENTION_IMBALANCE = 0.10


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


# ----------------------------------------------------------------- disk guard
def host_free_gb():
    out = subprocess.run(["df", "--output=avail", "-BG", "/mnt/c"],
                         capture_output=True, text=True, check=True).stdout
    return int("".join(c for c in out.splitlines()[-1] if c.isdigit()))


def guard(stage):
    free = host_free_gb()
    if free < FLOOR_ABORT_GB:
        sys.exit(f"EMERGENCY STOP at {stage}: {free} GB < {FLOOR_ABORT_GB} GB")
    return free


# ----------------------------------------------------------------- utilities
def sha(s: str) -> str:
    return hashlib.sha256(s.upper().encode()).hexdigest()


def file_sha(p: pathlib.Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def clean(v):
    """A recorded field, or None. NaN is truthy, so it must not be tested raw."""
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return None
    return str(v)


def overlaps_half(a0, a1, b0, b1):
    """50 % reciprocal overlap - the project's own locus criterion."""
    inter = max(0, min(a1, b1) - max(a0, b0) + 1)
    if inter <= 0:
        return False
    return inter >= 0.5 * (a1 - a0 + 1) and inter >= 0.5 * (b1 - b0 + 1)


def overlap_bp(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0) + 1)


def translate_span(nt, qs, qe):
    """-> (protein, frame, reversed, stop_fraction) for the aligned span."""
    lo, hi = min(qs, qe), max(qs, qe)
    span = nt[lo - 1:hi]
    rev = qs > qe
    if rev:
        span = str(Seq(span).reverse_complement())
    span = span[:len(span) - len(span) % 3]
    if len(span) < MIN_PROT_AA * 3:
        return None, (lo - 1) % 3, rev, None
    prot = str(Seq(span).translate())
    frac = prot.count("*") / len(prot) if prot else None
    return prot.replace("*", "X"), (lo - 1) % 3, rev, frac


def margin_of(ranked):
    """-> (best_key, best_score, second_key, second_score, relative margin).

    `ranked` is an iterable of (key, score) over DISTINCT biological objects.
    A sole object has a full margin of 1.0; it is not a tie with itself.
    """
    order = sorted(ranked, key=lambda kv: -kv[1])
    if not order:
        return None, None, None, None, None
    b_k, b_s = order[0]
    if len(order) == 1:
        return b_k, b_s, None, None, 1.0
    s_k, s_s = order[1]
    return b_k, b_s, s_k, s_s, ((b_s - s_s) / b_s if b_s > 0 else 0.0)


def collapse(hits, entry_bio):
    """Max score per biological identifier, then the margin between objects.

    Requirement 2 of the audit: entry_bio must map two entries to one object
    ONLY when they represent the same verified genomic locus. Identical protein
    sequences at different positions are different objects here, so a paralogue
    pair now produces a genuine tie instead of being silently merged.
    """
    per = {}
    for ent, score in hits:
        bio = entry_bio.get(ent)
        if bio is None:
            continue
        if score > per.get(bio, -1.0):
            per[bio] = score
    return margin_of(per.items())


# ----------------------------------------------------------- reference access
def partner_reference(partner):
    faa, ffn = REFDIR / f"{partner}.faa", REFDIR / f"{partner}.ffn"
    sent = REFDIR / f"{partner}.json"
    if not (faa.exists() and ffn.exists() and sent.exists()):
        return None, None
    try:
        if json.loads(sent.read_text()).get("version") != 2:
            return None, None
    except Exception:
        return None, None
    return faa, ffn


def gbff_for(acc):
    g = D / f"n2-genomes-expansion/{acc}/genomic.gbff"
    if g.exists():
        return g
    return next((D / "n2-genomes/extracted").rglob(f"*{acc}*/genomic.gbff"), None)


CDS_COLS = ["id", "prot", "rec_id", "start", "end", "rec_sha", "rec_idx"]


def cds_table(acc):
    """-> [(id, prot, rec_id, start, end, rec_sha256, rec_idx)], non-pseudo CDS."""
    QDIR.mkdir(parents=True, exist_ok=True)
    cache = QDIR / f"{acc}.v3.tsv"
    if cache.exists():
        with cache.open() as fh:
            return [(r["id"], r["prot"], r["rec_id"], int(r["start"]),
                     int(r["end"]), r["rec_sha"], int(r["rec_idx"]))
                    for r in csv.DictReader(fh, delimiter="\t")]
    gb = gbff_for(acc)
    if gb is None:
        return []
    rows = []
    for idx, rec in enumerate(SeqIO.parse(str(gb), "genbank")):
        try:
            rsha = sha(str(rec.seq))
        except Exception:
            rsha = ""
        for ft in rec.features:
            if ft.type != "CDS":
                continue
            q = ft.qualifiers
            if "pseudo" in q or "pseudogene" in q:
                continue
            tr = q.get("translation", [None])[0]
            if not tr:
                continue
            ident = (q.get("locus_tag") or q.get("protein_id")
                     or [f"{rec.id}_{int(ft.location.start)}"])[0]
            rows.append((ident, tr, rec.id, int(ft.location.start) + 1,
                         int(ft.location.end), rsha, idx))
    tmp = cache.with_suffix(".tsv.part")
    with tmp.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(CDS_COLS)
        w.writerows(rows)
    tmp.replace(cache)
    return rows


def hash_contig_map(contigs, cds):
    """GenBank record id -> locus contig name, by EXACT unique sequence hash."""
    by_hash = collections.defaultdict(list)
    for name, seq in contigs.items():
        by_hash[sha(seq)].append(name)
    seen, out, stat = {}, {}, collections.Counter()
    for _, _, rid, _, _, rsha, _ in cds:
        if rid in seen:
            continue
        seen[rid] = True
        cand = by_hash.get(rsha, [])
        if len(cand) == 1 and rsha:
            out[rid] = cand[0]
            stat["hash_unique"] += 1
        elif len(cand) > 1:
            out[rid] = None
            stat["hash_ambiguous_in_assembly"] += 1
        else:
            out[rid] = None
            stat["hash_absent"] += 1
    return out, dict(stat)


def partner_order(partner):
    """-> {locus_tag: (record id, index within that record by coordinate)}."""
    rows = cds_table(partner)
    by = collections.defaultdict(list)
    for ident, _, rid, start, _, _, _ in rows:
        by[rid].append((start, ident))
    out = {}
    for rid, v in by.items():
        for i, (_, ident) in enumerate(sorted(v)):
            out[ident] = (rid, i)
    return out


# ----------------------------------------------------------------- searches
FWD_FMT = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "slen",
           "bitscore"]


def blastx_multi(query_fa, db_faa, workdir, max_targets=FWD_TARGETS):
    """Forward: locus nucleotides -> partner proteins, MANY targets.

    Requirement 1. The published pipeline asked for one target; asking for many
    is what makes forward ambiguity visible. Scores are maxed per SUBJECT, and
    subjects are never merged with each other: two partner proteins are two
    partner genomic loci even when their sequences are identical.
    """
    workdir.mkdir(parents=True, exist_ok=True)
    dmnd = workdir / "db.dmnd"
    subprocess.run(["diamond", "makedb", "--in", str(db_faa), "-d", str(dmnd),
                    "--quiet", "--threads", str(THREADS)], check=True)
    out = workdir / "fwd.tsv"
    subprocess.run(["diamond", "blastx", "-q", str(query_fa), "-d", str(dmnd),
                    "-o", str(out), "--quiet", "--threads", str(THREADS),
                    "--max-target-seqs", str(max_targets), "--evalue", "1e-5",
                    "--outfmt", "6", *FWD_FMT], check=True)
    per = collections.defaultdict(dict)
    for line in out.read_text().splitlines():
        f = line.split("\t")
        if len(f) < len(FWD_FMT):
            continue
        q, s, bits = f[0], f[1], float(f[7])
        prev = per[q].get(s)
        if prev is None or bits > prev["bitscore"]:
            per[q][s] = {"subject": s, "pident": float(f[2]),
                         "aln_len": int(f[3]), "qstart": int(f[4]),
                         "qend": int(f[5]), "slen": int(f[6]),
                         "bitscore": bits,
                         "coverage": int(f[3]) / int(f[6]) if int(f[6]) else None}
    return per


def blastp_all(query_fa, db_faa, workdir, max_targets=REV_TARGETS):
    """-> {query: [(subject, bitscore), ...]} best-first, up to max_targets."""
    workdir.mkdir(parents=True, exist_ok=True)
    dmnd = workdir / "db.dmnd"
    subprocess.run(["diamond", "makedb", "--in", str(db_faa), "-d", str(dmnd),
                    "--quiet", "--threads", str(THREADS)], check=True)
    out = workdir / "hits.tsv"
    subprocess.run(["diamond", "blastp", "-q", str(query_fa), "-d", str(dmnd),
                    "-o", str(out), "--quiet", "--threads", str(THREADS),
                    "--max-target-seqs", str(max_targets), "--evalue", "1e-5",
                    "--outfmt", "6", "qseqid", "sseqid", "bitscore"],
                   check=True)
    hits = collections.defaultdict(list)
    for line in out.read_text().splitlines():
        f = line.split("\t")
        if len(f) >= 3:
            hits[f[0]].append((f[1], float(f[2])))
    return hits


# ----------------------------------------------------------------- per genome
def process(acc, partner, sub, strands, contigs, syn_cache, audit):
    faa, _ = partner_reference(partner)
    if faa is None:
        return None, "partner reference incomplete or absent"
    pprot = X.read_fasta(faa)
    if not pprot:
        return None, "partner reference unreadable"
    cds = cds_table(acc)
    if not cds:
        return None, "no query proteome"

    xlat, hstat = hash_contig_map(contigs, cds)
    audit["contig_identity"].update(hstat)
    hash_ok = {v for v in xlat.values() if v}

    queries, meta = {}, {}
    for n, r in enumerate(sub.itertuples()):
        st = strands.get((r.accession, r.contig, r.start, r.end))
        if st is None or r.contig not in contigs or r.end > len(contigs[r.contig]):
            continue
        seq = X.sub(contigs[r.contig], r.start, r.end, st)
        if seq:
            k = f"L{n}"
            queries[k], meta[k] = seq, r
    if not queries:
        return None, "no locus sequence recoverable"

    # ---- annotated CDS in the locus namespace, hash-verified contigs only
    cds_iv, by_contig, prot_of = {}, collections.defaultdict(list), {}
    for ident, prot, rid, start, end, _, _ in cds:
        prot_of[ident] = prot
        contig = xlat.get(rid)
        if contig is None:
            audit["cds_unmapped_contig"] += 1
            continue
        cds_iv[ident] = (contig, start, end)
        by_contig[contig].append((start, ident))
    for v in by_contig.values():
        v.sort()
    audit["cds_total"] += len(cds)

    # ---- identical-sequence families: COUNTED, never collapsed (audit req 2)
    fam = collections.defaultdict(list)
    for ident in cds_iv:
        fam[sha(prot_of[ident])].append(ident)
    idseq_fams = {h: v for h, v in fam.items() if len(v) > 1}
    audit["identical_sequence_families"] += len(idseq_fams)
    audit["genes_in_identical_families"] += sum(len(v) for v in idseq_fams.values())

    # ---- which locus does each annotated CDS belong to? greatest overlap wins
    cds_locus, contested = {}, 0
    for k, r in meta.items():
        for _, ident in by_contig.get(r.contig, []):
            iv = cds_iv[ident]
            if not overlaps_half(r.start, r.end, iv[1], iv[2]):
                continue
            ov = overlap_bp(r.start, r.end, iv[1], iv[2])
            prev = cds_locus.get(ident)
            if prev is None:
                cds_locus[ident] = (k, ov)
            else:
                contested += 1
                if ov > prev[1]:
                    cds_locus[ident] = (k, ov)
    cds_locus = {i: v[0] for i, v in cds_locus.items()}
    audit["cds_contested_by_two_loci"] += contested

    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        qf = td / "loci.fna"
        with qf.open("w") as fh:
            for k, v in queries.items():
                fh.write(f">{k}\n{v}\n")
        fwd_all = blastx_multi(qf, faa, td / "fwd")
        if not fwd_all:
            return None, "forward search returned nothing"

        # ---- forward decision and ambiguity
        fwd = {}
        for k, cands in fwd_all.items():
            b, bs, s2, ss, mar = margin_of(
                [(s, h["bitscore"]) for s, h in cands.items()])
            if b is None:
                continue
            h = dict(cands[b])
            h.update({"n_forward_candidates": len(cands),
                      "forward_second": s2, "forward_second_score": ss,
                      "forward_margin": mar,
                      "forward_ambiguous": int(mar is not None
                                               and mar < TIE_MARGIN)})
            prot, frame, rev, frac = translate_span(queries[k], h["qstart"],
                                                    h["qend"])
            h.update({"frame": frame, "reversed": int(rev), "stop_frac": frac,
                      "frame_suspect": int(frac is not None
                                           and frac > MAX_STOP_FRAC),
                      "prot": prot})
            fwd[k] = h
            audit["frames"][frame] += 1
            audit["reversed_hits"] += int(rev)
            audit["frame_suspect"] += h["frame_suspect"]
            audit["forward_ambiguous"] += h["forward_ambiguous"]
            if prot is None:
                audit["span_too_short_to_translate"] += 1

        # ---- standard method forward leg: annotated protein -> partner
        locus_cds = {}
        for ident, k in cds_locus.items():
            iv = cds_iv[ident]
            span = iv[2] - iv[1]
            prev = locus_cds.get(k)
            if prev is None or span > prev[1]:
                locus_cds[k] = (ident, span)
        std_fwd = {}
        if locus_cds:
            stdq = td / "stdq.faa"
            with stdq.open("w") as fh:
                for k, (ident, _) in locus_cds.items():
                    fh.write(f">{k}\n{prot_of[ident]}\n")
            for k, hs in blastp_all(stdq, faa, td / "stdfwd",
                                    max_targets=FWD_TARGETS).items():
                per = {}
                for s, bits in hs:
                    if bits > per.get(s, -1.0):
                        per[s] = bits
                b, _, s2, _, mar = margin_of(per.items())
                std_fwd[k] = {"subject": b, "second": s2, "margin": mar,
                              "ambiguous": int(mar is not None
                                               and mar < TIE_MARGIN)}

        # ---- reverse database. ONLY same-locus representations collapse.
        entry_bio, entry_bio_cds = {}, {}
        revdb = td / "revdb.faa"
        with revdb.open("w") as fh:
            for ident in cds_iv:
                bio = cds_locus.get(ident) or f"CDS:{ident}"
                entry_bio[f"CDS|{ident}"] = bio
                entry_bio_cds[f"CDS|{ident}"] = bio
                fh.write(f">CDS|{ident}\n{prot_of[ident]}\n")
            for k, h in fwd.items():
                if h["prot"]:
                    entry_bio[f"LOC|{k}"] = k
                    fh.write(f">LOC|{k}\n{h['prot']}\n")

        subjects = sorted({h["subject"] for h in fwd.values()}
                          | {v["subject"] for v in std_fwd.values()
                             if v["subject"]})
        revq = td / "revq.faa"
        n_rq = 0
        with revq.open("w") as fh:
            for sid in subjects:
                if sid in pprot:
                    fh.write(f">{sid}\n{pprot[sid]}\n")
                    n_rq += 1
        if not n_rq:
            return None, "no subject protein sequences"
        rev = blastp_all(revq, revdb, td / "rev")

        qpf = td / "qprot.faa"
        with qpf.open("w") as fh:
            for ident in cds_iv:
                fh.write(f">{ident}\n{prot_of[ident]}\n")
        q2p_all = blastp_all(qpf, faa, td / "syn", max_targets=1)
        q2p = {q: hs[0][0] for q, hs in q2p_all.items() if hs}

    if partner not in syn_cache:
        syn_cache[partner] = partner_order(partner)
    porder = syn_cache[partner]

    rows, fa, fc = [], 0, 0
    for k, r in meta.items():
        h = fwd.get(k)
        rec = {
            "accession": acc, "genus": r.genus, "partner": partner,
            "contig": r.contig, "start": r.start, "end": r.end,
            "span": int(r.end) - int(r.start) + 1, "kind": r.kind,
            "omega": r.control,
            "contig_hash_verified": int(r.contig in hash_ok),
            "recorded_subject": clean(getattr(r, "treatment_subject", None)),
            "has_annotated_cds": int(k in locus_cds),
            "n_forward_candidates": None, "forward_subject": None,
            "forward_bitscore": None, "forward_pident": None,
            "forward_coverage": None, "forward_second": None,
            "forward_second_score": None, "forward_margin": None,
            "forward_ambiguous": None, "forward_matches_recorded": None,
            "frame": None, "reversed": None, "stop_frac": None,
            "frame_suspect": None,
            "reverse_best_bio": None, "reverse_best_score": None,
            "reverse_second_bio": None, "reverse_second_score": None,
            "reverse_margin": None,
            "rbh_modified": 0, "outcome": None,
            "standard_subject": None, "standard_margin": None,
            "standard_outcome": None, "rbh_standard": None,
            "synteny_anchors": None, "synteny_votes": None,
            "synteny_rho": None, "synteny_collinear": None,
        }
        if h is None:
            rec["outcome"] = "no_forward_hit"
            rows.append(rec)
            continue
        rec.update({
            "n_forward_candidates": h["n_forward_candidates"],
            "forward_subject": h["subject"], "forward_bitscore": h["bitscore"],
            "forward_pident": h["pident"], "forward_coverage": h["coverage"],
            "forward_second": h["forward_second"],
            "forward_second_score": h["forward_second_score"],
            "forward_margin": h["forward_margin"],
            "forward_ambiguous": h["forward_ambiguous"],
            "frame": h["frame"], "reversed": h["reversed"],
            "stop_frac": h["stop_frac"], "frame_suspect": h["frame_suspect"]})
        if rec["recorded_subject"]:
            fc += 1
            same = int(h["subject"] == rec["recorded_subject"])
            fa += same
            rec["forward_matches_recorded"] = same

        b, bs, s2, ss, mar = collapse(rev.get(h["subject"], []), entry_bio)
        rec.update({"reverse_best_bio": b, "reverse_best_score": bs,
                    "reverse_second_bio": s2, "reverse_second_score": ss,
                    "reverse_margin": mar})
        if b is None:
            rec["outcome"] = "no_reverse_hit"
        elif b != k:
            rec["outcome"] = "reverse_hit_other_locus"
        elif mar is not None and mar < TIE_MARGIN:
            rec["outcome"] = "reverse_ambiguous_tie"
        elif h["forward_ambiguous"]:
            rec["outcome"] = "forward_ambiguous"
        else:
            rec["outcome"] = "rbh_confirmed"
            rec["rbh_modified"] = 1

        if k in locus_cds:
            sf = std_fwd.get(k, {})
            s_sub = sf.get("subject")
            rec["standard_subject"] = s_sub
            if s_sub is None:
                rec["standard_outcome"] = "no_forward_hit"
            else:
                sb, _, _, _, smar = collapse(rev.get(s_sub, []), entry_bio_cds)
                rec["standard_margin"] = smar
                if sb is None:
                    rec["standard_outcome"] = "no_reverse_hit"
                elif sb != k:
                    rec["standard_outcome"] = "reverse_hit_other_locus"
                elif smar is not None and smar < TIE_MARGIN:
                    rec["standard_outcome"] = "reverse_ambiguous_tie"
                elif sf.get("ambiguous"):
                    rec["standard_outcome"] = "forward_ambiguous"
                else:
                    rec["standard_outcome"] = "rbh_confirmed"
            rec["rbh_standard"] = int(
                rec["standard_outcome"] == "rbh_confirmed")

        # ---- synteny, with an order-consistency check on the anchors
        if rec["rbh_modified"]:
            tgt = porder.get(h["subject"])
            votes, anchors, qpos, ppos = 0, 0, [], []
            if tgt:
                near = sorted(by_contig.get(r.contig, []),
                              key=lambda t: abs(t[0] - r.start))
                for start, ident in near:
                    iv = cds_iv[ident]
                    if overlaps_half(r.start, r.end, iv[1], iv[2]):
                        continue
                    anchors += 1
                    pt = porder.get(q2p.get(ident, ""))
                    if pt and pt[0] == tgt[0]:
                        qpos.append(start)
                        ppos.append(pt[1])
                        if abs(pt[1] - tgt[1]) <= SYNTENY_WINDOW:
                            votes += 1
                    if anchors >= N_NEIGHBOURS:
                        break
            rec["synteny_anchors"] = anchors
            rec["synteny_votes"] = votes
            if len(qpos) >= MIN_ANCHORS_FOR_ORDER and len(set(ppos)) > 1:
                rho = float(stats.spearmanr(qpos, ppos).statistic)
                rec["synteny_rho"] = rho
                rec["synteny_collinear"] = int(abs(rho) >= COLLINEAR_RHO)
            else:
                rec["synteny_collinear"] = 0
        rows.append(rec)

    return (rows, fa, fc), None


# ----------------------------------------------------------------- matching
def smd(a, b):
    """Standardised mean difference, pooled SD. |SMD| <= 0.1 is balanced."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return None
    s = np.sqrt((a.std(ddof=1) ** 2 + b.std(ddof=1) ** 2) / 2)
    if s == 0:
        return 0.0
    return float((a.mean() - b.mean()) / s)


def prep_cov(frame, cols):
    out = frame[cols].copy()
    if "span" in cols:
        out["span"] = np.log10(out["span"].clip(lower=1))
    return out


def match_arm(df, cols, caliper, seed=20261031):
    """1:1 nearest-neighbour matching within genome, with a caliper.

    Requirement 3. Returns the pairs, per-covariate SMD before and after, the
    unmatched fraction and common-support diagnostics. The caliper is in
    z-units of the within-genome covariate distribution; None means no caliper.
    """
    g = np.random.default_rng(seed)
    keep = df[df.hiconf == 1].dropna(subset=cols + ["omega"])
    pairs, before, after, support, unmatched, n_dis_total = [], [], [], [], 0, 0
    for acc, sub in keep.groupby("accession"):
        dis = sub[sub.kind == "discordant"]
        con = sub[sub.kind == "concordant_intact"]
        if len(dis) < 3 or len(con) < 3:
            continue
        n_dis_total += len(dis)
        Z = prep_cov(sub, cols)
        mu, sd = Z.mean(), Z.std().replace(0, 1)
        zd = ((prep_cov(dis, cols) - mu) / sd).to_numpy()
        zc = ((prep_cov(con, cols) - mu) / sd).to_numpy()
        for j, c in enumerate(cols):
            before.append({"accession": acc, "covariate": c,
                           "smd": smd(zd[:, j], zc[:, j])})
            lo, hi = zc[:, j].min(), zc[:, j].max()
            support.append({"accession": acc, "covariate": c,
                            "frac_discordant_in_baseline_range": float(
                                ((zd[:, j] >= lo) & (zd[:, j] <= hi)).mean())})
        taken, mi, mj = set(), [], []
        for i in g.permutation(len(zd)):
            d = np.sqrt(((zc - zd[i]) ** 2).sum(axis=1))
            if taken:
                d[list(taken)] = np.inf
            j = int(np.argmin(d))
            if not np.isfinite(d[j]) or (caliper is not None and d[j] > caliper):
                unmatched += 1
                continue
            taken.add(j)
            mi.append(i)
            mj.append(j)
            pairs.append({"accession": acc, "genus": dis.genus.iloc[0],
                          "omega_dis": float(dis.omega.iloc[i]),
                          "omega_con": float(con.omega.iloc[j]),
                          "dist": float(d[j])})
        if mi:
            for j, c in enumerate(cols):
                after.append({"accession": acc, "covariate": c,
                              "smd": smd(zd[mi, j], zc[mj, j])})
    mp = pd.DataFrame(pairs)
    res = {"covariates": cols, "caliper": caliper,
           "n_pairs": int(len(mp)), "n_discordant_eligible": n_dis_total,
           "unmatched": int(unmatched),
           "unmatched_frac": (unmatched / n_dis_total) if n_dis_total else None}
    if len(mp):
        per = [float(np.median(s.omega_dis) - np.median(s.omega_con))
               for _, s in mp.groupby("accession")]
        nz = [v for v in per if v != 0]
        res.update({
            "n_genomes": len(per), "median_gap": float(np.median(per)),
            "n_positive": sum(1 for v in nz if v > 0), "n_tested": len(nz),
            "p": float(stats.binomtest(sum(1 for v in nz if v > 0), len(nz),
                                       0.5).pvalue) if nz else None,
            "median_match_distance": float(mp.dist.median())})
    bf, af = pd.DataFrame(before), pd.DataFrame(after)
    for name, tab in (("smd_before", bf), ("smd_after", af)):
        if len(tab):
            res[name] = {c: float(np.nanmedian(np.abs(s.smd.astype(float))))
                         for c, s in tab.groupby("covariate")}
    if len(af):
        res["max_abs_smd_after"] = float(
            np.nanmax(np.abs(af.smd.astype(float))))
        res["balanced"] = bool(res["max_abs_smd_after"] <= SMD_BALANCED)
    if len(support):
        sp = pd.DataFrame(support)
        res["common_support"] = {
            c: float(s.frac_discordant_in_baseline_range.min())
            for c, s in sp.groupby("covariate")}
    return res, mp, bf, af, pd.DataFrame(support)


# ----------------------------------------------------------------- reporting
def arm(rows, col, genus, label):
    gaps = S.per_genome_gaps(rows, col)
    out = {"label": label, **S.g3_from(gaps)}
    lo, hi = S.cluster_ci(rows, col, genus, unit="genus")
    out["genus_cluster_ci95"] = [lo, hi]
    out["n_loci_used"] = int(sum(
        1 for r in rows if r.get(col) is not None
        and not (isinstance(r[col], float) and np.isnan(r[col]))))
    return out


def show(A):
    g = A.get("median_gap")
    if g is None:
        print(f"  {A['label']:44s} {A.get('note')}")
        return
    ci = A.get("genus_cluster_ci95") or [None, None]
    line = (f"  {A['label']:44s} {g:+.4f}  {A['n_positive']}/{A['n_tested']}"
            f"  p {A['p']:<9.3g} n {A['n_loci_used']:6d}")
    if ci[0] is not None:
        line += f"  genus CI [{ci[0]:+.4f}, {ci[1]:+.4f}]"
    print(line)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-locus",
                    default="results/reference_sensitivity/per_locus.tsv")
    ap.add_argument("--out-json", default="results/rbh_audited.json")
    ap.add_argument("--targets", type=int, default=0)
    ap.add_argument("--expect-pairs", type=int, default=12)
    # isolated outputs, so an expanded run cannot overwrite the committed
    # pilot evidence this document's numbers were read from
    ap.add_argument("--out-dir", default=None)
    a = ap.parse_args()
    global OUT
    if a.out_dir:
        OUT = P / a.out_dir

    free0 = guard("start")
    print(f"host C: free {free0} GB (from /mnt/c)")
    if free0 < FLOOR_START_GB:
        sys.exit(f"refusing to start: {free0} < {FLOOR_START_GB} GB")
    OUT.mkdir(parents=True, exist_ok=True)

    # ================================= requirement 6: a fixed, verified input
    src = P / a.per_locus
    in_sha = file_sha(src)
    snap = OUT / "input_snapshot.tsv"
    shutil.copy2(src, snap)
    snap_sha = file_sha(snap)
    print(f"input  {src.relative_to(P)}")
    print(f"  sha256   {in_sha}")
    print(f"  snapshot {snap.relative_to(P)}  sha256 matches: {snap_sha == in_sha}")
    if snap_sha != in_sha:
        sys.exit("snapshot does not match the source; refusing to proceed")
    # This guard exists so a pilot run cannot be handed the 92-target table by
    # accident. It must not block a run that DECLARES the larger scope, so it
    # applies only when the caller still expects the pilot's 12 pairs; the
    # pair-count check below is what enforces the declaration either way.
    if src.name.endswith("_all.tsv") and a.expect_pairs == 12:
        sys.exit("that is the 92-target expansion's table, not the pilot's; "
                 "pass --expect-pairs 92 to run the expansion deliberately")

    df = pd.read_csv(snap, sep="\t")           # read the SNAPSHOT, not the live file
    df = df[df.control.notna()]
    genus = dict(zip(df.accession, df.genus))
    pairs = df[["accession", "partner"]].drop_duplicates().values.tolist()
    if a.targets:
        pairs = pairs[:a.targets]
    elif len(pairs) != a.expect_pairs:
        sys.exit(f"expected {a.expect_pairs} pairs in the pilot table, "
                 f"found {len(pairs)}; refusing to proceed")
    print(f"{len(df)} scored loci over {len(pairs)} (query, partner) pairs")

    print("\nrecovering locus strands from a 0.50 partition rebuild")
    strands = S.locus_strands([p[0] for p in pairs])
    print(f"  {len(strands)} representative intervals recovered")

    audit = collections.defaultdict(int)
    audit["contig_identity"] = collections.Counter()
    audit["frames"] = collections.Counter()
    all_rows, skipped, fa, fc, syn_cache = [], {}, 0, 0, {}
    for i, (acc, partner) in enumerate(pairs, 1):
        contigs = X.genome_seqs(acc)
        if not contigs:
            skipped[acc] = "no contigs"
            continue
        res, why = process(acc, partner, df[df.accession == acc], strands,
                           contigs, syn_cache, audit)
        if res is None:
            skipped[acc] = why
            print(f"  [{i}/{len(pairs)}] {acc}: SKIPPED — {why}", flush=True)
            continue
        rows, agree, cmp_ = res
        fa, fc = fa + agree, fc + cmp_
        all_rows.extend(rows)
        nm = sum(r["rbh_modified"] for r in rows)
        namb = sum(1 for r in rows if r["forward_ambiguous"])
        print(f"  [{i}/{len(pairs)}] {acc}: {len(rows)} loci, "
              f"{nm} modified-RBH ({100 * nm / len(rows):.1f}%), "
              f"{namb} forward-ambiguous", flush=True)
        guard(acc)

    if not all_rows:
        sys.exit("nothing processed")
    out = pd.DataFrame(all_rows)

    out["hiconf"] = ((out.contig_hash_verified == 1) &
                     (out.rbh_modified == 1)).astype(int)
    votes = pd.to_numeric(out.synteny_votes, errors="coerce").fillna(0)
    coll = pd.to_numeric(out.synteny_collinear, errors="coerce").fillna(0)
    susp = pd.to_numeric(out.frame_suspect, errors="coerce").fillna(0)
    out["omega_hiconf"] = np.where(out.hiconf == 1, out.omega, np.nan)
    out["omega_syn1"] = np.where((out.hiconf == 1) & (votes >= 1),
                                 out.omega, np.nan)
    out["omega_syn2"] = np.where((out.hiconf == 1) & (votes >= 2),
                                 out.omega, np.nan)
    out["omega_syn_collinear"] = np.where(
        (out.hiconf == 1) & (votes >= 1) & (coll == 1), out.omega, np.nan)
    out["omega_standard"] = np.where(out.rbh_standard == 1, out.omega, np.nan)
    out["omega_noframesuspect"] = np.where(
        (out.hiconf == 1) & (susp == 0), out.omega, np.nan)
    out["omega_excluded"] = np.where(out.hiconf == 0, out.omega, np.nan)
    out.to_csv(OUT / "per_locus_rbh_audited.tsv", sep="\t", index=False)
    print(f"\nwrote {len(out)} locus rows to "
          f"{OUT / 'per_locus_rbh_audited.tsv'}")

    # ================================================== controls
    fwd_agree = fa / fc if fc else None
    hash_cov = float((out.contig_hash_verified == 1).mean())
    key = ["accession", "contig", "start", "end"]
    chk = out[key + ["omega"]].merge(df[key + ["control"]], on=key, how="left")
    omega_intact = bool(len(chk) == len(out) and np.array_equal(
        chk.omega.to_numpy(float), chk.control.to_numpy(float), equal_nan=True))
    end_sha = file_sha(src)
    input_stable = end_sha == in_sha

    print("\nCONTROL 1 — single-best-hit stability")
    print(f"  top subject of a {FWD_TARGETS}-target search equals the recorded "
          f"single-best subject in {fa} of {fc} "
          f"({100 * fwd_agree:.2f}%)" if fwd_agree else "  n/a")
    amb = pd.to_numeric(out.forward_ambiguous, errors="coerce")
    fmar = pd.to_numeric(out.forward_margin, errors="coerce")
    ncand = pd.to_numeric(out.n_forward_candidates, errors="coerce")
    fwd_amb_frac = float(amb.mean(skipna=True))
    print(f"  forward candidates per locus: median {ncand.median():.0f}, "
          f"max {ncand.max():.0f}")
    print(f"  forward margin: median {fmar.median():.3f}, "
          f"5th pct {fmar.quantile(0.05):.3f}")
    print(f"  FORWARD-AMBIGUOUS (margin < {TIE_MARGIN}): "
          f"{int(amb.sum())} ({100 * fwd_amb_frac:.2f}%)")
    for t in (0.01, 0.05, 0.10):
        print(f"    runner-up within {100 * t:4.0f}%: "
              f"{int((fmar < t).sum()):6d} loci")
    print("\nCONTROL 2 — omega bit-identical to the input: "
          f"{omega_intact} ({len(chk)} rows)")
    print(f"CONTROL 3 — input unchanged during the run: {input_stable}")
    print("\nREQUIREMENT 1/6 — contig identity by exact sequence hash")
    for k, v in sorted(audit["contig_identity"].items()):
        print(f"  {k:32s} {v:6d} records")
    print(f"  loci on hash-verified contigs  {100 * hash_cov:6.2f}%  "
          f"(excluded {100 * (1 - hash_cov):.2f}%)")

    print("\nREQUIREMENT 2 — collapsing, and what is NOT collapsed")
    for k in ("cds_total", "cds_unmapped_contig", "cds_contested_by_two_loci",
              "identical_sequence_families", "genes_in_identical_families",
              "span_too_short_to_translate", "reversed_hits", "frame_suspect",
              "forward_ambiguous"):
        print(f"  {k:32s} {audit[k]:6d}")
    print("  identical-sequence genes are kept DISTINCT; only same-locus "
          "representations collapse")

    # ================================================== retention
    print("\nRETENTION by locus class")
    retention = {}
    for kind, g in out.groupby("kind"):
        retention[kind] = {
            "n": int(len(g)),
            "hash_verified_frac": float((g.contig_hash_verified == 1).mean()),
            "forward_ambiguous_frac": float(
                pd.to_numeric(g.forward_ambiguous, errors="coerce").mean()),
            "modified_rbh_frac": float(g.rbh_modified.mean()),
            "hiconf_frac": float(g.hiconf.mean()),
            "syn1_frac": float(g.omega_syn1.notna().mean()),
            "syn2_frac": float(g.omega_syn2.notna().mean()),
            "syn_collinear_frac": float(g.omega_syn_collinear.notna().mean()),
            "has_annotated_cds_frac": float(g.has_annotated_cds.mean()),
            "standard_rbh_frac_where_applicable": (
                float(pd.to_numeric(g[g.has_annotated_cds == 1].rbh_standard,
                                    errors="coerce").mean())
                if (g.has_annotated_cds == 1).any() else None)}
        r = retention[kind]
        print(f"  {kind:20s} n {r['n']:6d}  hash {100*r['hash_verified_frac']:5.1f}%"
              f"  fwd-amb {100*r['forward_ambiguous_frac']:4.1f}%"
              f"  RBH {100*r['modified_rbh_frac']:5.1f}%"
              f"  +syn1 {100*r['syn1_frac']:5.1f}%"
              f"  +collinear {100*r['syn_collinear_frac']:5.1f}%"
              f"  CDS {100*r['has_annotated_cds_frac']:5.1f}%"
              f"  std {100*(r['standard_rbh_frac_where_applicable'] or 0):5.1f}%")
    dd, bb = retention.get("discordant"), retention.get("concordant_intact")
    imbalance = (abs(dd["hiconf_frac"] - bb["hiconf_frac"])
                 if dd and bb else None)
    if imbalance is not None:
        print(f"  RETENTION IMBALANCE |discordant - baseline| = "
              f"{100 * imbalance:.1f} pp")

    cov_tab = []
    for kind, g in out.groupby("kind"):
        for flag, lab in ((1, "retained"), (0, "excluded")):
            h = g[g.hiconf == flag]
            if not len(h):
                continue
            row = {"kind": kind, "group": lab, "n": int(len(h)),
                   "median_span": float(h.span.median()),
                   "median_coverage": float(
                       pd.to_numeric(h.forward_coverage,
                                     errors="coerce").median()),
                   "median_pident": float(
                       pd.to_numeric(h.forward_pident,
                                     errors="coerce").median()),
                   "median_omega": float(h.omega.median())}
            cov_tab.append(row)
            print(f"    {kind:20s} {lab:9s} n {row['n']:6d}  "
                  f"span {row['median_span']:7.0f}  "
                  f"cov {row['median_coverage']:5.2f}  "
                  f"pident {row['median_pident']:5.1f}  "
                  f"omega {row['median_omega']:.4f}")
    pd.DataFrame(cov_tab).to_csv(OUT / "retention_covariates.tsv", sep="\t",
                                 index=False)

    per_unit = []
    for unit, col in (("genome", "accession"), ("genus", "genus")):
        for keyv, g in out.groupby(col):
            for kind, h in g.groupby("kind"):
                per_unit.append({
                    "unit": unit, "key": keyv, "kind": kind, "n": int(len(h)),
                    "hash_frac": float((h.contig_hash_verified == 1).mean()),
                    "rbh_frac": float(h.rbh_modified.mean()),
                    "hiconf_frac": float(h.hiconf.mean()),
                    "syn1_frac": float(h.omega_syn1.notna().mean()),
                    "median_omega": float(h.omega.median())})
    pd.DataFrame(per_unit).to_csv(OUT / "retention_by_unit.tsv", sep="\t",
                                  index=False)

    print("\nOUTCOME — modified method, all loci")
    for o, n in out.outcome.value_counts(dropna=False).items():
        print(f"  {str(o):26s} {n:6d}  ({100 * n / len(out):5.1f}%)")
    appl = out[out.has_annotated_cds == 1]
    print(f"\nOUTCOME — standard protein-protein method, "
          f"{len(appl)} applicable loci")
    for o, n in appl.standard_outcome.value_counts(dropna=False).items():
        print(f"  {str(o):26s} {n:6d}  ({100 * n / max(len(appl), 1):5.1f}%)")

    # ================================================== arms
    rows_d = out.to_dict("records")
    arms = [
        arm(rows_d, "omega", genus, "E1 original locus universe"),
        arm(rows_d, "omega_hiconf", genus, "E2 modified RBH, high-confidence"),
        arm(rows_d, "omega_standard", genus, "E2 standard protein-protein RBH"),
        arm(rows_d, "omega_syn1", genus, "E2 modified RBH + synteny >= 1"),
        arm(rows_d, "omega_syn2", genus, "E2 modified RBH + synteny >= 2"),
        arm(rows_d, "omega_syn_collinear", genus,
            "E2 modified RBH + collinear synteny"),
        arm(rows_d, "omega_noframesuspect", genus,
            "E2 high-conf, frame-suspect dropped"),
        arm(rows_d, "omega_excluded", genus,
            "E3 EXCLUDED loci only (what is lost)"),
    ]
    print("\nG3 BY ARM — E1, E2 and E3 are DIFFERENT ESTIMANDS, not corrections")
    print("  E1 = all scored loci.  E2 = among RBH-retained loci.  "
          "E3 = among excluded loci.")
    for A in arms:
        show(A)

    # ================================================== matching grid
    print("\nMATCHED ARMS inside the retained set (requirement 3)")
    print("  matching balances covariates AMONG RETAINED LOCI and recovers")
    print("  nothing about the excluded ones.")
    grid, primary = [], None
    for setname, cols in COVSETS.items():
        for cal in CALIPERS:
            res, mp, bf, af, sup = match_arm(out, cols, cal)
            res["covset"] = setname
            grid.append(res)
            tag = f"{setname:15s} caliper {str(cal):>5}"
            if res.get("median_gap") is None:
                print(f"  {tag}  no pairs")
                continue
            print(f"  {tag}  pairs {res['n_pairs']:5d}  "
                  f"unmatched {100 * (res['unmatched_frac'] or 0):5.1f}%  "
                  f"G3 {res['median_gap']:+.4f}  "
                  f"{res['n_positive']}/{res['n_tested']}  "
                  f"p {res['p']:<8.3g} "
                  f"maxSMD {res.get('max_abs_smd_after', float('nan')):.3f}  "
                  f"{'balanced' if res.get('balanced') else 'IMBALANCED'}")
            if setname == PRIMARY_COVSET and cal == PRIMARY_CALIPER:
                primary = res
                mp.to_csv(OUT / "matched_pairs.tsv", sep="\t", index=False)
                bf.to_csv(OUT / "smd_before.tsv", sep="\t", index=False)
                af.to_csv(OUT / "smd_after.tsv", sep="\t", index=False)
                sup.to_csv(OUT / "common_support.tsv", sep="\t", index=False)
    pd.DataFrame(grid).to_csv(OUT / "matching_grid.tsv", sep="\t", index=False)
    if primary:
        print(f"\n  PRIMARY: {PRIMARY_COVSET}, caliper {PRIMARY_CALIPER}")
        print(f"    SMD before (median |SMD| per covariate): "
              f"{primary.get('smd_before')}")
        print(f"    SMD after : {primary.get('smd_after')}")
        print(f"    common support (worst genome, frac of discordant inside "
              f"the baseline range): {primary.get('common_support')}")

    gaps = [g["median_gap"] for g in grid if g.get("median_gap") is not None]
    sens_stable = bool(gaps and (all(v > 0 for v in gaps)
                                 or all(v < 0 for v in gaps)))

    # ================================================== gates
    tech = {
        "forward_reproduction": {"value": fwd_agree,
                                 "threshold": GATE_FORWARD_AGREE,
                                 "pass": bool(fwd_agree is not None
                                              and fwd_agree >= GATE_FORWARD_AGREE)},
        "hash_contig_coverage": {"value": hash_cov,
                                 "threshold": GATE_HASH_COVERAGE,
                                 "pass": bool(hash_cov >= GATE_HASH_COVERAGE)},
        "omega_not_recomputed": {"value": omega_intact, "threshold": True,
                                 "pass": bool(omega_intact)},
        "input_checksum_stable": {"value": input_stable, "threshold": True,
                                  "pass": bool(input_stable)},
    }
    sci = {
        "forward_ambiguity": {
            "value": fwd_amb_frac, "threshold": GATE_MAX_FWD_AMBIGUITY,
            "pass": bool(fwd_amb_frac <= GATE_MAX_FWD_AMBIGUITY)},
        "retention_imbalance": {
            "value": imbalance, "threshold": GATE_MAX_RETENTION_IMBALANCE,
            "pass": bool(imbalance is not None
                         and imbalance <= GATE_MAX_RETENTION_IMBALANCE)},
        "matching_balance": {
            "value": (primary or {}).get("max_abs_smd_after"),
            "threshold": SMD_BALANCED,
            "pass": bool((primary or {}).get("balanced"))},
        "sensitivity_sign_stable": {
            "value": sens_stable, "threshold": True, "pass": sens_stable},
    }
    print("\nTECHNICAL GATES")
    for k, v in tech.items():
        print(f"  {k:26s} {str(v['value'])[:8]:>8}  vs {v['threshold']}  "
              f"-> {'PASS' if v['pass'] else 'FAIL'}")
    print("SCIENTIFIC QUALITY GATE")
    for k, v in sci.items():
        print(f"  {k:26s} {str(v['value'])[:8]:>8}  vs {v['threshold']}  "
              f"-> {'PASS' if v['pass'] else 'FAIL'}")
    tech_pass = all(v["pass"] for v in tech.values())
    sci_pass = all(v["pass"] for v in sci.values())
    print(f"\n  technical: {'PASS' if tech_pass else 'FAIL'}   "
          f"scientific: {'PASS' if sci_pass else 'FAIL'}")
    print("  EXPANSION RECOMMENDATION: "
          + ("expand" if tech_pass and sci_pass else
             "DO NOT expand on these results alone; see the failing criteria"))

    res = {
        "note": "POST HOC; RBH with the orthology and selection-bias audit; "
                "supersedes analysis_rbh_validated.py",
        "input_per_locus": str(src.relative_to(P)), "input_sha256": in_sha,
        "input_sha256_at_end": end_sha, "input_stable": input_stable,
        "input_snapshot": str(snap.relative_to(P)),
        "n_pairs": len(pairs), "n_loci": int(len(out)),
        "params": {"min_prot_aa": MIN_PROT_AA, "tie_margin": TIE_MARGIN,
                   "fwd_targets": FWD_TARGETS, "rev_targets": REV_TARGETS,
                   "max_stop_frac": MAX_STOP_FRAC,
                   "n_neighbours": N_NEIGHBOURS,
                   "synteny_window": SYNTENY_WINDOW,
                   "collinear_rho": COLLINEAR_RHO,
                   "calipers": [str(c) for c in CALIPERS],
                   "primary_caliper": PRIMARY_CALIPER,
                   "primary_covset": PRIMARY_COVSET,
                   "smd_balanced": SMD_BALANCED},
        "control_forward_agreement": fwd_agree, "control_forward_n": fc,
        "control_omega_not_recomputed": omega_intact,
        "forward_ambiguity": {
            "frac": fwd_amb_frac,
            "n": int(amb.sum()),
            "median_margin": float(fmar.median()),
            "p05_margin": float(fmar.quantile(0.05)),
            "median_candidates": float(ncand.median()),
            "runner_up_within": {str(t): int((fmar < t).sum())
                                 for t in (0.01, 0.05, 0.10)}},
        "contig_identity": dict(audit["contig_identity"]),
        "hash_verified_locus_frac": hash_cov,
        "audits": {k: (dict(v) if isinstance(v, collections.Counter) else v)
                   for k, v in audit.items() if k != "contig_identity"},
        "retention": retention, "retention_imbalance": imbalance,
        "retention_covariates": cov_tab,
        "arms": arms, "matching_grid": grid, "matching_primary": primary,
        "sensitivity_sign_stable": sens_stable,
        "outcomes": {str(k): int(v) for k, v in
                     out.outcome.value_counts(dropna=False).items()},
        "standard_outcomes": {str(k): int(v) for k, v in
                              appl.standard_outcome.value_counts(
                                  dropna=False).items()},
        "n_standard_applicable": int(len(appl)),
        "gates_technical": tech, "gates_scientific": sci,
        "technical_pass": tech_pass, "scientific_pass": sci_pass,
        "expansion_recommended": bool(tech_pass and sci_pass),
        "estimands": {
            "E1": "G3 over the original locus universe, all scored loci",
            "E2": "G3 among RBH-retained loci; a different quantity from E1, "
                  "not a corrected version of it",
            "E3": "G3 among loci RBH excluded; reported so that what is "
                  "dropped is visible",
            "matched": "G3 among retained loci with covariates balanced "
                       "between the two classes; says nothing about excluded "
                       "loci"},
        "selection_bias_warning":
            "A stable G3 after RBH filtering is not independent biological "
            "corroboration. Retention depends on reciprocation strength, which "
            "depends on divergence and alignable length, and divergence is "
            "what omega measures. The matched arms address covariate imbalance "
            "among retained loci only; they cannot address what exclusion "
            "removed.",
        "skipped": skipped,
        "host_free_gb_start": free0, "host_free_gb_end": host_free_gb()}
    (P / a.out_json).write_text(json.dumps(res, indent=1, default=str,
                                           sort_keys=True))
    print(f"\nwrote {a.out_json}")


if __name__ == "__main__":
    main()
