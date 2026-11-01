#!/usr/bin/env python3
"""Adversarially validated reciprocal-best-hit orthology.

POST HOC. Not pre-registered. No frozen plan is touched or reinterpreted.
Supersedes analysis_rbh_orthology.py, whose result is kept as a baseline and
cross-check rather than deleted.

WHY A SECOND VERSION. The first version produced a number (97.1 % reciprocation
on the smoke target) that looked good, and a number that looks good is the most
dangerous output this project can produce. Five specific ways it could have been
wrong were raised, and each is now answered by construction rather than by
argument:

  1 CONTIG IDENTITY. Version 1 reconciled GenBank records to locus contigs by
    sequence LENGTH, with record ORDER as a fall-back. Length collisions are
    real in multi-replicon assemblies and order is not evidence of anything.
    Here the mapping is by SHA-256 of the sequence itself, required to be unique
    on both sides. There is no order fall-back into the high-confidence cohort;
    loci whose contig cannot be hash-identified are reported as an excluded
    fraction and carried only in the exploratory arm.

  2 STANDARD VERSUS MODIFIED RBH. Textbook RBH is protein against protein.
    Version 1 ran a MODIFIED procedure in which a locus is represented by the
    translation of the span its forward hit aligned, because a pseudogene locus
    has no annotated protein. Both are computed here and never conflated:
      STANDARD   available only where the locus has an overlapping annotated
                 CDS; both directions are annotated proteins.
      MODIFIED   available for every locus, including those with no CDS.
    The applicability of the standard method is itself asymmetric between
    discordant and baseline loci, and that asymmetry is a reported result.

  3 BIOLOGICAL LOCI, NOT FASTA ENTRIES. A reverse search over a database that
    contains both a gene's annotated protein and the translated span of the
    locus covering it has two entries for one biological object. Collapsing
    them after the search - as version 1 did, with an overlap clause - is not
    the same as collapsing before the decision. Here every database entry
    carries a biological identifier, the best score per identifier is taken,
    and the decision compares identifiers. Duplicate annotated proteins
    (identical sequence) are collapsed the same way.

  4 TIES AND PARALOGUES. Version 1 took the first line DIAMOND printed. Here
    the reverse search returns up to 25 targets, the best and runner-up
    biological identifiers are both recorded with their bitscores, and a
    relative margin below 1 % is declared AMBIGUOUS and kept out of the
    high-confidence cohort instead of being resolved arbitrarily.

  5 FRAME SELECTION. The translated span inherits whatever frame the forward
    hit chose. Frame, orientation, internal-stop load and alignment coverage
    are recorded per locus, and a span whose translation is more than 10 %
    stop codons is flagged.

THE SELECTION-BIAS WARNING, WHICH IS THE POINT OF THE MATCHED ARM. A change in
G3 after RBH filtering is NOT evidence that the retained orthologues are the
correct ones. RBH retention depends on how well a locus reciprocates, which
depends on its divergence and its alignable length - and divergence is what
omega measures. Filtering on reciprocation is therefore filtering on a
correlate of the outcome, and it will move G3 whether or not the orthologue
assignments improve. That is why the fourth arm matches retained discordant
loci to retained baseline loci on span length, alignment coverage and percent
identity within each genome, so that the comparison is between loci the filter
treated alike. Nothing in this script licenses the inference "G3 survived
filtering, therefore the orthologues are right".

OMEGA IS NEVER RECOMPUTED. Every arm reuses the omega in the input per-locus
table, whose control column reproduced results/loci.tsv exactly (6360/6360).
Only the locus set changes between arms, and the script asserts this.

ISOLATION. Reads results/reference_sensitivity/per_locus.tsv (the 12-target
pilot) read-only. Writes only results/rbh2/ and its own JSON. Never touches
results/loci.tsv, results/locus.json, per_locus_all.tsv or anything else the
92-target reference expansion owns, and reads the shared partner reference only
behind its completion sentinel.

DISK. Start floor 10 GB, abort floor 8 GB, read from /mnt/c, never `df /`.

    python3 analysis_rbh_validated.py [--per-locus T] [--out-json J] [--targets N]
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
from Bio import SeqIO
from Bio.Seq import Seq
from scipy import stats

P = pathlib.Path(__file__).resolve().parent
R = P / "results"
OUT = R / "rbh2"
D = pathlib.Path.home() / "research/datasets"
REFDIR = D / "n2-refseq-reference"          # READ ONLY
QDIR = D / "n2-query-proteomes"             # this script's own cache

MIN_PROT_AA = 20
MAX_STOP_FRAC = 0.10        # above this the translated span is frame_suspect
TIE_MARGIN = 0.01           # relative bitscore margin below which it is a tie
REV_TARGETS = 25
N_NEIGHBOURS = 5
SYNTENY_WINDOW = 10
THREADS = 4
FLOOR_START_GB, FLOOR_ABORT_GB = 10, 8

# validation gates; expansion is allowed only if all pass
GATE_FORWARD_AGREE = 0.99
GATE_HASH_COVERAGE = 0.95


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


def translate_span(nt, qs, qe):
    """-> (protein, frame, reversed, stop_fraction) for the aligned span.

    Mirrors analysis_dnds.codon_pairs: slice on min/max of the query bounds and
    reverse-complement when the hit runs backwards. Stops become X so the entry
    is a legal DIAMOND database sequence, and their share is reported.
    """
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


# ----------------------------------------------------------- reference access
def partner_reference(partner):
    """(faa, ffn) from the shared cache, only behind its completion sentinel."""
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
    """-> [(id, prot, rec_id, start, end, rec_sha256, rec_idx)], non-pseudo CDS.

    rec_sha256 is the SHA-256 of the GenBank record's own sequence. It is what
    identifies the contig; length and record order are recorded nowhere here,
    because version 1 used them and that was the defect.
    """
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
    """GenBank record id -> locus contig name, by EXACT unique sequence hash.

    Requirement 1. No positional or order fall-back is offered: a record whose
    hash is absent, or ambiguous on either side, maps to None and every locus on
    it is excluded from the high-confidence cohort.
    """
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
    """-> {locus_tag: (record id, index within that record by coordinate)}.

    Both sides of the synteny test are in the partner's own GenBank namespace,
    so no cross-namespace reconciliation arises here.
    """
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


def blastx_best(query_fa, db_faa, workdir):
    """Forward: locus nucleotides -> partner proteins.

    Same parameters as analysis_dnds.best_hits (one target, e-value 1e-5, first
    line per query), with bitscore and coverage added, so the subject it returns
    must reproduce what the published pipeline recorded. That reproduction is
    control 1 and a gate.
    """
    workdir.mkdir(parents=True, exist_ok=True)
    dmnd = workdir / "db.dmnd"
    subprocess.run(["diamond", "makedb", "--in", str(db_faa), "-d", str(dmnd),
                    "--quiet", "--threads", str(THREADS)], check=True)
    out = workdir / "fwd.tsv"
    subprocess.run(["diamond", "blastx", "-q", str(query_fa), "-d", str(dmnd),
                    "-o", str(out), "--quiet", "--threads", str(THREADS),
                    "--max-target-seqs", "1", "--evalue", "1e-5",
                    "--outfmt", "6", *FWD_FMT], check=True)
    best = {}
    for line in out.read_text().splitlines():
        f = line.split("\t")
        if len(f) < len(FWD_FMT) or f[0] in best:
            continue
        best[f[0]] = {"subject": f[1], "pident": float(f[2]),
                      "aln_len": int(f[3]), "qstart": int(f[4]),
                      "qend": int(f[5]), "slen": int(f[6]),
                      "bitscore": float(f[7]),
                      "coverage": int(f[3]) / int(f[6]) if int(f[6]) else None}
    return best


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


# ------------------------------------------------- the decision, collapsed
def collapse(hits, entry_bio):
    """-> (best_bio, best_score, second_bio, second_score, margin).

    Requirement 3. Scores are maxed per BIOLOGICAL identifier before anything is
    compared, so two database entries for one gene cannot out-compete or tie
    with each other. Requirement 4: the runner-up is a DIFFERENT biological
    identifier, and the relative margin against it is what decides a tie.
    """
    per = {}
    for ent, score in hits:
        bio = entry_bio.get(ent)
        if bio is None:
            continue
        if score > per.get(bio, -1.0):
            per[bio] = score
    if not per:
        return None, None, None, None, None
    order = sorted(per.items(), key=lambda kv: -kv[1])
    b_bio, b_sc = order[0]
    if len(order) == 1:
        return b_bio, b_sc, None, None, 1.0
    s_bio, s_sc = order[1]
    margin = (b_sc - s_sc) / b_sc if b_sc > 0 else 0.0
    return b_bio, b_sc, s_bio, s_sc, margin


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

    # ---- requirement 1: contigs by exact hash, no order fall-back
    xlat, hstat = hash_contig_map(contigs, cds)
    audit["contig_identity"].update(hstat)

    # ---- query loci
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

    # ---- annotated CDS, in the locus contig namespace where hash-identified
    cds_iv, by_contig, prot_of = {}, collections.defaultdict(list), {}
    unmapped_cds = 0
    for ident, prot, rid, start, end, _, _ in cds:
        prot_of[ident] = prot
        contig = xlat.get(rid)
        if contig is None:
            unmapped_cds += 1
            continue
        cds_iv[ident] = (contig, start, end)
        by_contig[contig].append((start, ident))
    for v in by_contig.values():
        v.sort()
    audit["cds_unmapped_contig"] += unmapped_cds
    audit["cds_total"] += len(cds)

    # ---- requirement 2 audit: duplicate annotated proteins
    dup = collections.defaultdict(list)
    for ident, prot in prot_of.items():
        dup[sha(prot)].append(ident)
    dup_family = {i: fam for fam in dup.values() if len(fam) > 1 for i in fam}
    audit["duplicate_protein_entries"] += len(dup_family)

    # ---- which locus does each annotated CDS belong to, if any
    cds_locus = {}
    overlap_multi = 0
    for k, r in meta.items():
        for start, ident in by_contig.get(r.contig, []):
            iv = cds_iv[ident]
            if overlaps_half(r.start, r.end, iv[1], iv[2]):
                if ident in cds_locus:
                    overlap_multi += 1
                cds_locus[ident] = k
    audit["cds_claimed_by_two_loci"] += overlap_multi

    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        qf = td / "loci.fna"
        with qf.open("w") as fh:
            for k, v in queries.items():
                fh.write(f">{k}\n{v}\n")
        fwd = blastx_best(qf, faa, td / "fwd")
        if not fwd:
            return None, "forward search returned nothing"

        # ---- translated spans, with the frame audit
        cand = {}
        for k, h in fwd.items():
            prot, frame, rev, frac = translate_span(queries[k], h["qstart"],
                                                    h["qend"])
            h.update({"frame": frame, "reversed": int(rev),
                      "stop_frac": frac,
                      "frame_suspect": int(frac is not None
                                           and frac > MAX_STOP_FRAC)})
            if prot:
                cand[k] = prot
                audit["frames"][frame] += 1
                audit["reversed_hits"] += int(rev)
                audit["frame_suspect"] += h["frame_suspect"]
            else:
                audit["span_too_short_to_translate"] += 1

        # ---- the locus's own annotated protein, where one exists: the longest
        #      CDS overlapping it. This is the standard method's query.
        locus_cds = {}
        for ident, k in cds_locus.items():
            iv = cds_iv[ident]
            span = iv[2] - iv[1]
            prev = locus_cds.get(k)
            if prev is None or span > prev[1]:
                locus_cds[k] = (ident, span)

        # ---- standard RBH, forward leg: annotated protein -> partner proteins
        std_fwd = {}
        if locus_cds:
            stdq = td / "stdq.faa"
            with stdq.open("w") as fh:
                for k, (ident, _) in locus_cds.items():
                    fh.write(f">{k}\n{prot_of[ident]}\n")
            for k, hs in blastp_all(stdq, faa, td / "stdfwd",
                                    max_targets=1).items():
                if hs:
                    std_fwd[k] = hs[0][0]

        # ---- reverse database, every entry carrying a biological identifier.
        #      entry_bio collapses FASTA entries to biological objects, and
        #      entry_bio_cds is the same map without the translated spans, so
        #      the standard method stays protein-against-protein throughout.
        entry_bio, entry_bio_cds = {}, {}
        revdb = td / "revdb.faa"
        with revdb.open("w") as fh:
            for ident, prot in prot_of.items():
                if ident not in cds_iv:
                    continue            # contig not hash-identified
                bio = cds_locus.get(ident)
                if bio is None:
                    # its own gene; a duplicate family is one biological object
                    bio = (f"DUP:{sha(prot)[:12]}" if ident in dup_family
                           else f"CDS:{ident}")
                entry_bio[f"CDS|{ident}"] = bio
                entry_bio_cds[f"CDS|{ident}"] = bio
                fh.write(f">CDS|{ident}\n{prot}\n")
            for k, prot in cand.items():
                entry_bio[f"LOC|{k}"] = k
                fh.write(f">LOC|{k}\n{prot}\n")

        # ---- one reverse search serves both methods: its query set is the
        #      union of the subjects either forward leg produced.
        subjects = sorted({h["subject"] for h in fwd.values()}
                          | set(std_fwd.values()))
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

        # ---- synteny map: whole annotated proteome -> partner proteome
        qpf = td / "qprot.faa"
        with qpf.open("w") as fh:
            for ident in cds_iv:
                fh.write(f">{ident}\n{prot_of[ident]}\n")
        q2p_all = blastp_all(qpf, faa, td / "syn", max_targets=1)
        q2p = {q: hs[0][0] for q, hs in q2p_all.items() if hs}

    if partner not in syn_cache:
        syn_cache[partner] = partner_order(partner)
    porder = syn_cache[partner]

    # ---- decide, locus by locus
    hash_ok = {v for v in xlat.values() if v}
    rows, fa, fc = [], 0, 0
    for k, r in meta.items():
        h = fwd.get(k)
        rec = {
            "accession": acc, "genus": r.genus, "partner": partner,
            "contig": r.contig, "start": r.start, "end": r.end,
            "span": int(r.end) - int(r.start) + 1, "kind": r.kind,
            "omega": r.control,
            "contig_hash_verified": int(r.contig in hash_ok),
            "forward_subject": None, "forward_bitscore": None,
            "forward_pident": None, "forward_coverage": None,
            "frame": None, "reversed": None, "stop_frac": None,
            "frame_suspect": None,
            "recorded_subject": clean(getattr(r, "treatment_subject", None)),
            "has_annotated_cds": int(k in locus_cds),
            "reverse_best_bio": None, "reverse_best_score": None,
            "reverse_second_bio": None, "reverse_second_score": None,
            "reverse_margin": None,
            "rbh_modified": 0, "rbh_standard": None, "ambiguous_tie": 0,
            "outcome": None, "synteny_votes": None,
            "standard_subject": std_fwd.get(k),
            "standard_outcome": None, "standard_margin": None,
        }
        if h is None:
            rec["outcome"] = "no_forward_hit"
            rows.append(rec)
            continue
        rec.update({"forward_subject": h["subject"],
                    "forward_bitscore": h["bitscore"],
                    "forward_pident": h["pident"],
                    "forward_coverage": h["coverage"],
                    "frame": h["frame"], "reversed": h["reversed"],
                    "stop_frac": h["stop_frac"],
                    "frame_suspect": h["frame_suspect"]})
        if rec["recorded_subject"]:
            fc += 1
            fa += int(h["subject"] == rec["recorded_subject"])

        b, bs, s2, ss, margin = collapse(rev.get(h["subject"], []), entry_bio)
        rec.update({"reverse_best_bio": b, "reverse_best_score": bs,
                    "reverse_second_bio": s2, "reverse_second_score": ss,
                    "reverse_margin": margin})
        if b is None:
            rec["outcome"] = "no_reverse_hit"
        elif b != k:
            rec["outcome"] = "reverse_hit_other_locus"
        elif margin is not None and margin < TIE_MARGIN:
            rec["outcome"] = "ambiguous_tie"
            rec["ambiguous_tie"] = 1
        else:
            rec["outcome"] = "rbh_confirmed"
            rec["rbh_modified"] = 1

        # ---- the standard method, decided on its own terms: annotated protein
        #      forward, annotated proteins only in reverse, same tie policy.
        if k in locus_cds:
            s_sub = std_fwd.get(k)
            if s_sub is None:
                rec["standard_outcome"] = "no_forward_hit"
                rec["rbh_standard"] = 0
            else:
                sb, _, _, _, smar = collapse(rev.get(s_sub, []),
                                             entry_bio_cds)
                rec["standard_margin"] = smar
                if sb is None:
                    rec["standard_outcome"] = "no_reverse_hit"
                elif sb != k:
                    rec["standard_outcome"] = "reverse_hit_other_locus"
                elif smar is not None and smar < TIE_MARGIN:
                    rec["standard_outcome"] = "ambiguous_tie"
                else:
                    rec["standard_outcome"] = "rbh_confirmed"
                rec["rbh_standard"] = int(
                    rec["standard_outcome"] == "rbh_confirmed")

        if rec["rbh_modified"]:
            tgt = porder.get(h["subject"])
            votes = 0
            if tgt:
                near = sorted(by_contig.get(r.contig, []),
                              key=lambda t: abs(t[0] - r.start))
                taken = 0
                for start, ident in near:
                    iv = cds_iv[ident]
                    if overlaps_half(r.start, r.end, iv[1], iv[2]):
                        continue
                    pt = porder.get(q2p.get(ident, ""))
                    taken += 1
                    if pt and pt[0] == tgt[0] and \
                            abs(pt[1] - tgt[1]) <= SYNTENY_WINDOW:
                        votes += 1
                    if taken >= N_NEIGHBOURS:
                        break
            rec["synteny_votes"] = votes
        rows.append(rec)

    return (rows, fa, fc), None


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
        print(f"  {A['label']:46s} {A.get('note')}")
        return
    ci = A.get("genus_cluster_ci95") or [None, None]
    line = (f"  {A['label']:46s} {g:+.4f}  {A['n_positive']}/{A['n_tested']}"
            f"  p {A['p']:<9.3g} n_loci {A['n_loci_used']:6d}")
    if ci[0] is not None:
        line += f"  genus CI [{ci[0]:+.4f}, {ci[1]:+.4f}]"
    print(line)


def matched_pairs(df, seed=20261031):
    """Requirement 5, fourth arm: nearest-neighbour matching within a genome.

    Retained discordant loci are matched 1:1, without replacement, to retained
    baseline loci on z-scored log span, alignment coverage and percent identity.
    Both sides have passed the same filter, so differential retention cannot
    drive the comparison, and the covariates the filter selects on are balanced
    by construction. Matching is deterministic given the seed.
    """
    g = np.random.default_rng(seed)
    keep = df[(df.hiconf == 1)].copy()
    cols = ["span", "forward_coverage", "forward_pident"]
    keep = keep.dropna(subset=cols + ["omega"])
    pairs, balance = [], []
    for acc, sub in keep.groupby("accession"):
        dis = sub[sub.kind == "discordant"]
        con = sub[sub.kind == "concordant_intact"]
        if len(dis) < 3 or len(con) < 3:
            continue
        Z = sub[cols].copy()
        Z["span"] = np.log10(Z["span"].clip(lower=1))
        mu, sd = Z.mean(), Z.std().replace(0, 1)
        zd = ((dis[cols].assign(span=np.log10(dis["span"].clip(lower=1)))
               - mu) / sd).to_numpy()
        zc = ((con[cols].assign(span=np.log10(con["span"].clip(lower=1)))
               - mu) / sd).to_numpy()
        taken = set()
        order = g.permutation(len(zd))      # order of matching, not the match
        for i in order:
            d = ((zc - zd[i]) ** 2).sum(axis=1)
            d[list(taken)] = np.inf
            j = int(np.argmin(d))
            if not np.isfinite(d[j]):
                break
            taken.add(j)
            pairs.append({"accession": acc, "genus": dis.genus.iloc[0],
                          "omega_dis": float(dis.omega.iloc[i]),
                          "omega_con": float(con.omega.iloc[j]),
                          "dist": float(np.sqrt(d[j]))})
        if taken:
            sel = con.iloc[sorted(taken)]
            for c in cols:
                balance.append({"accession": acc, "covariate": c,
                                "discordant_median": float(dis[c].median()),
                                "baseline_all_median": float(con[c].median()),
                                "baseline_matched_median": float(sel[c].median())})
    return pd.DataFrame(pairs), pd.DataFrame(balance)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-locus",
                    default="results/reference_sensitivity/per_locus.tsv")
    ap.add_argument("--out-json", default="results/rbh_validated.json")
    ap.add_argument("--targets", type=int, default=0)
    a = ap.parse_args()

    free0 = guard("start")
    print(f"host C: free {free0} GB (from /mnt/c)")
    if free0 < FLOOR_START_GB:
        sys.exit(f"refusing to start: {free0} < {FLOOR_START_GB} GB")
    OUT.mkdir(parents=True, exist_ok=True)

    src = P / a.per_locus
    raw = src.read_bytes()
    in_sha = hashlib.sha256(raw).hexdigest()
    print(f"input {src.name}  sha256 {in_sha[:16]}…  {len(raw)} bytes")
    df = pd.read_csv(src, sep="\t")
    df = df[df.control.notna()]
    genus = dict(zip(df.accession, df.genus))
    pairs = df[["accession", "partner"]].drop_duplicates().values.tolist()
    if a.targets:
        pairs = pairs[:a.targets]
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
        n = sum(r["rbh_modified"] for r in rows)
        print(f"  [{i}/{len(pairs)}] {acc}: {len(rows)} loci, "
              f"{n} modified-RBH ({100 * n / len(rows):.1f}%), "
              f"{sum(r['ambiguous_tie'] for r in rows)} ties", flush=True)
        guard(acc)

    if not all_rows:
        sys.exit("nothing processed")
    out = pd.DataFrame(all_rows)

    # ---- high-confidence cohort: hash-verified contig, no tie, confirmed
    out["hiconf"] = ((out.contig_hash_verified == 1) &
                     (out.rbh_modified == 1) &
                     (out.ambiguous_tie == 0)).astype(int)
    votes = pd.to_numeric(out.synteny_votes, errors="coerce").fillna(0)
    out["omega_hiconf"] = np.where(out.hiconf == 1, out.omega, np.nan)
    out["omega_syn1"] = np.where((out.hiconf == 1) & (votes >= 1),
                                 out.omega, np.nan)
    out["omega_syn2"] = np.where((out.hiconf == 1) & (votes >= 2),
                                 out.omega, np.nan)
    out["omega_standard"] = np.where(out.rbh_standard == 1, out.omega, np.nan)
    out.to_csv(OUT / "per_locus_rbh_validated.tsv", sep="\t", index=False)
    print(f"\nwrote {len(out)} locus rows to "
          f"{OUT / 'per_locus_rbh_validated.tsv'}")

    # ================================================== controls and gates
    fwd_agree = fa / fc if fc else None
    hash_cov = float((out.contig_hash_verified == 1).mean())
    print("\nCONTROL 1 — re-run forward search vs recorded treatment subject")
    print(f"  {fa} of {fc} identical "
          f"({100 * fwd_agree:.2f}%)" if fwd_agree is not None else "  n/a")
    print("\nREQUIREMENT 1 — contig identity by exact sequence hash")
    for k, v in sorted(audit["contig_identity"].items()):
        print(f"  {k:30s} {v:6d} records")
    print(f"  loci on hash-verified contigs   {100 * hash_cov:6.2f}%")
    print(f"  EXCLUDED from high-confidence   {100 * (1 - hash_cov):6.2f}% "
          f"({int((1 - hash_cov) * len(out))} loci)")
    print("  no positional or order fall-back was used anywhere")

    # ---- control 2: omega must be bit-identical to the input's control column
    key = ["accession", "contig", "start", "end"]
    chk = out[key + ["omega"]].merge(df[key + ["control"]], on=key, how="left")
    omega_intact = bool(len(chk) == len(out) and np.array_equal(
        chk.omega.to_numpy(float), chk.control.to_numpy(float),
        equal_nan=True))
    print(f"\nCONTROL 2 — omega bit-identical to the input, never recomputed: "
          f"{omega_intact}  ({len(chk)} rows compared)")

    print("\nREQUIREMENT 2 — audits")
    for k in ("cds_total", "cds_unmapped_contig", "duplicate_protein_entries",
              "cds_claimed_by_two_loci", "span_too_short_to_translate",
              "reversed_hits", "frame_suspect"):
        print(f"  {k:30s} {audit[k]:6d}")
    print(f"  forward frame distribution      {dict(audit['frames'])}")
    print(f"  loci with an annotated CDS      "
          f"{int(out.has_annotated_cds.sum())} of {len(out)}")
    for kind, g in out.groupby("kind"):
        print(f"    {kind:20s} {100 * g.has_annotated_cds.mean():5.1f}% "
              f"have one  -> standard RBH is only applicable there")

    print("\nREQUIREMENT 4 — retention by locus class")
    retention = {}
    for kind, g in out.groupby("kind"):
        retention[kind] = {
            "n": int(len(g)),
            "hash_verified_frac": float((g.contig_hash_verified == 1).mean()),
            "modified_rbh_frac": float(g.rbh_modified.mean()),
            "ambiguous_tie_frac": float(g.ambiguous_tie.mean()),
            "hiconf_frac": float(g.hiconf.mean()),
            "syn1_frac": float(g.omega_syn1.notna().mean()),
            "syn2_frac": float(g.omega_syn2.notna().mean()),
            "standard_rbh_frac_where_applicable": (
                float(pd.to_numeric(
                    g[g.has_annotated_cds == 1].rbh_standard,
                    errors="coerce").mean())
                if (g.has_annotated_cds == 1).any() else None)}
        r = retention[kind]
        print(f"  {kind:20s} n {r['n']:6d}  hash {100*r['hash_verified_frac']:5.1f}%"
              f"  RBH {100*r['modified_rbh_frac']:5.1f}%"
              f"  tie {100*r['ambiguous_tie_frac']:4.1f}%"
              f"  hi-conf {100*r['hiconf_frac']:5.1f}%"
              f"  +syn1 {100*r['syn1_frac']:5.1f}%"
              f"  +syn2 {100*r['syn2_frac']:5.1f}%")
    d, b = retention.get("discordant"), retention.get("concordant_intact")
    if d and b:
        print(f"  retention imbalance, discordant minus baseline: "
              f"{100 * (d['hiconf_frac'] - b['hiconf_frac']):+.1f} pp")

    print("\n  WHY: covariates of retained vs rejected, by class")
    cov_tab = []
    for kind, g in out.groupby("kind"):
        for flag, lab in ((1, "retained"), (0, "rejected")):
            h = g[g.hiconf == flag]
            if not len(h):
                continue
            row = {"kind": kind, "group": lab, "n": int(len(h)),
                   "median_span": float(h.span.median()),
                   "median_coverage": float(h.forward_coverage.median()),
                   "median_pident": float(h.forward_pident.median()),
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
        for key, g in out.groupby(col):
            for kind, h in g.groupby("kind"):
                per_unit.append({"unit": unit, "key": key, "kind": kind,
                                 "n": int(len(h)),
                                 "hiconf_frac": float(h.hiconf.mean())})
    pd.DataFrame(per_unit).to_csv(OUT / "retention_by_unit.tsv", sep="\t",
                                  index=False)
    print(f"  per-genome and per-genus retention -> "
          f"{OUT / 'retention_by_unit.tsv'}")

    print("\nOUTCOME PER LOCUS — modified method (all loci)")
    for o, n in out.outcome.value_counts(dropna=False).items():
        print(f"  {str(o):26s} {n:6d}  ({100 * n / len(out):5.1f}%)")
    appl = out[out.has_annotated_cds == 1]
    print(f"\nOUTCOME PER LOCUS — standard protein-protein method "
          f"({len(appl)} applicable loci)")
    for o, n in appl.standard_outcome.value_counts(dropna=False).items():
        print(f"  {str(o):26s} {n:6d}  ({100 * n / max(len(appl), 1):5.1f}%)")

    # ================================================== requirement 5, arms
    rows_d = out.to_dict("records")
    arms = [arm(rows_d, "omega", genus, "1 original locus universe"),
            arm(rows_d, "omega_hiconf", genus, "2 high-confidence RBH"),
            arm(rows_d, "omega_syn1", genus, "3 high-conf RBH + synteny >= 1"),
            arm(rows_d, "omega_syn2", genus, "3b high-conf RBH + synteny >= 2"),
            arm(rows_d, "omega_standard", genus,
                "aux standard protein-protein RBH only")]
    print("\nG3 BY ARM (omega never recomputed; only the locus set changes)")
    for A in arms:
        show(A)

    mp, bal = matched_pairs(out)
    matched = {"n_pairs": int(len(mp))}
    if len(mp):
        mp.to_csv(OUT / "matched_pairs.tsv", sep="\t", index=False)
        bal.to_csv(OUT / "matched_balance.tsv", sep="\t", index=False)
        per = [float(np.median(g.omega_dis) - np.median(g.omega_con))
               for _, g in mp.groupby("accession")]
        nz = [v for v in per if v != 0]
        pos = sum(1 for v in nz if v > 0)
        matched.update({
            "n_genomes": int(len(per)),
            "median_gap": float(np.median(per)),
            "n_positive": pos, "n_tested": len(nz),
            "p": float(stats.binomtest(pos, len(nz), 0.5).pvalue) if nz else None,
            "median_match_distance": float(mp.dist.median())})
        print("\n4 COVARIATE-MATCHED ARM (retained loci only, 1:1 within genome)")
        print(f"  {len(mp)} matched pairs over {len(per)} genomes, "
              f"median match distance {mp.dist.median():.3f} z-units")
        print(f"  G3 {matched['median_gap']:+.4f}  "
              f"{pos}/{len(nz)} genomes positive  p {matched['p']:.3g}")
        print("  covariate balance -> " + str(OUT / "matched_balance.tsv"))

    gates = {
        "forward_reproduction": {"value": fwd_agree,
                                 "threshold": GATE_FORWARD_AGREE,
                                 "pass": bool(fwd_agree is not None
                                              and fwd_agree >= GATE_FORWARD_AGREE)},
        "hash_contig_coverage": {"value": hash_cov,
                                 "threshold": GATE_HASH_COVERAGE,
                                 "pass": bool(hash_cov >= GATE_HASH_COVERAGE)},
        "omega_not_recomputed": {"value": omega_intact, "threshold": True,
                                 "pass": bool(omega_intact)},
    }
    print("\nVALIDATION GATES (expansion is allowed only if all pass)")
    for k, v in gates.items():
        print(f"  {k:26s} {str(v['value']):>8}  vs {v['threshold']}  "
              f"-> {'PASS' if v['pass'] else 'FAIL'}")
    all_pass = all(v["pass"] for v in gates.values())
    print(f"  overall: {'ALL GATES PASS' if all_pass else 'GATE FAILURE'}")

    res = {"note": "POST HOC; adversarially validated RBH; supersedes "
                   "analysis_rbh_orthology.py",
           "input_per_locus": str(src.relative_to(P)), "input_sha256": in_sha,
           "n_pairs": len(pairs), "n_loci": int(len(out)),
           "params": {"min_prot_aa": MIN_PROT_AA, "tie_margin": TIE_MARGIN,
                      "rev_targets": REV_TARGETS,
                      "max_stop_frac": MAX_STOP_FRAC,
                      "n_neighbours": N_NEIGHBOURS,
                      "synteny_window": SYNTENY_WINDOW},
           "control_forward_agreement": fwd_agree, "control_forward_n": fc,
           "control_omega_not_recomputed": omega_intact,
           "contig_identity": dict(audit["contig_identity"]),
           "hash_verified_locus_frac": hash_cov,
           "excluded_locus_frac": 1 - hash_cov,
           "audits": {k: (dict(v) if isinstance(v, collections.Counter) else v)
                      for k, v in audit.items() if k != "contig_identity"},
           "retention": retention, "retention_covariates": cov_tab,
           "arms": arms, "matched_arm": matched,
           "outcomes": {str(k): int(v) for k, v in
                        out.outcome.value_counts(dropna=False).items()},
           "standard_outcomes": {
               str(k): int(v) for k, v in
               out[out.has_annotated_cds == 1].standard_outcome
               .value_counts(dropna=False).items()},
           "n_standard_applicable": int((out.has_annotated_cds == 1).sum()),
           "gates": gates, "all_gates_pass": all_pass,
           "skipped": skipped,
           "selection_bias_warning":
               "A change in G3 between arms is a change in the sample, not "
               "evidence that retained orthologue assignments are correct. RBH "
               "retention depends on reciprocation strength, which depends on "
               "divergence and alignable length, and divergence is what omega "
               "measures. The matched arm exists because of this.",
           "host_free_gb_start": free0, "host_free_gb_end": host_free_gb()}
    (P / a.out_json).write_text(json.dumps(res, indent=1, default=str,
                                           sort_keys=True))
    print(f"\nwrote {a.out_json}")


if __name__ == "__main__":
    main()
