#!/usr/bin/env python3
"""Reciprocal-best-hit orthology, with synteny, as a filter on G3.

POST HOC. Not pre-registered. No frozen plan is touched or reinterpreted.

THE ASSUMPTION UNDER TEST. Every omega in this project comes from a ONE-WAY
DIAMOND best hit: the locus's nucleotides are searched against the partner's
proteins and the top protein is called the orthologue. A one-way best hit is not
orthology. If the locus's true counterpart is absent from the partner (lost, or
filtered out of the reference) the search still returns something - most often a
paralogue - and a paralogous pair inflates omega for reasons that have nothing
to do with pseudogenisation. FINAL_INDEPENDENT_VALIDATION.md ranks this as the
largest standing assumption in the analysis.

THE TEST. Reciprocal best hits. The locus L maps forward to partner protein S;
S must map back to L and to nothing better in L's own genome. Then, separately,
synteny: L's genomic neighbours should map to S's genomic neighbours.

THE TRAP, AND HOW IT IS AVOIDED. The obvious implementation - reverse-search S
against the query genome's annotated proteome - is INVALID here, and would have
manufactured the result it was meant to check. A discordant locus is a candidate
pseudogene; RefSeq annotates pseudogenes as `gene` features with no CDS and no
/translation, so such a locus is ABSENT from its own genome's proteome. The
reverse search could never return it, every discordant locus would be rejected,
the concordant-intact baseline would survive, and the surviving comparison would
be structurally rigged. The asymmetry would look like a finding.

So the reverse database here is the query proteome PLUS the loci themselves,
each represented by the translation of the span that the forward search actually
aligned, in that hit's frame. Every locus is therefore present in the reverse
database on the same footing whether or not its genome's annotation grants it a
protein. A locus is accepted when the reverse best hit is its own translated
span, OR an annotated CDS that overlaps it by at least 50 % reciprocally - the
second clause matters because an intact locus is usually also an annotated gene,
so the two entries are near-duplicates and a tie between them must not count as
a rejection.

SYNTENY. A second search maps the whole query proteome to the partner proteome,
giving a gene-order correspondence. For an accepted locus, its 5 nearest
annotated neighbours that have partners of their own are checked against the
partner's gene order: a neighbour whose partner lies within 10 genes of S on the
same partner contig is a syntenic vote. Both >= 1 and >= 2 votes are reported;
neither threshold is privileged.

WHAT IS MEASURED, AND WHAT IS NOT ASSUMED. G3 on all scored loci, on
RBH-confirmed loci, and on RBH-plus-synteny loci, with the attrition reported
SEPARATELY FOR DISCORDANT AND BASELINE LOCI, because asymmetric attrition would
bias the comparison whichever way it fell. No assumption is made about whether
verified orthology should raise, lower or preserve G3. Direction and magnitude
are measured. A filter that removes most of one side and little of the other
invalidates its own comparison, and that is reported as a result rather than
explained away.

CONTROLS. Two, both internal. (1) The forward search is re-run here rather than
read from disk, so its subjects must reproduce the treatment subjects recorded
by analysis_reference_sensitivity.py; the agreement rate is reported and a low
one invalidates the run. (2) omega is never recomputed - every arm reuses the
omega already in the per-locus table, whose control column reproduced
results/loci.tsv exactly (6360/6360). So any change between arms is the filter
and nothing else.

RACE SAFETY. A full-corpus reference-sensitivity run may hold
~/research/datasets/n2-refseq-reference. That directory is read ONLY here, and
only for a partner whose `{partner}.json` sentinel exists with version 2 - that
file is written last by the builder, so its presence means the .faa and .ffn are
complete. Nothing else is shared.

DISK. Start floor 10 GB, abort floor 8 GB, read from /mnt/c and never from
`df /`, which reports WSL's virtual capacity. Measured footprint for this
analysis: analysis_rbh_footprint.py, results/rbh_footprint.json.

    python3 analysis_rbh_orthology.py [--manifest M] [--out-json J]
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
OUT = R / "rbh"
D = pathlib.Path.home() / "research/datasets"
REFDIR = D / "n2-refseq-reference"          # READ ONLY here
QDIR = D / "n2-query-proteomes"             # built here, outside the repo

MIN_PROT_AA = 20
N_NEIGHBOURS = 5
SYNTENY_WINDOW = 10
THREADS = 4
FLOOR_START_GB, FLOOR_ABORT_GB = 10, 8


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


# ----------------------------------------------------------------- partner side
def partner_reference(partner):
    """(faa, ffn) from the shared cache, only if its sentinel says complete."""
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


CDS_COLS = ["id", "prot", "contig", "start", "end", "reclen", "recidx"]


def cds_table(acc):
    """-> [(id, protein, gbff_contig, start, end, record_length, record_index)].

    Cached as a TSV, because both the reverse database and the synteny gene
    order need it and GenBank parsing is the slow step in this script.

    record_length and record_index are carried because the GenBank record ids
    are NOT the contig names the locus table uses - the loci are named in
    Pseudofinder's pf_contigs.fasta namespace. Reconciling the two is the
    project's existing problem (analysis_dnds.contig_key) and the first version
    of this script ignored it, which silently disabled both the annotated-CDS
    acceptance clause and the whole synteny test. Length, with record order as
    a fall-back, is the same reconciliation the project already uses.
    """
    QDIR.mkdir(parents=True, exist_ok=True)
    cache = QDIR / f"{acc}.v2.tsv"
    if cache.exists():
        with cache.open() as fh:
            return [(r["id"], r["prot"], r["contig"], int(r["start"]),
                     int(r["end"]), int(r["reclen"]), int(r["recidx"]))
                    for r in csv.DictReader(fh, delimiter="\t")]
    gb = gbff_for(acc)
    if gb is None:
        return []
    rows = []
    for idx, rec in enumerate(SeqIO.parse(str(gb), "genbank")):
        reclen = len(rec.seq)
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
                         int(ft.location.end), reclen, idx))
    tmp = cache.with_suffix(".tsv.part")
    with tmp.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(CDS_COLS)
        w.writerows(rows)
    tmp.replace(cache)          # atomic, so a partial file is never read
    return rows


def contig_translation(contigs, cds):
    """GenBank record id -> the contig name the locus table uses.

    By unique sequence length first, then by record order, which is what
    analysis_dnds.contig_key does for the caller annotations.
    """
    by_len = collections.defaultdict(list)
    for name, seq in contigs.items():
        by_len[len(seq)].append(name)
    order = list(contigs)
    out, how = {}, collections.Counter()
    for _, _, rid, _, _, reclen, recidx in cds:
        if rid in out:
            continue
        cand = by_len.get(reclen, [])
        if len(cand) == 1:
            out[rid], tag = cand[0], "length"
        elif recidx < len(order):
            out[rid], tag = order[recidx], "order"
        else:
            out[rid], tag = None, "unmapped"
        how[tag] += 1
    return out, dict(how)


def partner_order(partner):
    """-> {locus_tag: (gbff contig, index within that contig by coordinate)}.

    Both sides of the synteny test live in the partner's own GenBank namespace,
    so no reconciliation is needed here.
    """
    rows = cds_table(partner)
    by = collections.defaultdict(list)
    for ident, _, contig, start, _, _, _ in rows:
        by[contig].append((start, ident))
    out = {}
    for contig, v in by.items():
        for i, (_, ident) in enumerate(sorted(v)):
            out[ident] = (contig, i)
    return out


# ----------------------------------------------------------------- searches
def blastp(query_fa, db_faa, workdir):
    """-> {query: subject} single best hit."""
    workdir.mkdir(parents=True, exist_ok=True)
    dmnd = workdir / "db.dmnd"
    subprocess.run(["diamond", "makedb", "--in", str(db_faa), "-d", str(dmnd),
                    "--quiet", "--threads", str(THREADS)], check=True)
    out = workdir / "hits.tsv"
    subprocess.run(["diamond", "blastp", "-q", str(query_fa), "-d", str(dmnd),
                    "-o", str(out), "--quiet", "--threads", str(THREADS),
                    "--max-target-seqs", "1", "--evalue", "1e-5",
                    "--outfmt", "6", "qseqid", "sseqid", "bitscore"],
                   check=True)
    best = {}
    for line in out.read_text().splitlines():
        f = line.split("\t")
        if len(f) >= 2 and f[0] not in best:
            best[f[0]] = f[1]
    return best


def translate_hit_span(nt, qs, qe):
    """The protein the forward blastx actually aligned, in that hit's frame.

    Mirrors analysis_dnds.codon_pairs: slice on min/max of the query bounds and
    reverse-complement when the hit runs backwards. Stops become X so the
    sequence is a legal DIAMOND database entry.
    """
    span = nt[min(qs, qe) - 1:max(qs, qe)]
    if qs > qe:
        span = str(Seq(span).reverse_complement())
    span = span[:len(span) - len(span) % 3]
    if len(span) < MIN_PROT_AA * 3:
        return None
    return str(Seq(span).translate()).replace("*", "X")


def clean(v):
    """A recorded field, or None. NaN is truthy, so it must not be tested raw."""
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return None
    return str(v)


def overlaps_half(a0, a1, b0, b1):
    """50 % reciprocal overlap, the project's own locus criterion."""
    inter = max(0, min(a1, b1) - max(a0, b0) + 1)
    if inter <= 0:
        return False
    return inter >= 0.5 * (a1 - a0 + 1) and inter >= 0.5 * (b1 - b0 + 1)


# ----------------------------------------------------------------- per genome
def process(acc, partner, sub, strands, contigs, syn_cache):
    faa, _ = partner_reference(partner)
    if faa is None:
        return None, "partner reference incomplete or absent"
    pprot = X.read_fasta(faa)
    if not pprot:
        return None, "partner reference unreadable"

    # ---- query nucleotides for each locus
    queries, meta = {}, {}
    for n, r in enumerate(sub.itertuples()):
        st = strands.get((r.accession, r.contig, r.start, r.end))
        if st is None or r.contig not in contigs:
            continue
        if r.end > len(contigs[r.contig]):
            continue
        seq = X.sub(contigs[r.contig], r.start, r.end, st)
        if not seq:
            continue
        k = f"L{n}"
        queries[k], meta[k] = seq, r
    if not queries:
        return None, "no locus sequence recoverable"

    cds = cds_table(acc)
    if not cds:
        return None, "no query proteome"

    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        # ---- forward, re-run rather than read, as a control
        qf = td / "loci.fna"
        with qf.open("w") as fh:
            for k, v in queries.items():
                fh.write(f">{k}\n{v}\n")
        fwd = X.best_hits(qf, faa)
        if not fwd:
            return None, "forward search returned nothing"

        # ---- reverse database: annotated proteome + the loci themselves
        cand = {}
        for k, (sid, qs, qe) in fwd.items():
            pr = translate_hit_span(queries[k], qs, qe)
            if pr:
                cand[k] = pr
        revdb = td / "revdb.faa"
        with revdb.open("w") as fh:
            for ident, prot, *_ in cds:
                fh.write(f">CDS|{ident}\n{prot}\n")
            for k, pr in cand.items():
                fh.write(f">LOC|{k}\n{pr}\n")

        subjects = sorted({sid for sid, _, _ in fwd.values()})
        revq = td / "revq.faa"
        n_rq = 0
        with revq.open("w") as fh:
            for sid in subjects:
                if sid in pprot:
                    fh.write(f">{sid}\n{pprot[sid]}\n")
                    n_rq += 1
        if not n_rq:
            return None, "no subject protein sequences"
        rev = blastp(revq, revdb, td / "rev")

        # ---- synteny: whole query proteome against the partner proteome
        qpf = td / "qprot.faa"
        with qpf.open("w") as fh:
            for ident, prot, *_ in cds:
                fh.write(f">{ident}\n{prot}\n")
        q2p = blastp(qpf, faa, td / "syn")

    # ---- gene order on both sides
    if partner not in syn_cache:
        syn_cache[partner] = partner_order(partner)
    porder = syn_cache[partner]
    xlat, how = contig_translation(contigs, cds)
    by_contig = collections.defaultdict(list)
    cds_iv = {}
    for ident, _, rid, start, end, _, _ in cds:
        contig = xlat.get(rid)
        if contig is None:
            continue
        by_contig[contig].append((start, ident))
        cds_iv[ident] = (contig, start, end)
    for v in by_contig.values():
        v.sort()

    # ---- decide each locus
    rows, n_fwd_agree, n_fwd_cmp = [], 0, 0
    for k, r in meta.items():
        f = fwd.get(k)
        rec = {"accession": acc, "genus": r.genus, "partner": partner,
               "contig": r.contig, "start": r.start, "end": r.end,
               "kind": r.kind, "omega": r.control,
               "forward_subject": f[0] if f else None,
               "recorded_subject": clean(getattr(r, "treatment_subject", None)),
               "rbh": 0, "reverse_best": None, "accept_reason": None,
               "synteny_votes": None}
        if f and rec["recorded_subject"]:
            n_fwd_cmp += 1
            n_fwd_agree += int(f[0] == rec["recorded_subject"])
        if not f:
            rec["accept_reason"] = "no_forward_hit"
            rows.append(rec)
            continue
        best = rev.get(f[0])
        rec["reverse_best"] = best
        if best == f"LOC|{k}":
            rec["rbh"], rec["accept_reason"] = 1, "own_translated_span"
        elif best and best.startswith("CDS|"):
            ident = best[4:]
            iv = cds_iv.get(ident)
            if iv and iv[0] == r.contig and overlaps_half(
                    r.start, r.end, iv[1], iv[2]):
                rec["rbh"], rec["accept_reason"] = 1, "overlapping_annotated_cds"
            else:
                rec["accept_reason"] = "reverse_hit_elsewhere"
        else:
            rec["accept_reason"] = "no_reverse_hit"

        if rec["rbh"]:
            tgt = porder.get(f[0])
            votes = 0
            if tgt:
                near = by_contig.get(r.contig, [])
                mid = r.start
                order = sorted(near, key=lambda t: abs(t[0] - mid))
                taken = 0
                for start, ident in order:
                    iv = cds_iv[ident]
                    if overlaps_half(r.start, r.end, iv[1], iv[2]):
                        continue            # the locus itself, not a neighbour
                    pt = porder.get(q2p.get(ident, ""))
                    taken += 1
                    if pt and pt[0] == tgt[0] and abs(pt[1] - tgt[1]) <= SYNTENY_WINDOW:
                        votes += 1
                    if taken >= N_NEIGHBOURS:
                        break
            rec["synteny_votes"] = votes
        rows.append(rec)

    return (rows, n_fwd_agree, n_fwd_cmp, how), None


# ----------------------------------------------------------------- reporting
def arm(rows, col, genus, label):
    gaps = S.per_genome_gaps(rows, col)
    out = {"label": label, **S.g3_from(gaps)}
    lo, hi = S.cluster_ci(rows, col, genus, unit="genus")
    out["genus_cluster_ci95"] = [lo, hi]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-locus",
                    default="results/reference_sensitivity/per_locus.tsv")
    ap.add_argument("--out-json", default="results/rbh_orthology.json")
    ap.add_argument("--targets", type=int, default=0)
    a = ap.parse_args()

    free0 = guard("start")
    print(f"host C: free {free0} GB (from /mnt/c)")
    if free0 < FLOOR_START_GB:
        sys.exit(f"refusing to start: {free0} < {FLOOR_START_GB} GB")
    OUT.mkdir(parents=True, exist_ok=True)

    src = P / a.per_locus
    raw = src.read_bytes()
    print(f"input {src.name}  sha256 {hashlib.sha256(raw).hexdigest()[:16]}…  "
          f"{len(raw)} bytes")
    df = pd.read_csv(src, sep="\t")
    df = df[df.control.notna()]
    genus = dict(zip(df.accession, df.genus))
    pairs = df[["accession", "partner"]].drop_duplicates().values.tolist()
    if a.targets:
        pairs = pairs[:a.targets]
    print(f"{len(df)} scored loci over {len(pairs)} (query, partner) pairs")

    accs = [p[0] for p in pairs]
    print("\nrecovering locus strands from a 0.50 partition rebuild")
    strands = S.locus_strands(accs)
    print(f"  {len(strands)} representative intervals recovered")

    all_rows, skipped, fa, fc, syn_cache = [], {}, 0, 0, {}
    xlat_how = collections.Counter()
    for i, (acc, partner) in enumerate(pairs, 1):
        contigs = X.genome_seqs(acc)
        if not contigs:
            skipped[acc] = "no contigs"
            continue
        res, why = process(acc, partner, df[df.accession == acc], strands,
                           contigs, syn_cache)
        if res is None:
            skipped[acc] = why
            print(f"  [{i}/{len(pairs)}] {acc}: SKIPPED — {why}", flush=True)
            continue
        rows, agree, cmp_, how = res
        fa, fc = fa + agree, fc + cmp_
        xlat_how.update(how)
        all_rows.extend(rows)
        n_rbh = sum(r["rbh"] for r in rows)
        print(f"  [{i}/{len(pairs)}] {acc}: {len(rows)} loci, "
              f"{n_rbh} RBH-confirmed ({100 * n_rbh / len(rows):.1f}%)",
              flush=True)
        guard(acc)

    if not all_rows:
        sys.exit("nothing processed")

    out = pd.DataFrame(all_rows)
    out["omega_rbh"] = np.where(out.rbh == 1, out.omega, np.nan)
    votes = pd.to_numeric(out.synteny_votes, errors="coerce").fillna(0)
    out["omega_syn1"] = np.where((out.rbh == 1) & (votes >= 1), out.omega, np.nan)
    out["omega_syn2"] = np.where((out.rbh == 1) & (votes >= 2), out.omega, np.nan)
    out.to_csv(OUT / "per_locus_rbh.tsv", sep="\t", index=False)
    print(f"\nwrote {len(out)} locus rows to {OUT / 'per_locus_rbh.tsv'}")

    # ---- control 1: the re-run forward search must reproduce what was recorded
    fwd_agree = fa / fc if fc else None
    print("\nCONTROL 1 — re-run forward search vs recorded treatment subject")
    print(f"  {fa} of {fc} identical"
          + (f"  ({100 * fwd_agree:.1f}%)" if fwd_agree is not None else ""))
    print("\nCONTROL 2 — GenBank record to locus-contig reconciliation")
    tot = sum(xlat_how.values()) or 1
    for k, v in sorted(xlat_how.items()):
        print(f"  by {k:9s} {v:6d}  ({100 * v / tot:5.1f}%)")

    # ---- attrition, separately per side
    print("\nATTRITION, by locus class (the asymmetry that would rig the test)")
    attr = {}
    for kind, g in out.groupby("kind"):
        attr[kind] = {
            "n": int(len(g)),
            "rbh": int(g.rbh.sum()),
            "rbh_frac": float(g.rbh.mean()),
            "syn1": int(g.omega_syn1.notna().sum()),
            "syn1_frac": float(g.omega_syn1.notna().mean()),
            "syn2": int(g.omega_syn2.notna().sum()),
            "syn2_frac": float(g.omega_syn2.notna().mean())}
        print(f"  {kind:20s} n {len(g):6d}  RBH {g.rbh.sum():6d} "
              f"({100 * g.rbh.mean():5.1f}%)  "
              f"+syn>=1 {g.omega_syn1.notna().sum():6d} "
              f"({100 * g.omega_syn1.notna().mean():5.1f}%)  "
              f"+syn>=2 {g.omega_syn2.notna().sum():6d} "
              f"({100 * g.omega_syn2.notna().mean():5.1f}%)")
    da, ba = attr.get("discordant", {}), attr.get("concordant_intact", {})
    gap_rbh = (da.get("rbh_frac"), ba.get("rbh_frac"))
    if None not in gap_rbh:
        print(f"  retention difference, discordant minus baseline: "
              f"{100 * (gap_rbh[0] - gap_rbh[1]):+.1f} percentage points")

    print("\nOUTCOME PER LOCUS (first two are acceptances)")
    for reason, n in out.accept_reason.value_counts(dropna=False).items():
        print(f"  {str(reason):26s} {n:6d}  ({100 * n / len(out):5.1f}%)")

    rows_d = out.to_dict("records")
    arms = [arm(rows_d, "omega", genus, "all scored loci (= published)"),
            arm(rows_d, "omega_rbh", genus, "RBH-confirmed"),
            arm(rows_d, "omega_syn1", genus, "RBH + synteny >= 1"),
            arm(rows_d, "omega_syn2", genus, "RBH + synteny >= 2")]
    print("\nG3 BY ARM (omega never recomputed; only the locus set changes)")
    for A in arms:
        g = A.get("median_gap")
        ci = A.get("genus_cluster_ci95") or [None, None]
        if g is None:
            print(f"  {A['label']:32s} {A.get('note')}")
            continue
        line = (f"  {A['label']:32s} {g:+.4f}  "
                f"{A['n_positive']}/{A['n_tested']}  p {A['p']:.3g}")
        if ci[0] is not None:
            line += f"  genus CI [{ci[0]:+.4f}, {ci[1]:+.4f}]"
        print(line)

    # ---- paired: does the filter move the per-genome gap?
    g_all = S.per_genome_gaps(rows_d, "omega")
    g_rbh = S.per_genome_gaps(rows_d, "omega_rbh")
    shared = sorted(set(g_all) & set(g_rbh))
    deltas = [g_rbh[k] - g_all[k] for k in shared]
    paired = {"n_genomes": len(shared)}
    if len(deltas) >= 5:
        pos = sum(1 for d in deltas if d > 0)
        nz = [d for d in deltas if d != 0]
        paired.update({
            "median_delta": float(np.median(deltas)),
            "n_positive": pos,
            "sign_test_p": float(stats.binomtest(
                sum(1 for d in nz if d > 0), len(nz), 0.5).pvalue) if nz else None})
        print(f"\nPAIRED per-genome gap change, RBH minus all, "
              f"{len(shared)} genomes")
        print(f"  median {paired['median_delta']:+.5f}  "
              f"{pos}/{len(shared)} positive  "
              f"sign test p {paired['sign_test_p']}")

    res = {"note": "POST HOC; reciprocal best hit plus synteny as a filter",
           "input_per_locus": str(src.relative_to(P)),
           "input_sha256": hashlib.sha256(raw).hexdigest(),
           "n_pairs": len(pairs), "n_loci": int(len(out)),
           "min_prot_aa": MIN_PROT_AA, "n_neighbours": N_NEIGHBOURS,
           "synteny_window": SYNTENY_WINDOW,
           "control_forward_agreement": fwd_agree,
           "control_forward_n": fc,
           "contig_reconciliation": dict(xlat_how),
           "attrition": attr, "arms": arms, "paired_rbh_vs_all": paired,
           "rejection_reasons": {str(k): int(v) for k, v in
                                 out.accept_reason.value_counts(
                                     dropna=False).items()},
           "skipped": skipped,
           "host_free_gb_start": free0, "host_free_gb_end": host_free_gb()}
    (P / a.out_json).write_text(json.dumps(res, indent=1, default=str,
                                           sort_keys=True))
    print(f"\nwrote {a.out_json}")


if __name__ == "__main__":
    main()
