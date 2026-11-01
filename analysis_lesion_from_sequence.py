#!/usr/bin/env python3
"""Classify discordant loci as frameshift or internal stop from SEQUENCE alone.

POST HOC. Not registered. Writes only to results/lesion_from_sequence*.

WHY. The lesion-stratified G3 covers 11.0 % of discordant loci because DFAST is
the only caller that publishes a per-call lesion class: PGAP's RefSeq pseudogene
features carry no /note at all, and Pseudofinder's reasons are length- and
fragmentation-based. Pseudofinder is on the pseudogene side of 84.1 % of
discordant loci and is unclassifiable from its own output. A caller-independent
classifier is the only way to lift that ceiling.

THE SIGNAL. DIAMOND `blastx` with `--max-hsps 0` reports every high-scoring
segment of a nucleotide query against a protein. For a frameshifted gene the
alignment breaks into two segments in DIFFERENT reading frames; the frame of an
HSP is `(qstart - 1) % 3`. Validated on synthetic controls before use: a gene
with one base deleted mid-sequence gives two HSPs in frames 0 and 2, and the
same gene undamaged gives a single HSP in frame 0.

An internal stop is read directly from the query's own frame.

FOUR OUTCOMES, kept apart as the brief requires:
  confirmed_frameshift      >=2 substantial HSPs in >=2 frames, no premature stop
  confirmed_internal_stop   premature in-frame stop, single frame
  ambiguous                 both signals, or frame evidence that is marginal
  unclassifiable            no orthologue hit at all

EXTERNAL VALIDATION. DFAST publishes a class for some of these loci. Agreement
between this classifier and DFAST's annotation is reported, which tests the
classifier against an independent label rather than asserting it.

NOT ASSUMED: that a frameshift-aware aligner removes alignment or selection
bias. This classifies lesions; it does not re-estimate omega.

    python3 analysis_lesion_from_sequence.py [--targets N]
"""
from __future__ import annotations

import argparse
import ast
import collections
import csv
import json
import pathlib
import re
import subprocess
import tempfile
import types

import pandas as pd
from Bio.Seq import Seq

P = pathlib.Path(__file__).resolve().parent
R = P / "results"
D = pathlib.Path.home() / "research/datasets"
MIN_HSP_AA = 20          # an HSP must be this long to count as frame evidence
LOC = re.compile(r"^([^:]+):<?(\d+)\.\.>?(\d+)\(([+-])\)")
DFAST_DISRUPTED = {"frameshift", "internal_stop_codon",
                   "frameshift,internal_stop_codon"}


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


L = borrow(P / "analysis_locus.py")
X = L.X
CALLERS = L.CALLERS


def hsps(queries, faa, threads=8):
    """-> {query: [(subject, qstart, qend, aa_len, bits)]} with every HSP."""
    dmnd = X.__dict__.get("DIAMOND") or "diamond"
    out = collections.defaultdict(list)
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        db = td / "db.dmnd"
        qf = td / "q.fna"
        with qf.open("w") as fh:
            for k, v in queries.items():
                fh.write(f">{k}\n{v}\n")
        subprocess.run(["diamond", "makedb", "--in", str(faa), "-d", str(db),
                        "--quiet", "--threads", str(threads)], check=True)
        res = td / "h.tsv"
        subprocess.run(["diamond", "blastx", "-q", str(qf), "-d", str(db),
                        "-o", str(res), "--quiet", "--threads", str(threads),
                        "--max-target-seqs", "1", "--max-hsps", "0",
                        "--evalue", "1e-5", "--outfmt", "6", "qseqid",
                        "sseqid", "pident", "length", "qstart", "qend",
                        "bitscore"], check=True)
        for line in res.read_text().splitlines():
            f = line.split("\t")
            if len(f) < 7:
                continue
            out[f[0]].append((f[1], int(f[4]), int(f[5]), int(f[3]),
                              float(f[6])))
    return out


def premature_stop(nt):
    """Is there an in-frame stop before the final codon, in the query's frame?"""
    trimmed = nt[:len(nt) - len(nt) % 3]
    if len(trimmed) < 6:
        return False
    prot = str(Seq(trimmed).translate())
    return "*" in prot[:-1]


def classify(nt, hs):
    """-> (class, n_frames, n_hsps, detail)"""
    if not hs:
        return "unclassifiable", 0, 0, "no orthologue hit"
    best = max(hs, key=lambda h: h[4])[0]
    same = [h for h in hs if h[0] == best]
    strong = [h for h in same if h[3] >= MIN_HSP_AA]
    frames = {(h[1] - 1) % 3 for h in strong}
    stop = premature_stop(nt)
    if len(frames) >= 2 and stop:
        return "ambiguous", len(frames), len(strong), "frame change AND premature stop"
    if len(frames) >= 2:
        return "confirmed_frameshift", len(frames), len(strong), \
            f"{len(strong)} HSPs in frames {sorted(frames)}"
    if stop:
        return "confirmed_internal_stop", len(frames), len(strong), \
            "single frame, premature in-frame stop"
    if len({(h[1] - 1) % 3 for h in same}) >= 2:
        return "ambiguous", len({(h[1] - 1) % 3 for h in same}), len(same), \
            f"frame change only in HSPs shorter than {MIN_HSP_AA} aa"
    return "ambiguous", len(frames), len(strong), "no disruption signal found"


def dfast_labels(acc):
    """DFAST's own class per interval, in reference coordinates, for validation."""
    contigs = X.genome_seqs(acc)
    if not contigs:
        return []
    resolve = X.contig_key(contigs)
    reglen = X.gff_region_lengths(D / "n2-dfast" / acc / "genome.gff")
    p = D / "n2-dfast" / acc / "pseudogene_summary.tsv"
    if not p.exists():
        return []
    out, seen = [], {}
    for line in p.read_text(errors="replace").splitlines()[1:]:
        f = line.split("\t")
        if len(f) < 4 or f[3] not in DFAST_DISRUPTED:
            continue
        m = LOC.match(f[2])
        if not m:
            continue
        cg = m.group(1)
        if cg not in seen:
            seen[cg] = len(seen)
        real, _ = resolve(cg, reglen.get(cg), seen[cg])
        if real is None:
            continue
        out.append((real, int(m.group(2)) - 1, int(m.group(3)), f[3]))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", type=int, default=12)
    a = ap.parse_args()

    pilot = pd.read_csv(P / "manifests/refsens_pilot.tsv", sep="\t").head(a.targets)
    loci = pd.read_csv(R / "loci.tsv", sep="\t")
    accs = list(pilot.accession)
    partner = dict(zip(pilot.accession, pilot.partner))
    dis = loci[(loci.kind == "discordant") & loci.accession.isin(accs)]
    print(f"{len(dis)} discordant loci across {len(accs)} pilot genomes")

    L.RECIPROCAL = 0.50
    rows = []
    for i, acc in enumerate(accs, 1):
        sub = dis[dis.accession == acc]
        if not len(sub):
            continue
        contigs = X.genome_seqs(acc)
        if not contigs:
            continue
        # strands, from the 0.50 partition that reproduces the validated set
        resolve = X.contig_key(contigs)
        reglen = {
            "Bakta": X.gff_region_lengths(D / "n2-bakta" / acc / f"{acc}.gff3"),
            "DFAST": X.gff_region_lengths(D / "n2-dfast" / acc / "genome.gff"),
        }
        ivs = []
        for caller in CALLERS:
            seen = {}
            for klass, want in (("pseudo", True), ("intact", False)):
                for (c, s, e, st) in X.calls_for(acc, caller, want):
                    if c not in seen:
                        seen[c] = len(seen)
                    real, _ = resolve(c, reglen.get(caller, {}).get(c), seen[c])
                    if real is None or e > len(contigs[real]) or e - s < 30:
                        continue
                    ivs.append((real, s, e, st, caller, klass))
        strands = {}
        for comp in L.components(ivs):
            if len(comp) < 2:
                continue
            lg = max(comp, key=lambda j: ivs[j][2] - ivs[j][1])
            c, s, e, st = ivs[lg][:4]
            strands[(c, s, e)] = st

        dl = dfast_labels(acc)
        queries, meta = {}, {}
        for n, r in enumerate(sub.itertuples()):
            st = strands.get((r.contig, r.start, r.end))
            if st is None or r.contig not in contigs:
                continue
            seq = X.sub(contigs[r.contig], r.start, r.end, st)
            if not seq:
                continue
            k = f"L{n}"
            queries[k] = seq
            # DFAST's own label for this locus, by >=50% reciprocal overlap
            lab, bov = None, 0
            qlen = r.end - r.start
            for cg, s2, e2, cl in dl:
                if cg != r.contig:
                    continue
                ov = min(r.end, e2) - max(r.start, s2)
                if ov > 0 and ov >= 0.5 * qlen and ov >= 0.5 * (e2 - s2) and ov > bov:
                    lab, bov = cl, ov
            meta[k] = (r, lab)
        if not queries:
            continue
        faa = D / f"n2-pseudofinder/{partner[acc]}/pf_intact.faa"
        h = hsps(queries, faa)
        for k, (r, lab) in meta.items():
            cls, nf, nh, detail = classify(queries[k], h.get(k, []))
            rows.append({"accession": acc, "genus": r.genus,
                         "contig": r.contig, "start": r.start, "end": r.end,
                         "pseudo_callers": r.pseudo_callers,
                         "omega": r.omega, "n_frames": nf, "n_hsps": nh,
                         "sequence_class": cls, "dfast_label": lab,
                         "detail": detail})
        print(f"  [{i}/{len(accs)}] {acc}: {len(queries)} loci classified",
              flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(R / "lesion_from_sequence.tsv", sep="\t", index=False)
    print(f"\n{len(df)} discordant loci classified from sequence")
    dist = df.sequence_class.value_counts()
    for k, v in dist.items():
        print(f"  {k:26s} {v:5d}  ({100 * v / len(df):.1f}%)")

    # ---- external validation against DFAST's own labels
    val = df[df.dfast_label.notna()].copy()
    print(f"\nEXTERNAL VALIDATION against DFAST's published class, "
          f"{len(val)} loci")
    out = {"note": "POST HOC; caller-independent lesion classification",
           "n_classified": int(len(df)),
           "distribution": {str(k): int(v) for k, v in dist.items()},
           "min_hsp_aa": MIN_HSP_AA}
    if len(val):
        val["dfast_has_fs"] = val.dfast_label.str.contains("frameshift")
        val["seq_says_fs"] = val.sequence_class == "confirmed_frameshift"
        val["seq_says_stop"] = val.sequence_class == "confirmed_internal_stop"
        tab = pd.crosstab(val.dfast_label, val.sequence_class)
        print(tab.to_string())
        decided = val[val.sequence_class.isin(
            ["confirmed_frameshift", "confirmed_internal_stop"])]
        if len(decided):
            agree = int((decided.dfast_has_fs == decided.seq_says_fs).sum())
            print(f"\n  on {len(decided)} loci where BOTH decide: "
                  f"{agree} agree ({100 * agree / len(decided):.1f}%)")
            out["validation"] = {
                "n_with_dfast_label": int(len(val)),
                "n_both_decide": int(len(decided)),
                "n_agree": agree,
                "agreement": agree / len(decided)}
        out["crosstab"] = {str(i): {str(c): int(v) for c, v in row.items()}
                           for i, row in tab.iterrows()}
    # ---- coverage gain
    rr = json.loads((R / "review_response.json").read_text())
    n_dis_all = rr["R6_exposure"]["n_discordant"]
    decided_n = int(df.sequence_class.isin(
        ["confirmed_frameshift", "confirmed_internal_stop"]).sum())
    print(f"\nCOVERAGE: {decided_n} of {len(df)} pilot discordant loci decided "
          f"({100 * decided_n / len(df):.1f}%)")
    print(f"  DFAST's published class covered 11.0% of all "
          f"{n_dis_all} discordant loci")
    out["coverage"] = {"pilot_loci": int(len(df)), "decided": decided_n,
                       "frac_decided": decided_n / max(len(df), 1)}
    (R / "lesion_from_sequence.json").write_text(
        json.dumps(out, indent=2, sort_keys=True) + "\n")
    print("\nwrote results/lesion_from_sequence.{tsv,json}")


if __name__ == "__main__":
    main()
