#!/usr/bin/env python3
"""Does G3 depend on the orthologue reference being Pseudofinder's intact set?

POST HOC. Not registered. No frozen plan and no validated output is modified;
everything is written under results/reference_sensitivity/.

THE PROBLEM. Every dN/dS value in this project - PGAP's, Bakta's and DFAST's
included - is computed against the PARTNER genome's pf_intact.ffn/.faa, which
is Pseudofinder's INTACT set for that partner. Genes Pseudofinder called
pseudogenes in the partner are absent from the reference and can never serve as
an orthologue. Across the 58 partners used that is 15,777 of 143,129 genes,
11.0 % overall, median 6.7 %, and up to 35.0 % for one partner
(results/pf_residual_dependence.json). That is a Pseudofinder filter sitting
under the measuring instrument for all four callers.

THE INDEPENDENT ALTERNATIVE. The partner's own RefSeq GenBank annotation. Every
CDS carrying no /pseudo qualifier, with its /translation for the protein search
and the extracted nucleotides for the codon alignment. This is PGAP-derived
rather than caller-neutral - stated as a limitation, not hidden - but it is
independent of Pseudofinder and, crucially, it CONTAINS the genes Pseudofinder
withheld.

THREE ARMS PER LOCUS, so that a difference can be attributed:
  published   the omega already in loci.tsv
  control     recomputed here against pf_intact - must reproduce `published`,
              or nothing else in this script means anything
  treatment   recomputed here against the RefSeq-derived reference

No assumption is made about which way including more genes should push G3. The
direction and magnitude are measured.

    python3 analysis_reference_sensitivity.py [--targets N] [--baseline-cap N]
"""
from __future__ import annotations

import argparse
import ast
import collections
import csv
import json
import pathlib
import tempfile
import types

import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy import stats

P = pathlib.Path(__file__).resolve().parent
R = P / "results"
OUT = R / "reference_sensitivity"
D = pathlib.Path.home() / "research/datasets"
REFDIR = D / "n2-refseq-reference"      # built here, outside the repo
MIN_LOCI = 3


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


# ------------------------------------------------------- the new reference
def build_refseq_reference(partner: str):
    """-> (faa, ffn) of every non-pseudo CDS in the partner's own annotation.

    Protein from /translation, nucleotides extracted from the record, keyed
    identically in both files so a DIAMOND hit can be looked up for the codon
    alignment. Cached, so a partner shared by several targets is built once.
    """
    REFDIR.mkdir(parents=True, exist_ok=True)
    faa, ffn = REFDIR / f"{partner}.faa", REFDIR / f"{partner}.ffn"
    stats_f = REFDIR / f"{partner}.json"
    # v2 keys on locus_tag. Pseudofinder's pf_intact headers are locus_tags
    # (">ABU_RS00005 NC_009850 [1:1317](+)"), so keying on protein_id made the
    # orthologue-identity comparison read 0 % agreement for trivial reasons.
    VERSION = 2
    if faa.exists() and ffn.exists() and stats_f.exists():
        st = json.loads(stats_f.read_text())
        if st.get("version") == VERSION:
            return faa, ffn, st
    gb = X.D / f"n2-genomes-expansion/{partner}/genomic.gbff"
    if not gb.exists():
        gb = next((X.D / "n2-genomes/extracted").rglob(
            f"*{partner}*/genomic.gbff"), None)
    if gb is None or not gb.exists():
        return None, None, None
    n_cds = n_pseudo = n_kept = n_no_tr = 0
    with faa.open("w") as fa, ffn.open("w") as fn:
        for rec in SeqIO.parse(str(gb), "genbank"):
            for ft in rec.features:
                if ft.type != "CDS":
                    continue
                n_cds += 1
                q = ft.qualifiers
                if "pseudo" in q or "pseudogene" in q:
                    n_pseudo += 1
                    continue
                tr = q.get("translation", [None])[0]
                if not tr:
                    n_no_tr += 1
                    continue
                # locus_tag FIRST, to match Pseudofinder's own key
                ident = (q.get("locus_tag") or q.get("protein_id")
                         or [f"{rec.id}_{int(ft.location.start)}"])[0]
                try:
                    nt = str(ft.extract(rec).seq).upper()
                except Exception:
                    continue
                if len(nt) < 90:
                    continue
                fa.write(f">{ident}\n{tr}\n")
                fn.write(f">{ident}\n{nt}\n")
                n_kept += 1
    st = {"version": VERSION, "partner": partner, "cds_total": n_cds,
          "cds_pseudo": n_pseudo, "cds_no_translation": n_no_tr,
          "reference_genes": n_kept}
    stats_f.write_text(json.dumps(st, indent=2) + "\n")
    return faa, ffn, st


# ------------------------------------------------------- locus strands
def locus_strands(accs):
    """loci.tsv carries no strand, so rebuild the partition to recover it.

    analysis_overlap_sensitivity.py showed the rebuild at RECIPROCAL = 0.50
    reproduces the validated locus set exactly (103 224 of 103 224), so the
    strand of the representative interval can be recovered by matching on
    (accession, contig, start, end).
    """
    L.RECIPROCAL = 0.50
    out = {}
    for acc in accs:
        contigs = X.genome_seqs(acc)
        if not contigs:
            continue
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
        if not ivs:
            continue
        for comp in L.components(ivs):
            if len(comp) < 2:
                continue
            longest = max(comp, key=lambda i: ivs[i][2] - ivs[i][1])
            c, s, e, st = ivs[longest][:4]
            out[(acc, c, s, e)] = st
    return out


# ------------------------------------------------------- omega
def omegas(queries, faa, ffn):
    """-> {name: (omega, codons, subject)} through the project's own pipeline."""
    pcds = X.read_fasta(ffn)
    out = {}
    with tempfile.TemporaryDirectory() as td:
        qf = pathlib.Path(td) / "q.fna"
        with qf.open("w") as fh:
            for k, v in queries.items():
                fh.write(f">{k}\n{v}\n")
        hits = X.best_hits(qf, faa)
    for name, (sid, qs, qe) in hits.items():
        pnt = pcds.get(sid)
        if pnt is None:
            continue
        cp = X.codon_pairs(queries[name], pnt, qs, qe)
        if cp is None:
            continue
        dn, ds, used = X.ng86(*cp)
        if dn is None or ds is None or ds <= 0:
            continue
        out[name] = (dn / ds, used, sid)
    return out


def g3_from(gaps):
    nz = [v for v in gaps.values() if v != 0.0]
    if len(nz) < 5:
        return {"n_genomes": len(gaps), "p": None, "median_gap": None,
                "note": f"only {len(nz)} non-tied genomes"}
    pos = sum(1 for v in nz if v > 0)
    return {"n_genomes": len(gaps), "n_tested": len(nz),
            "median_gap": float(np.median(list(gaps.values()))),
            "n_positive": pos,
            "p": float(stats.binomtest(pos, len(nz), 0.5).pvalue)}


def per_genome_gaps(rows, col):
    by = collections.defaultdict(lambda: {"discordant": [], "concordant_intact": []})
    for r in rows:
        v = r.get(col)
        if v is not None and not (isinstance(v, float) and np.isnan(v)):
            by[r["accession"]][r["kind"]].append(v)
    out = {}
    for acc, d in by.items():
        a, b = d["discordant"], d["concordant_intact"]
        if len(a) >= MIN_LOCI and len(b) >= MIN_LOCI:
            out[acc] = float(np.median(a) - np.median(b))
    return out


def cluster_ci(rows, col, genus, unit="genus", B=2000, seed=20261108):
    gaps = per_genome_gaps(rows, col)
    buckets = collections.defaultdict(list)
    for acc, v in gaps.items():
        buckets[genus.get(acc, "?") if unit == "genus" else acc].append(v)
    keys = list(buckets)
    if len(keys) < 3:
        return None, None
    g = np.random.default_rng(seed)
    meds = []
    for _ in range(B):
        v = [x for k in g.choice(keys, len(keys), replace=True) for x in buckets[k]]
        meds.append(float(np.median(v)))
    return float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", type=int, default=0)
    ap.add_argument("--baseline-cap", type=int, default=400)
    ap.add_argument("--manifest", default="manifests/refsens_pilot.tsv")
    ap.add_argument("--out-json", default="results/reference_sensitivity.json")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    pilot = pd.read_csv(P / a.manifest, sep="\t")
    if a.targets:
        pilot = pilot.head(a.targets)
    accs = list(pilot.accession)
    genus = dict(zip(pilot.accession, pilot.genus))
    partner = dict(zip(pilot.accession, pilot.partner))
    print(f"pilot: {len(accs)} targets, {pilot.genus.nunique()} genera, "
          f"{pilot.stratum.nunique()} strata, exclusion "
          f"{pilot.frac_excluded.min():.3f}-{pilot.frac_excluded.max():.3f}")

    loci = pd.read_csv(R / "loci.tsv", sep="\t")
    loci = loci[loci.accession.isin(accs)]
    print(f"loci in these targets: {len(loci)}")

    print("\nrecovering locus strands from a 0.50 partition rebuild")
    strands = locus_strands(accs)
    print(f"  {len(strands)} representative intervals recovered")

    print("\nbuilding the RefSeq-derived reference for each partner")
    refs, refstats = {}, []
    for pt in sorted(set(partner.values())):
        faa, ffn, st = build_refseq_reference(pt)
        if faa is None:
            print(f"  {pt}: NO GenBank; skipped")
            continue
        refs[pt] = (faa, ffn)
        pf_n = sum(1 for _ in open(D / f"n2-pseudofinder/{pt}/pf_intact.ffn")
                   if _.startswith(">"))
        st["pf_intact_genes"] = pf_n
        st["extra_vs_pf"] = st["reference_genes"] - pf_n
        refstats.append(st)
        print(f"  {pt}: refseq {st['reference_genes']} genes vs pf_intact "
              f"{pf_n}  ({st['reference_genes'] - pf_n:+d}), "
              f"{st['cds_pseudo']} CDS excluded as pseudo")
    pd.DataFrame(refstats).to_csv(OUT / "reference_build.tsv", sep="\t",
                                 index=False)

    rows = []
    for i, acc in enumerate(accs, 1):
        pt = partner[acc]
        if pt not in refs:
            continue
        contigs = X.genome_seqs(acc)
        if not contigs:
            continue
        sub = loci[loci.accession == acc]
        dis = sub[sub.kind == "discordant"]
        con = sub[sub.kind == "concordant_intact"]
        if len(con) > a.baseline_cap:          # evenly spaced, no RNG
            step = len(con) / a.baseline_cap
            con = con.iloc[[int(j * step) for j in range(a.baseline_cap)]]
        use = pd.concat([dis, con])
        queries, meta = {}, {}
        for n, r in enumerate(use.itertuples()):
            st = strands.get((r.accession, r.contig, r.start, r.end))
            if st is None:
                continue
            if r.contig not in contigs or r.end > len(contigs[r.contig]):
                continue
            seq = X.sub(contigs[r.contig], r.start, r.end, st)
            if not seq:
                continue
            k = f"L{n}"
            queries[k] = seq
            meta[k] = r
        if not queries:
            continue
        pf_faa = D / f"n2-pseudofinder/{pt}/pf_intact.faa"
        pf_ffn = D / f"n2-pseudofinder/{pt}/pf_intact.ffn"
        o_ctl = omegas(queries, pf_faa, pf_ffn)
        o_new = omegas(queries, refs[pt][0], refs[pt][1])
        for k, r in meta.items():
            c = o_ctl.get(k)
            n_ = o_new.get(k)
            rows.append({
                "accession": acc, "genus": genus[acc], "partner": pt,
                "contig": r.contig, "start": r.start, "end": r.end,
                "kind": r.kind, "published": float(r.omega),
                "control": c[0] if c else None,
                "treatment": n_[0] if n_ else None,
                "control_subject": c[2] if c else None,
                "treatment_subject": n_[2] if n_ else None,
                "same_subject": (int(c[2] == n_[2]) if c and n_ else None)})
        print(f"  [{i}/{len(accs)}] {acc}: {len(queries)} queries, "
              f"control {len(o_ctl)}, treatment {len(o_new)}", flush=True)

    tag = "all" if "all" in a.manifest else "pilot"
    df = pd.DataFrame(rows)
    df.to_csv(OUT / f"per_locus_{tag}.tsv", sep="\t", index=False)
    print(f"\nwrote {len(df)} locus rows")

    # ---------------------------------------------------------- control
    both = df[df.control.notna()]
    agree = both[np.isclose(both.published, both.control, atol=1e-6, rtol=0)]
    print("\nCONTROL: recomputing against pf_intact must reproduce loci.tsv")
    print(f"  {len(agree)} of {len(both)} recomputed omegas identical "
          f"({100 * len(agree) / max(len(both), 1):.1f}%)")
    ctl_ok = len(agree) / max(len(both), 1) >= 0.95

    # ---------------------------------------------------------- coverage
    n = len(df)
    cov = {
        "n_loci": n,
        "control_scored": int(df.control.notna().sum()),
        "treatment_scored": int(df.treatment.notna().sum()),
        "both_scored": int((df.control.notna() & df.treatment.notna()).sum()),
        "treatment_only": int((df.control.isna() & df.treatment.notna()).sum()),
        "control_only": int((df.control.notna() & df.treatment.isna()).sum()),
    }
    print("\nCOVERAGE and MISSINGNESS")
    for k, v in cov.items():
        print(f"  {k:20s} {v:6d}" + (f"  ({100 * v / n:.1f}%)" if k != "n_loci" else ""))
    same = df[df.same_subject.notna()]
    frac_same = float(same.same_subject.mean()) if len(same) else None
    print(f"  orthologue identity agrees on {100 * frac_same:.1f}% of loci "
          f"scored by both" if frac_same is not None else "  no overlap")

    # ---------------------------------------------------------- effect
    print("\nG3 PER GENOME, by reference")
    res = {}
    for col in ("published", "control", "treatment"):
        gaps = per_genome_gaps(df.to_dict("records"), col)
        r = g3_from(gaps)
        lo, hi = cluster_ci(df.to_dict("records"), col, genus, "genus")
        r["genus_cluster_ci95"] = [lo, hi]
        res[col] = r
        if r.get("p") is None:
            print(f"  {col:10s} NOT TESTABLE ({r.get('note')})")
        else:
            ci = (f"[{lo:+.4f}, {hi:+.4f}]" if lo is not None else "n/a")
            print(f"  {col:10s} {r['median_gap']:+.4f}  "
                  f"{r['n_positive']}/{r['n_tested']}  p {r['p']:.3g}  "
                  f"genus CI {ci}")

    # paired, per locus, on loci scored by both
    pair = df[df.control.notna() & df.treatment.notna()].copy()
    pair["delta"] = pair.treatment - pair.control
    print(f"\nPAIRED per-locus change, {len(pair)} loci scored by both")
    for kind in ("discordant", "concordant_intact"):
        s = pair[pair.kind == kind]
        if len(s) < 5:
            continue
        print(f"  {kind:20s} n {len(s):5d}  median delta "
              f"{np.median(s.delta):+.5f}  mean {np.mean(s.delta):+.5f}")
    # and the per-genome gap change, paired by genome
    gp_c = per_genome_gaps(pair.to_dict("records"), "control")
    gp_t = per_genome_gaps(pair.to_dict("records"), "treatment")
    common = sorted(set(gp_c) & set(gp_t))
    d = [gp_t[k] - gp_c[k] for k in common]
    pp = None
    if len(d) >= 5:
        nz = [x for x in d if x != 0]
        pp = float(stats.binomtest(sum(1 for x in nz if x > 0), len(nz),
                                   0.5).pvalue) if len(nz) >= 5 else None
    print(f"  per-genome gap change: {len(common)} genomes, median "
          f"{np.median(d):+.5f}" + (f", sign test p {pp:.3g}" if pp else ""))

    out = {"note": "POST HOC; independent orthologue reference sensitivity",
           "pilot": {"n_targets": len(accs), "n_genera": int(pilot.genus.nunique()),
                     "n_strata": int(pilot.stratum.nunique()),
                     "exclusion_min": float(pilot.frac_excluded.min()),
                     "exclusion_max": float(pilot.frac_excluded.max())},
           "control_reproduces_published": bool(ctl_ok),
           "control_identical_frac": len(agree) / max(len(both), 1),
           "coverage": cov, "orthologue_identity_agreement": frac_same,
           "g3": res,
           "paired": {"n_loci": int(len(pair)),
                      "median_delta_discordant": float(np.median(
                          pair[pair.kind == "discordant"].delta))
                      if (pair.kind == "discordant").sum() else None,
                      "median_delta_baseline": float(np.median(
                          pair[pair.kind == "concordant_intact"].delta))
                      if (pair.kind == "concordant_intact").sum() else None,
                      "n_genomes": len(common),
                      "median_gap_change": float(np.median(d)) if d else None,
                      "sign_test_p": pp}}
    (P / a.out_json).write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    # name the files actually written: the pilot and the expansion must never
    # be confused for each other, and a hard-coded message once said
    # "reference_sensitivity.json" even when the expansion wrote _all.json
    print(f"\nwrote {a.out_json}, {OUT / f'per_locus_{tag}.tsv'} and "
          f"{OUT / 'reference_build.tsv'}")
    print("results/loci.tsv and results/locus.json were NOT touched")


if __name__ == "__main__":
    main()
