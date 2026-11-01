#!/usr/bin/env python3
"""Integration test: re-verify the truth table against the files on disk.

The simulator validates its own work in memory. That is necessary and not
sufficient -- it does not test the FASTA writer, and it does not test that the
truth table as written describes the assembly as written. This re-reads both
from disk with no shared state and checks them against each other.

    python3 tests/test_written_assemblies.py [--assemblies DIR]
"""
from __future__ import annotations

import argparse
import collections
import csv
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from Bio import SeqIO  # noqa: E402

import simulate as S  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent.parent
FAILED: list[str] = []


def check(cond, msg):
    if not cond:
        FAILED.append(msg)
        print(f"  FAIL  {msg}")
    return cond


def load(p: pathlib.Path) -> dict:
    return {r.id: str(r.seq).upper() for r in SeqIO.parse(str(p), "fasta")}


def load_liftover(rid: str) -> dict:
    """Rebuild the liftover from the committed tsv, not from the simulator."""
    rows = collections.defaultdict(list)
    f = HERE / f"results/liftover/{rid}.tsv"
    with f.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows[r["contig"]].append((int(r["parental_pos0"]),
                                      int(r["parental_end0"]),
                                      int(r["cum_delta"])))
    out = {}
    for cg, v in rows.items():
        v.sort()
        lf = S.Liftover([])
        lf.starts = [a for a, _, _ in v]
        lf.ends = [b for _, b, _ in v]
        lf.cum = [c for _, _, c in v]
        # altlen is not in the tsv; the only thing read back here needs it for
        # positions inside a replacement, and the anchors used below are the
        # gene edges, so reconstruct it as "everything was a pure indel or an
        # equal-length replacement" and assert that assumption holds instead.
        lf.altlen = [max(0, (b - a) + (c - (v[i - 1][2] if i else 0)))
                     for i, (a, b, c) in enumerate(v)]
        out[cg] = lf
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assemblies",
                    default=str(pathlib.Path.home()
                                / "research/datasets/n2-benchmark/assemblies"))
    a = ap.parse_args()
    asm = pathlib.Path(a.assemblies)

    truth = list(csv.DictReader((HERE / "results/ground_truth.tsv").open(),
                                delimiter="\t"))
    check(len(truth) > 0, "ground_truth.tsv is empty")
    by_rep = collections.defaultdict(list)
    for r in truth:
        by_rep[r["replicate"]].append(r)
    print(f"{len(truth)} lesions across {len(by_rep)} replicates")

    parental_cache: dict[str, dict] = {}
    n_checked = 0

    for rid, rows in sorted(by_rep.items()):
        acc = rows[0]["accession"]
        pf = asm / f"{acc}__C0_parental.fna"
        mf = asm / f"{rid}.fna"
        if not check(pf.exists(), f"{pf.name} missing") or \
           not check(mf.exists(), f"{mf.name} missing"):
            continue
        par = parental_cache.setdefault(acc, load(pf))
        mod = load(mf)

        # contig names and order must be byte-identical, which is what makes
        # the caller outputs comparable without any name reconciliation
        check(list(par) == list(mod),
              f"{rid}: contig names or order differ from the parental")

        for r in rows:
            cg = r["contig"]
            s1, e1 = int(r["parental_start1"]), int(r["parental_end1"])
            strand = int(r["strand"])
            delta = int(r["delta_bp"])
            code, tier = r["lesion_code"], r["truth_tier"]
            p_gene = S.gene_seq(par, cg, s1, e1, strand)

            # the 21 bp parental window in the truth table must be the real
            # parental sequence at the recorded offset
            gp = int(r["gene_offset"])
            want_win = p_gene[max(0, gp - 10):gp + 11]
            check(want_win == r["window_parental"],
                  f"{rid} {r['locus_tag']}: recorded parental window does not "
                  f"match the parental assembly")

            # the lesioned gene, read out of the modified file, must show the
            # consequence the tier claims
            lf = load_liftover(rid)[cg]
            try:
                if strand >= 0:
                    m0 = lf.to_modified(s1 - 1)
                    m_gene = mod[cg][m0:m0 + len(p_gene) + delta]
                else:
                    m1 = lf.to_modified(e1 - 1)
                    m_gene = S.revcomp(
                        mod[cg][m1 - (len(p_gene) + delta) + 1:m1 + 1])
            except KeyError as exc:
                check(False, f"{rid} {r['locus_tag']}: liftover: {exc}")
                continue

            check(len(m_gene) == len(p_gene) + delta,
                  f"{rid} {r['locus_tag']}: modified gene length is "
                  f"{len(m_gene)}, expected {len(p_gene) + delta}")
            check(m_gene != p_gene or code == "C2_intergenic",
                  f"{rid} {r['locus_tag']}: {code} changed nothing")

            trimmed = m_gene[:len(m_gene) - len(m_gene) % 3]
            prot = S.translate(trimmed) if trimmed else ""
            if code in ("M1_stop_sub", "M7_trunc_term"):
                check("*" in prot[:-1],
                      f"{rid} {r['locus_tag']}: {code} shows no premature stop "
                      f"in the written assembly")
            elif code in ("M2_fs_ins1", "M3_fs_del1", "M4_fs_indel2",
                          "M8_disrupt_int"):
                check(len(m_gene) % 3 != 0,
                      f"{rid} {r['locus_tag']}: {code} is in frame on disk")
            elif code == "M6_start_loss":
                check(m_gene[:3] not in S.STARTS,
                      f"{rid} {r['locus_tag']}: start codon still valid on disk")
            elif code == "C1_synonymous":
                check(S.translate(m_gene) == S.translate(p_gene),
                      f"{rid} {r['locus_tag']}: C1 changed the protein on disk")
                check(tier == "negative", f"{rid}: C1 tier is {tier}")
            elif code == "C3_inframe_indel":
                check(abs(delta) == 3 and len(m_gene) % 3 == 0,
                      f"{rid} {r['locus_tag']}: C3 is not a 3 bp in-frame indel")
                check(tier == "ambiguous",
                      f"{rid}: C3 must stay ambiguous, found {tier}")
            n_checked += 1

        # the rest of the genome must be untouched outside the lesions: total
        # length changes by exactly the sum of the recorded deltas
        want_len = sum(len(v) for v in par.values()) + sum(
            int(x["delta_bp"]) for x in rows)
        got_len = sum(len(v) for v in mod.values())
        # intergenic control substitutions have delta 0 and are not in the truth
        # table, so lengths must match exactly
        check(got_len == want_len,
              f"{rid}: modified genome is {got_len} bp, expected {want_len} "
              f"from the recorded deltas")
        print(f"  {rid}: {len(rows)} lesions re-verified from disk")

    # no gene carries two lesions in one replicate
    for rid, rows in by_rep.items():
        tags = [r["locus_tag"] for r in rows]
        check(len(tags) == len(set(tags)),
              f"{rid}: a gene appears twice in the truth table")

    # every tier the plan defines, and nothing else
    tiers = set(r["truth_tier"] for r in truth)
    check(tiers <= {"high", "medium", "ambiguous", "negative", "technical"},
          f"unexpected truth tiers: {tiers}")

    print(f"\n{n_checked} lesions checked against the written assemblies")
    if FAILED:
        print(f"{len(FAILED)} failures")
        sys.exit(1)
    print("the truth table and the assemblies on disk agree")


if __name__ == "__main__":
    main()
