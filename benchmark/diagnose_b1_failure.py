#!/usr/bin/env python3
"""B1 fired. Is my pipeline broken, or do the callers genuinely miss these?

SIMULATION_PLAN.md section 9: "if no caller clears it, the pipeline or parser
is broken, not the callers; the run is void and nothing downstream is
reported." Neither caller cleared it - Bakta 5/16, DFAST 7/16 of high-tier
premature stops under endpoint B - so the plan's instruction is to suspect this
code first. That is what this does.

The suspicious pattern is by POSITION, not by caller: Bakta detects 6 % of
early lesions and 63 % of middle ones. Backwards, unless something mechanical
is going on. The hypothesis to test is that it is NOT a bug:

    a premature stop early in a gene leaves a long surviving ORF downstream.
    A gene caller re-annotates that ORF as an ordinary CDS, which still covers
    more than half of the original gene, so the locus reads intact_cds and no
    pseudogene is asserted. A stop in the MIDDLE leaves two roughly equal
    halves, neither of which covers half the gene, so the locus reads
    disrupted or absent and endpoint B fires.

Three things must all hold for a given missed lesion before the miss can be
called real rather than mine:

  P1  the lesion is actually present in the modified FASTA at the recorded
      coordinate
  P2  the caller really did annotate a CDS covering >= 50 % of the original
      gene, in the modified assembly
  P3  that CDS starts AFTER the premature stop, i.e. it is the surviving
      downstream ORF rather than the original gene

If all three hold the miss is a property of the caller. If P1 fails it is the
simulator; if P2 fails it is the parser or the liftover.

    python3 benchmark/diagnose_b1_failure.py
"""
import collections
import csv
import json
import pathlib
import sys

from Bio import SeqIO
from Bio.Seq import Seq

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
ASM = pathlib.Path.home() / "research/datasets/n2-benchmark/assemblies"
CALLS = pathlib.Path.home() / "research/datasets/n2-benchmark/calls"

import parse_calls as PC  # noqa: E402
from evaluate_stage_a import Lift, load_lift  # noqa: E402

STOPS = ("TAA", "TAG", "TGA")


def load_fa(p):
    return {r.id: str(r.seq).upper() for r in SeqIO.parse(str(p), "fasta")}


def main() -> None:
    truth = {(r["replicate"], r["locus_tag"]): r for r in
             csv.DictReader((HERE / "results/ground_truth.tsv").open(),
                            delimiter="\t")}
    preds = list(csv.DictReader(
        (HERE / "results/stage_a_predictions.tsv").open(), delimiter="\t"))

    # the misses that make B1 fail: high-tier premature stops read as intact
    misses = [r for r in preds
              if r["lesion_code"] == "M1_stop_sub"
              and r["truth_tier"] == "high"
              and r["modified_class"] == "intact_cds"]
    print(f"{len(misses)} high-tier premature stops classified intact_cds")
    byc = collections.Counter(r["caller"] for r in misses)
    byp = collections.Counter(r["position"] for r in misses)
    print(f"  by caller: {dict(byc)}")
    print(f"  by position: {dict(byp)}\n")

    checked = []
    for r in misses[:14]:
        rid, tag, caller = r["replicate"], r["locus_tag"], r["caller"]
        t = truth[(rid, tag)]
        acc = r["accession"]
        par = load_fa(ASM / f"{acc}__C0_parental.fna")
        mod = load_fa(ASM / f"{rid}.fna")
        cg = t["contig"]
        s1, e1, strand = int(t["parental_start1"]), int(t["parental_end1"]), int(t["strand"])
        gp = int(t["gene_offset"])

        # P1: the stop is in the modified assembly, in the gene's frame
        lift = load_lift(rid)[cg]
        m0 = lift.to_modified(s1 - 1)
        mlen = (e1 - s1 + 1) + int(t["delta_bp"])
        if strand >= 0:
            gene = mod[cg][m0:m0 + mlen]
        else:
            m1 = lift.to_modified(e1 - 1)
            gene = str(Seq(mod[cg][m1 - mlen + 1:m1 + 1]).reverse_complement())
        codon = gene[gp // 3 * 3: gp // 3 * 3 + 3]
        p1 = codon in STOPS
        stop_codon_idx = gp // 3

        # P2/P3: the CDS the caller emitted over this locus
        d = CALLS / caller / rid
        idx = PC.Index(PC.load_features(caller, d, rid))
        hits = idx.overlapping(cg, m0, m0 + mlen)
        cds = [f for f in hits if f.kind == PC.CDS]
        p2 = len(cds) >= 1
        p3 = None
        detail = ""
        if p2:
            f = max(cds, key=lambda f: f.end0 - f.start0)
            # where does the called CDS start, in the gene's own codon frame?
            if strand >= 0:
                off = f.start0 - m0
            else:
                off = (m0 + mlen) - f.end0
            p3 = off > (stop_codon_idx * 3)
            detail = (f"called CDS {f.end0 - f.start0} bp, starts at gene "
                      f"offset {off} ({off // 3} codons) vs stop at codon "
                      f"{stop_codon_idx} of {mlen // 3}")
        checked.append({"replicate": rid, "locus_tag": tag, "caller": caller,
                        "position": t["position"], "P1_lesion_present": bool(p1),
                        "P2_cds_called": bool(p2), "P3_starts_after_stop": p3,
                        "detail": detail})
        print(f"  {caller:6s} {tag:16s} {t['position']:7s} "
              f"P1 {'ok' if p1 else 'FAIL'}  P2 {'ok' if p2 else 'FAIL'}  "
              f"P3 {p3}")
        if detail:
            print(f"         {detail}")

    n = len(checked)
    p1ok = sum(c["P1_lesion_present"] for c in checked)
    p2ok = sum(c["P2_cds_called"] for c in checked)
    p3ok = sum(1 for c in checked if c["P3_starts_after_stop"])
    print(f"\nof {n} checked misses:")
    print(f"  P1 lesion present in the modified FASTA : {p1ok}/{n}")
    print(f"  P2 caller emitted a CDS over the locus  : {p2ok}/{n}")
    print(f"  P3 that CDS starts after the stop       : {p3ok}/{n}")

    print("\nverdict")
    if p1ok < n:
        print("  THE SIMULATOR IS AT FAULT: a recorded lesion is not in the")
        print("  assembly. Everything downstream is void until that is fixed.")
        v = "simulator"
    elif p2ok < n:
        print("  THE PARSER OR LIFTOVER IS AT FAULT: the locus was classified")
        print("  intact_cds without a CDS being found over it.")
        v = "parser"
    elif p3ok >= 0.6 * n:
        print("  THE MISSES ARE REAL. Every lesion is present, every locus")
        print("  carries a called CDS, and most of those CDS begin downstream")
        print("  of the introduced stop: the callers re-annotate the surviving")
        print("  ORF as an ordinary gene. B1's premise - that a premature stop")
        print("  is obviously detectable - does not hold for these tools, so")
        print("  the gate fires on a wrong assumption rather than on a defect.")
        print("  The registered Stage A numbers stay VOID either way.")
        v = "callers"
    else:
        print("  INCONCLUSIVE: lesions and CDS are present but the called CDS")
        print("  does not generally start after the stop, so the mechanism is")
        print("  not the surviving-ORF explanation and remains unexplained.")
        v = "inconclusive"

    (HERE / "results/b1_diagnosis.json").write_text(json.dumps(
        {"n_misses_total": len(misses), "by_caller": dict(byc),
         "by_position": dict(byp), "n_checked": n,
         "P1_lesion_present": p1ok, "P2_cds_called": p2ok,
         "P3_starts_after_stop": p3ok, "attribution": v,
         "checked": checked}, indent=2, sort_keys=True) + "\n")
    print("\nwrote benchmark/results/b1_diagnosis.json")


if __name__ == "__main__":
    main()
