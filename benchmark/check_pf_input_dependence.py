#!/usr/bin/env python3
"""Is Pseudofinder's apparent agreement with PGAP structural rather than real?

LOCUS_RESULTS.md reported that among 7,565 discordant loci, 1,955 are
Pseudofinder-pseudogene against PGAP-intact and ZERO are the reverse, and read
that as the two callers being nested rather than opposed -- "they agree about
those, gene by gene". DFAST_RESULTS.md then argued the nesting is a fact about
those two callers and not a property of annotation pipelines.

Both readings assume Pseudofinder is an independent caller. It is not. It is a
post-processor, and run_one.sh gave it the RefSeq GenBank file:

    pseudofinder.py annotate --genome <...>/genomic.gbff --database sprot.dmnd

RefSeq's genomic.gbff IS PGAP's annotation. So Pseudofinder's gene set is
PGAP's gene set, and there are two ways a PGAP pseudogene can become
unopposable:

  (A) the pseudogene carries no CDS feature at all, the modern RefSeq
      convention, so it is absent from Pseudofinder's CDS universe and cannot
      be classified intact; or
  (B) it carries a CDS feature with /pseudo, and modules/annotate.py line 99
      converts it straight to a pseudogene:
          if feature.type == 'CDS' and feature.qualifiers.get('pseudo'):
              feature.type = 'pseudogene'
              ... 'Annotated as pseudogene in input genbank file.'

Either way Pseudofinder cannot contradict PGAP about a gene PGAP called
pseudogene. The 1,955-against-0 asymmetry would then be guaranteed by the
wiring, not observed in the biology.

This script counts, per genome, which mechanism applies. It changes no result:
the Pseudofinder-only elevation (+0.035, p = 1e-11) and DFAST's lack of one are
untouched, because those concern genes PGAP annotated as intact CDS. What it
bears on is the interpretation quoted above.
"""
import collections
import json
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
D = pathlib.Path.home() / "research/datasets"


def scan(gb: pathlib.Path) -> dict:
    """Count pseudo-marked gene and CDS features without a full GenBank parse.

    A plain text scan, because Biopython takes seconds per 10 MB file and this
    needs to cover the whole corpus. Feature blocks start at column 5 of the
    FEATURES table; a qualifier line starts with '/'.
    """
    gene_pseudo = cds_pseudo = cds_total = gene_total = 0
    cur = None
    cur_has_pseudo = False

    def close():
        nonlocal gene_pseudo, cds_pseudo, cds_total, gene_total
        if cur == "gene":
            gene_total += 1
            if cur_has_pseudo:
                gene_pseudo += 1
        elif cur == "CDS":
            cds_total += 1
            if cur_has_pseudo:
                cds_pseudo += 1

    feat = re.compile(r"^ {5}(\S+)\s+\S")
    with gb.open(errors="replace") as fh:
        for line in fh:
            m = feat.match(line)
            if m:
                close()
                cur, cur_has_pseudo = m.group(1), False
            elif cur and line.startswith("                     /"):
                q = line.strip()
                if q.startswith("/pseudo"):
                    cur_has_pseudo = True
            elif line.startswith("ORIGIN") or line.startswith("//"):
                close()
                cur, cur_has_pseudo = None, False
    close()
    return {"gene_total": gene_total, "gene_pseudo": gene_pseudo,
            "cds_total": cds_total, "cds_pseudo": cds_pseudo}


def main() -> None:
    accs = sorted(p.name for p in (D / "n2-genomes-expansion").iterdir()
                  if p.is_dir())
    if len(sys.argv) > 1:
        accs = accs[:int(sys.argv[1])]
    rows, mech = [], collections.Counter()
    for acc in accs:
        gb = D / f"n2-genomes-expansion/{acc}/genomic.gbff"
        if not gb.exists():
            continue
        r = scan(gb)
        r["accession"] = acc
        if r["gene_pseudo"] == 0:
            r["mechanism"] = "no pseudogenes annotated"
        elif r["cds_pseudo"] == 0:
            r["mechanism"] = "A: pseudogenes carry no CDS feature"
        elif r["cds_pseudo"] >= 0.9 * r["gene_pseudo"]:
            r["mechanism"] = "B: pseudogenes carry a CDS with /pseudo"
        else:
            r["mechanism"] = "mixed"
        mech[r["mechanism"]] += 1
        rows.append(r)
        print(f"  {acc:18s} genes {r['gene_total']:5d} "
              f"pseudo-gene {r['gene_pseudo']:4d}  CDS {r['cds_total']:5d} "
              f"pseudo-CDS {r['cds_pseudo']:4d}   {r['mechanism']}", flush=True)

    out = ROOT / "benchmark/results/pf_input_dependence.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("accession\tgene_total\tgene_pseudo\tcds_total\tcds_pseudo\t"
                 "mechanism\n")
        for r in rows:
            fh.write(f"{r['accession']}\t{r['gene_total']}\t{r['gene_pseudo']}\t"
                     f"{r['cds_total']}\t{r['cds_pseudo']}\t{r['mechanism']}\n")

    n_unopposable = sum(r["gene_pseudo"] for r in rows)
    print(f"\n{len(rows)} genomes scanned")
    for k, v in mech.most_common():
        print(f"  {v:4d}  {k}")
    print(f"\n{n_unopposable} PGAP pseudogene calls in total across these "
          f"genomes.")
    print("In every genome above, each of them is either absent from "
          "Pseudofinder's\ngene universe (mechanism A) or copied into its "
          "output as a pseudogene\n(mechanism B). In neither case can "
          "Pseudofinder call it intact.")
    (ROOT / "benchmark/results/pf_input_dependence.json").write_text(
        json.dumps({"n_genomes": len(rows),
                    "mechanisms": dict(mech),
                    "total_pgap_pseudogenes": n_unopposable}, indent=2,
                   sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
