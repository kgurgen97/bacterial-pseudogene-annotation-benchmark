#!/usr/bin/env python3
"""Re-derive the eligible gene set, which the simulator counted but never listed.

simulate.py wrote results/eligibility.tsv -- exclusion counts by reason -- but
not the genes themselves. Without them there is no negative universe, so no
false-positive rate and no specificity: plan section 7 defines a false positive
as a gene called positive in a modified replicate that carried NO lesion and
was not already positive in the parental.

Re-deriving is safe because eligibility is fully deterministic: the eight
criteria in plan section 3 involve no RNG. The lesion DRAW uses a seed; the
eligible set does not. This imports simulate.py and calls its own functions
rather than reimplementing them, so the set cannot drift from the one the
lesions were drawn from. The row count is checked against eligibility.tsv's
recorded ELIGIBLE figure, which is the independent cross-check.

    python3 benchmark/dump_eligible.py
"""
import csv
import pathlib
import sys
import tempfile

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import simulate as S  # noqa: E402


def recorded_counts():
    """accession -> the ELIGIBLE count simulate.py recorded at the time."""
    out = {}
    p = HERE / "results/eligibility.tsv"
    if not p.exists():
        return out
    with p.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["reason"] == "ELIGIBLE":
                out[r["accession"]] = int(r["n"])
    return out


def main() -> None:
    S.DIAMOND = S.find_diamond()
    want = recorded_counts()
    truth = list(csv.DictReader((HERE / "results/ground_truth.tsv").open(),
                                delimiter="\t"))
    accs = sorted({r["accession"] for r in truth})
    print(f"re-deriving the eligible set for {len(accs)} genomes")
    print(f"diamond: {S.DIAMOND}")

    rows, bad = [], []
    for acc in accs:
        contigs, cds, feats, fna, notes = S.load_parental(acc)
        if notes.get("reject"):
            bad.append(f"{acc}: {notes['reject']}")
            continue
        with tempfile.TemporaryDirectory() as td:
            elig, why = S.eligible(acc, contigs, cds, feats,
                                   set(notes.get("plasmids", [])),
                                   pathlib.Path(td))
        n_want = want.get(acc)
        ok = n_want is None or len(elig) == n_want
        print(f"  {acc:20s} eligible {len(elig):5d}"
              f"  recorded {n_want}  {'ok' if ok else 'MISMATCH'}")
        if not ok:
            bad.append(f"{acc}: re-derived {len(elig)} but simulate.py "
                       f"recorded {n_want}")
        for c in elig:
            rows.append({"accession": acc, "locus_tag": c.locus_tag,
                         "protein_id": c.protein_id, "contig": c.contig,
                         "start1": c.start1, "end1": c.end1,
                         "strand": c.strand, "n_codons": c.n_codons})

    out = HERE / "results/eligible_genes.tsv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {out.relative_to(HERE.parent)} with {len(rows)} genes")

    # every lesioned gene must be in the eligible set, or the truth table
    # describes genes the simulator could not have drawn
    keyed = {(r["accession"], r["locus_tag"]) for r in rows}
    missing = [(r["accession"], r["locus_tag"]) for r in truth
               if (r["accession"], r["locus_tag"]) not in keyed]
    print(f"lesioned genes absent from the re-derived eligible set: "
          f"{len(missing)}")
    for m in missing[:5]:
        print(f"  {m}")
    if missing:
        bad.append(f"{len(missing)} lesioned genes are not in the eligible set")

    if bad:
        print("\nPROBLEMS")
        for b in bad:
            print(f"  {b}")
        sys.exit(1)
    print("the re-derived set matches simulate.py's recorded counts exactly")


if __name__ == "__main__":
    main()
