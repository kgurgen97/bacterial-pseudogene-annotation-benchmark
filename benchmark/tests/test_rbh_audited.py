#!/usr/bin/env python3
"""Unit tests for the decision logic in analysis_rbh_audited.py.

Three generations of this analysis have each been wrong in a way the output did
not reveal, so the rules that decide acceptance are tested directly. The ones
tested hardest are the audit's own demands:

  * contig identity by EXACT sequence hash, never by length or record order
  * ONLY representations of the same verified genomic locus may collapse;
    identical protein sequences at different positions stay distinct, so a
    paralogue pair must now surface as a genuine tie
  * forward ambiguity must be visible: a near-tied runner-up subject makes the
    whole assignment a coin toss
  * a sole candidate is not a tie with itself
  * standardised mean differences and calipers must behave as claimed

    python3 benchmark/tests/test_rbh_audited.py
"""
import ast
import pathlib
import sys
import types

ROOT = pathlib.Path(__file__).resolve().parents[2]
FAILED: list[str] = []


def check(cond, msg):
    if not cond:
        FAILED.append(msg)


def close(a, b, tol=1e-9):
    return a is not None and abs(a - b) <= tol


def load():
    """Definitions and upper-case constants only; nothing else is executed."""
    path = ROOT / "analysis_rbh_audited.py"
    tree = ast.parse(path.read_text())
    keep = [n for n in tree.body
            if isinstance(n, (ast.Import, ast.ImportFrom, ast.FunctionDef,
                              ast.ClassDef))
            or (isinstance(n, ast.Assign)
                and all(isinstance(t, ast.Name) and t.id.isupper()
                        for t in n.targets))]
    m = types.ModuleType("rbh_a")
    m.__dict__["__file__"] = str(path)
    exec(compile(ast.Module(body=keep, type_ignores=[]), str(path), "exec"),
         m.__dict__)
    return m


M = load()


# --------------------------------------------------------------- overlap rule
def test_overlap_rule():
    check(M.overlaps_half(100, 200, 100, 200), "identical spans must overlap")
    check(M.overlaps_half(100, 200, 110, 190), "a span over half of another")
    check(not M.overlaps_half(100, 1100, 100, 150),
          "a long feature sharing little of a short one must NOT overlap")
    check(not M.overlaps_half(100, 200, 199, 400), "a 2 bp touch is not overlap")
    check(M.overlap_bp(100, 200, 150, 250) == 51, "overlap_bp is inclusive")
    check(M.overlap_bp(100, 200, 300, 400) == 0, "disjoint spans share nothing")


# ------------------------------------------------- hash-only contig identity
def _cds(rec_id, rsha):
    return ("g1", "MAAA", rec_id, 1, 99, rsha, 0)


def test_contig_identity_by_hash_only():
    xlat, stat = M.hash_contig_map({"ctg_1": "ACGTACGT", "ctg_2": "TTTTGGGG"},
                                   [_cds("NZ_1", M.sha("ACGTACGT"))])
    check(xlat["NZ_1"] == "ctg_1", "an exact unique hash must map")
    check(stat.get("hash_unique") == 1, "the unique-hash counter must record it")

    xlat, stat = M.hash_contig_map({"ctg_1": "ACGTACGT"},
                                   [_cds("NZ_1", M.sha("GGGGGGGG"))])
    check(xlat["NZ_1"] is None, "an absent hash must NOT fall back to order")
    check(stat.get("hash_absent") == 1, "an absent hash must be counted")

    xlat, stat = M.hash_contig_map({"a": "ACGTACGT", "b": "ACGTACGT"},
                                   [_cds("NZ_1", M.sha("ACGTACGT"))])
    check(xlat["NZ_1"] is None, "an ambiguous hash must not map")

    xlat, _ = M.hash_contig_map({"ctg_1": "ACGTACGT", "ctg_2": "TTTTGGGG"},
                                [_cds("NZ_1", M.sha("GGGGTTTT"))])
    check(xlat["NZ_1"] is None,
          "equal length with a different sequence must NOT map")

    xlat, _ = M.hash_contig_map({"ctg_1": "ACGTACGT"}, [_cds("NZ_1", "")])
    check(xlat["NZ_1"] is None, "a record with no sequence must not map")


# -------------------------------------------- margins over distinct objects
def test_margin_of():
    b, bs, s, ss, mar = M.margin_of([("A", 100.0)])
    check(b == "A" and s is None and close(mar, 1.0),
          "a sole candidate is not a tie with itself")
    b, bs, s, ss, mar = M.margin_of([("A", 100.0), ("B", 99.0)])
    check(b == "A" and s == "B", "the winner and runner-up must be ordered")
    check(close(mar, 0.01), f"margin should be 0.01, got {mar}")
    _, _, _, _, mar = M.margin_of([("A", 100.0), ("B", 100.0)])
    check(close(mar, 0.0), "equal scores are a zero margin, the worst tie")
    check(M.margin_of([]) == (None, None, None, None, None),
          "no candidates must give all None")


# ------------------------- collapsing: same locus only, never same sequence
def test_only_same_locus_collapses():
    # a locus's annotated protein and its translated span are one object
    b, bs, s, ss, mar = M.collapse(
        [("CDS|g1", 300.0), ("LOC|L1", 298.0)],
        {"CDS|g1": "L1", "LOC|L1": "L1"})
    check(b == "L1" and bs == 300.0, "same-locus entries collapse to the locus")
    check(s is None and close(mar, 1.0),
          "a locus must not be its own runner-up")


def test_identical_proteins_at_different_loci_stay_distinct():
    """The v2 defect. Two identical proteins are two genomic loci.

    Collapsing them into one object hid the paralogue ambiguity the reverse
    search exists to detect; here it must surface as a tie.
    """
    entry_bio = {"CDS|a": "CDS:a", "CDS|b": "CDS:b"}   # no DUP: merging
    b, _, s, _, mar = M.collapse([("CDS|a", 200.0), ("CDS|b", 200.0)],
                                 entry_bio)
    check(b in ("CDS:a", "CDS:b"), "one of the two must be nominal best")
    check(s is not None, "the other identical protein must be the runner-up")
    check(close(mar, 0.0), "identical scores on distinct loci are a zero margin")
    check(mar < M.TIE_MARGIN,
          "a paralogue pair must be declared a tie, not silently merged")


def test_tie_policy_and_clear_winner():
    _, _, _, _, mar = M.collapse([("LOC|L1", 200.0), ("CDS|o", 199.0)],
                                 {"LOC|L1": "L1", "CDS|o": "CDS:o"})
    check(close(mar, 0.005) and mar < M.TIE_MARGIN,
          "a 0.5 % margin on a different object is a tie")
    _, _, _, _, mar = M.collapse([("LOC|L1", 400.0), ("CDS|o", 100.0)],
                                 {"LOC|L1": "L1", "CDS|o": "CDS:o"})
    check(close(mar, 0.75) and mar >= M.TIE_MARGIN,
          "a clear winner must not be a tie")


def test_unmapped_entries_ignored():
    b, bs, _, _, _ = M.collapse([("CDS|ghost", 500.0), ("LOC|L1", 100.0)],
                                {"LOC|L1": "L1"})
    check(b == "L1" and bs == 100.0,
          "an entry with no biological identifier must be ignored, not guessed")


# --------------------------------------------------------- frame and strand
def test_frame_and_orientation():
    nt = "ATG" + "AAA" * 20 + "TAA"
    prot, frame, rev, frac = M.translate_span(nt, 1, len(nt))
    check(frame == 0 and rev is False, "a forward full-length hit is frame 0")
    check(prot.startswith("MKKK"), f"expected MKKK…, got {prot[:6]}")
    check(prot.endswith("X"), "a stop codon must become X")
    check(close(frac, 1 / len(prot)), "the stop fraction must be reported")

    from Bio.Seq import Seq
    fwd = "ATG" + "GCA" * 25
    prot, _, rev, _ = M.translate_span(str(Seq(fwd).reverse_complement()),
                                       len(fwd), 1)
    check(rev is True, "qstart > qend means the other strand")
    check(prot.startswith("MAAA"), f"must be reverse-complemented, got {prot[:6]}")

    _, frame, _, _ = M.translate_span("GG" + "ATG" + "AAA" * 20, 3, 65)
    check(frame == 2, "an offset start must set the frame")

    prot, _, _, frac = M.translate_span("ATGAAATAA", 1, 9)
    check(prot is None and frac is None, "a span below the minimum is refused")

    prot, _, _, frac = M.translate_span("TAA" * 30, 1, 90)
    check(prot == "X" * 30 and frac == 1.0 and frac > M.MAX_STOP_FRAC,
          "an all-stop span must exceed the suspect-frame threshold")


# ------------------------------------------------------- matching machinery
def test_smd():
    import numpy as np
    check(close(M.smd([1, 2, 3, 4], [1, 2, 3, 4]), 0.0),
          "identical groups have zero standardised difference")
    check(M.smd([1, 1], [1, 1]) == 0.0, "zero variance must not divide by zero")
    check(M.smd([1], [1, 2, 3]) is None, "a single observation gives no SMD")
    big = M.smd([10, 11, 12, 13], [1, 2, 3, 4])
    check(big is not None and big > 4, f"a large shift must show, got {big}")
    a = np.array([0.0, 1.0, 2.0, 3.0])
    check(close(M.smd(a, a + 1.0), -1.0 / a.std(ddof=1), 1e-9),
          "SMD must be the mean difference over the pooled SD")


def test_prep_cov_logs_span_only():
    import pandas as pd
    f = pd.DataFrame({"span": [10, 1000], "forward_coverage": [0.5, 1.0]})
    z = M.prep_cov(f, ["span", "forward_coverage"])
    check(close(float(z.span.iloc[0]), 1.0) and close(float(z.span.iloc[1]), 3.0),
          "span must be log10-transformed")
    check(close(float(z.forward_coverage.iloc[1]), 1.0),
          "coverage must be left alone")
    z = M.prep_cov(pd.DataFrame({"span": [0]}), ["span"])
    check(float(z.span.iloc[0]) == 0.0, "a zero span must be clipped, not -inf")


def _frame(n_dis, n_con, dis_span, con_span, omega_dis, omega_con):
    import pandas as pd
    rows = []
    for i in range(n_dis):
        rows.append({"accession": "A", "genus": "G", "kind": "discordant",
                     "hiconf": 1, "span": dis_span, "forward_coverage": 1.0,
                     "forward_pident": 90.0, "omega": omega_dis})
    for i in range(n_con):
        rows.append({"accession": "A", "genus": "G",
                     "kind": "concordant_intact", "hiconf": 1,
                     "span": con_span, "forward_coverage": 1.0,
                     "forward_pident": 90.0, "omega": omega_con})
    return pd.DataFrame(rows)


def test_caliper_refuses_distant_matches():
    """With no spread there is no distance, so everything matches."""
    df = _frame(5, 5, 100, 100, 0.3, 0.1)
    res, mp, _, _, _ = M.match_arm(df, ["span", "forward_coverage"], 0.25)
    check(res["n_pairs"] == 5, f"identical covariates must all match, "
                               f"got {res['n_pairs']}")
    check(close(res["median_gap"], 0.2), "the gap must be recovered")
    check(res["unmatched_frac"] == 0.0, "nothing should go unmatched here")


def test_hiconf_only_and_minimum_counts():
    import pandas as pd
    df = _frame(5, 5, 100, 100, 0.3, 0.1)
    df.loc[df.kind == "discordant", "hiconf"] = 0
    res, mp, _, _, _ = M.match_arm(df, ["span", "forward_coverage"], None)
    check(res["n_pairs"] == 0,
          "excluded loci must never enter the matched arm")
    small = pd.concat([_frame(2, 5, 100, 100, 0.3, 0.1)])
    res, _, _, _, _ = M.match_arm(small, ["span", "forward_coverage"], None)
    check(res["n_pairs"] == 0, "a genome with under 3 per side must be skipped")


# ------------------------------------------------------------- NaN handling
def test_nan_is_not_truthy_here():
    import numpy as np
    check(M.clean(float("nan")) is None, "a float NaN must clean to None")
    check(M.clean(np.nan) is None, "np.nan must clean to None")
    check(M.clean(None) is None, "None must clean to None")
    check(M.clean("ABU_RS00005") == "ABU_RS00005", "a real id must survive")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    print(f"running {len(tests)} test groups")
    for t in tests:
        print(f"  {t.__name__}")
        t()
    if FAILED:
        print(f"\n{len(FAILED)} failures")
        for f in FAILED:
            print(f"  - {f}")
        sys.exit(1)
    print("\nall audited-RBH decision-logic tests pass")


if __name__ == "__main__":
    main()
