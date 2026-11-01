#!/usr/bin/env python3
"""Unit tests for the decision logic in analysis_rbh_validated.py.

These test the parts that decide whether a locus is accepted, because that is
where version 1 of the RBH analysis was wrong in a way no output inspection
would have caught: a contig-namespace mismatch silently disabled two of the
three acceptance paths and the whole synteny test, and the script still printed
a plausible-looking number (32.7 % reciprocation, zero synteny votes).

The rules tested hardest are the ones an adversarial reviewer asked for:

  * contig identity by EXACT sequence hash, with no positional fall-back, and
    equal length refused as evidence
  * two FASTA entries for one biological object must not compete or tie
  * a close score on a DIFFERENT object must be declared a tie, not resolved
  * the translated span must inherit the forward hit's frame and orientation
  * NaN must never pass a truthiness test

    python3 benchmark/tests/test_rbh_validated.py
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
    """Load the module's own definitions without executing anything else.

    The script's top level borrows sibling modules and prints plan hashes; this
    keeps definitions and upper-case constants only, the idiom the project uses
    elsewhere for the same reason.
    """
    path = ROOT / "analysis_rbh_validated.py"
    tree = ast.parse(path.read_text())
    keep = [n for n in tree.body
            if isinstance(n, (ast.Import, ast.ImportFrom, ast.FunctionDef,
                              ast.ClassDef))
            or (isinstance(n, ast.Assign)
                and all(isinstance(t, ast.Name) and t.id.isupper()
                        for t in n.targets))]
    m = types.ModuleType("rbh_v")
    m.__dict__["__file__"] = str(path)
    exec(compile(ast.Module(body=keep, type_ignores=[]), str(path), "exec"),
         m.__dict__)
    return m


M = load()


# --------------------------------------------------------------- overlap rule
def test_overlap_rule():
    check(M.overlaps_half(100, 200, 100, 200), "identical spans must overlap")
    check(M.overlaps_half(100, 200, 110, 190),
          "a span inside another, over half of it, must overlap")
    check(not M.overlaps_half(100, 1100, 100, 150),
          "a long feature sharing a little of a short one must NOT overlap")
    check(not M.overlaps_half(100, 200, 199, 400), "a 2 bp touch is not overlap")
    check(not M.overlaps_half(100, 200, 300, 400), "disjoint is not overlap")


# ------------------------------------------------- hash-only contig identity
def _cds(rec_id, rsha):
    return ("g1", "MAAA", rec_id, 1, 99, rsha, 0)


def test_contig_identity_by_hash():
    contigs = {"ctg_1": "ACGTACGT", "ctg_2": "TTTTGGGG"}
    xlat, stat = M.hash_contig_map(contigs, [_cds("NZ_1", M.sha("ACGTACGT"))])
    check(xlat["NZ_1"] == "ctg_1", "an exact unique hash must map")
    check(stat.get("hash_unique") == 1, "the unique-hash counter must record it")

    # the version 1 defect: no positional guess is ever offered
    xlat, stat = M.hash_contig_map({"ctg_1": "ACGTACGT"},
                                   [_cds("NZ_1", M.sha("GGGGGGGG"))])
    check(xlat["NZ_1"] is None, "an absent hash must NOT fall back to order")
    check(stat.get("hash_absent") == 1, "an absent hash must be counted")

    # two identical contigs cannot identify either one
    xlat, stat = M.hash_contig_map({"a": "ACGTACGT", "b": "ACGTACGT"},
                                   [_cds("NZ_1", M.sha("ACGTACGT"))])
    check(xlat["NZ_1"] is None, "an ambiguous hash must not map")
    check(stat.get("hash_ambiguous_in_assembly") == 1,
          "an ambiguous hash must be counted")

    # length matching, which version 1 used, would have accepted this
    xlat, _ = M.hash_contig_map({"ctg_1": "ACGTACGT", "ctg_2": "TTTTGGGG"},
                                [_cds("NZ_1", M.sha("GGGGTTTT"))])
    check(xlat["NZ_1"] is None,
          "equal length with a different sequence must NOT map")

    # an empty hash (no sequence in the record) is not a match either
    xlat, _ = M.hash_contig_map({"ctg_1": "ACGTACGT"}, [_cds("NZ_1", "")])
    check(xlat["NZ_1"] is None, "a record with no sequence must not map")


# ------------------------------------------- collapsing entries to organisms
def test_entries_collapse_to_biological_objects():
    # the annotated protein and the translated span are one object; without
    # collapsing, the runner-up would be the same gene and read as a tie
    best, bs, second, ss, margin = M.collapse(
        [("CDS|g1", 300.0), ("LOC|L1", 298.0)],
        {"CDS|g1": "L1", "LOC|L1": "L1"})
    check(best == "L1" and bs == 300.0, "the locus must win with its best score")
    check(second is None and ss is None,
          "its own second entry must not be the runner-up")
    check(close(margin, 1.0), "a sole object must have a full margin")

    best, _, second, _, margin = M.collapse(
        [("CDS|a", 200.0), ("CDS|b", 200.0)],
        {"CDS|a": "DUP:x", "CDS|b": "DUP:x"})
    check(best == "DUP:x" and second is None,
          "identical duplicated proteins are one object")
    check(close(margin, 1.0), "a duplicate family must not tie with itself")


def test_tie_policy():
    best, _, second, _, margin = M.collapse(
        [("LOC|L1", 200.0), ("CDS|other", 199.0)],
        {"LOC|L1": "L1", "CDS|other": "CDS:other"})
    check(best == "L1" and second == "CDS:other",
          "a different object must be the runner-up")
    check(close(margin, 0.005), f"margin should be 0.005, got {margin}")
    check(margin < M.TIE_MARGIN, "a 0.5 % margin must count as a tie")

    _, _, _, _, margin = M.collapse(
        [("LOC|L1", 400.0), ("CDS|other", 100.0)],
        {"LOC|L1": "L1", "CDS|other": "CDS:other"})
    check(close(margin, 0.75), f"margin should be 0.75, got {margin}")
    check(margin >= M.TIE_MARGIN, "a clear winner must not be a tie")


def test_unmapped_and_empty_hits():
    best, bs, _, _, _ = M.collapse([("CDS|ghost", 500.0), ("LOC|L1", 100.0)],
                                   {"LOC|L1": "L1"})
    check(best == "L1" and bs == 100.0,
          "an entry with no biological identifier must be ignored, not guessed")
    check(M.collapse([], {}) == (None, None, None, None, None),
          "no hits must give all None")


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
    check(rev is True, "qstart > qend means the hit is on the other strand")
    check(prot.startswith("MAAA"),
          f"a reverse hit must be reverse-complemented, got {prot[:6]}")

    _, frame, _, _ = M.translate_span("GG" + "ATG" + "AAA" * 20, 3, 65)
    check(frame == 2, "an offset start must set the frame")

    prot, _, _, frac = M.translate_span("ATGAAATAA", 1, 9)
    check(prot is None and frac is None,
          "a span below the minimum must be refused, not translated")

    prot, _, _, frac = M.translate_span("TAA" * 30, 1, 90)
    check(prot == "X" * 30, "an all-stop span must translate to all X")
    check(frac == 1.0 and frac > M.MAX_STOP_FRAC,
          "an all-stop span must exceed the suspect-frame threshold")


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
    print("\nall RBH decision-logic tests pass")


if __name__ == "__main__":
    main()
