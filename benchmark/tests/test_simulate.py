#!/usr/bin/env python3
"""Unit tests for the simulator, with the strand arithmetic done by hand.

Runs with plain python3 (no pytest needed):  python3 tests/test_simulate.py

The test that earns its keep is the minus-strand one. A lesion placed in a
minus-strand gene with the sign wrong still lands on real sequence and still
breaks a gene -- just not the gene the truth table names. Nothing downstream
would notice, and every accuracy number would be quietly wrong for roughly half
the genes. So the expected genomic coordinates below are written out longhand
rather than derived from the function under test.
"""
import pathlib
import random
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import simulate as S  # noqa: E402

FAILED: list[str] = []


def check(cond, msg):
    if not cond:
        FAILED.append(msg)
        print(f"  FAIL  {msg}")
    return cond


# --------------------------------------------------------------------------
def test_apply_edits():
    s = "AAAACCCCGGGGTTTT"
    check(S.apply_edits(s, [S.Edit("c", 0, 1, "G")]) == "GAAACCCCGGGGTTTT",
          "substitution at position 0")
    check(S.apply_edits(s, [S.Edit("c", 4, 0, "TT")]) == "AAAATTCCCCGGGGTTTT",
          "insertion of 2 bp at position 4")
    check(S.apply_edits(s, [S.Edit("c", 4, 4, "")]) == "AAAAGGGGTTTT",
          "deletion of 4 bp at position 4")
    check(S.apply_edits(s, [S.Edit("c", 1, 1, "T"), S.Edit("c", 12, 0, "AA")])
          == "ATAACCCCGGGGAATTTT", "two edits applied in coordinate order")
    check(S.apply_edits(s, []) == s, "no edits leaves the sequence alone")


def test_liftover():
    # delete 4 bp at 4, insert 2 bp at 12
    edits = [S.Edit("c", 4, 4, ""), S.Edit("c", 12, 0, "TT")]
    lf = S.Liftover(edits)
    check(lf.to_modified(0) == 0, "position before any edit is unmoved")
    check(lf.to_modified(3) == 3, "position just before a deletion is unmoved")
    check(lf.to_modified(8) == 4, "position after a 4 bp deletion shifts by -4")
    # the inserted bases precede parental[12], so parental[12] carries the
    # insertion's +2 as well as the deletion's -4: 12 - 4 + 2 = 10.
    check(lf.to_modified(12) == 10,
          "position at an insertion point is pushed right by the insertion")
    check(lf.to_modified(13) == 11, "position after both edits shifts by -2")
    try:
        lf.to_modified(5)
        check(False, "a position inside a deletion must raise, not map silently")
    except KeyError:
        pass

    # an equal-length replacement is NOT a deletion: every base in it has an
    # image. Three bases replaced by three at position 20 is what a start-codon
    # lesion does, and refusing it silently broke every M6 lesion.
    lf2 = S.Liftover([S.Edit("c", 20, 3, "GGG")])
    for k in range(3):
        check(lf2.to_modified(20 + k) == 20 + k,
              f"replaced base {20 + k} must map to itself, not raise")
    # a shrinking replacement keeps images only for the bases it still covers
    lf3 = S.Liftover([S.Edit("c", 20, 3, "G")])
    check(lf3.to_modified(20) == 20, "first base of a shrinking replacement maps")
    for k in (1, 2):
        try:
            lf3.to_modified(20 + k)
            check(False, f"base {20 + k} was deleted and must raise")
        except KeyError:
            pass


def test_liftover_matches_apply_edits():
    """The liftover must agree with the sequence the edits actually produce."""
    rng = random.Random(7)
    for trial in range(200):
        n = rng.randrange(60, 200)
        seq = "".join(rng.choice("ACGT") for _ in range(n))
        edits, pos = [], 5
        while pos < n - 10:
            kind = rng.choice(("sub", "ins", "del", "skip"))
            if kind == "sub":
                edits.append(S.Edit("c", pos, 1, rng.choice("ACGT")))
                pos += 1
            elif kind == "ins":
                edits.append(S.Edit("c", pos, 0, rng.choice("ACGT")))
            elif kind == "del":
                k = rng.randrange(1, 4)
                edits.append(S.Edit("c", pos, k, ""))
                pos += k
            pos += rng.randrange(3, 12)
        mod = S.apply_edits(seq, edits)
        lf = S.Liftover(edits)
        # a substituted position maps to the NEW base, so check against the
        # replacement rather than skipping it -- skipping it is what let an
        # equal-length-replacement bug through.
        subbed = {e.pos: e.alt for e in edits if e.ref_len == 1 and e.alt}
        for p in range(n):
            try:
                m = lf.to_modified(p)
            except KeyError:
                continue
            want = subbed[p] if p in subbed else seq[p]
            if not check(mod[m] == want,
                         f"trial {trial}: liftover sent parental {p} to "
                         f"{mod[m]!r}, expected {want!r}"):
                return


# --------------------------------------------------------------------------
def build_contig():
    """A contig holding one plus-strand and one minus-strand gene, by hand."""
    gene = "ATG" + "AAACCCGGGTTTAAACCCGGG" + "TAA"   # 9 codons, no internal stop
    assert len(gene) % 3 == 0
    assert "*" not in S.translate(gene)[:-1]
    left = "CGCGCGCGCG"                               # 10 bp
    mid = "TTTTTTTTTT"                                 # 10 bp
    right = "GCGCGCGCGC"                               # 10 bp
    seq = left + gene + mid + S.revcomp(gene) + right
    plus = dict(contig="c1", start1=len(left) + 1, end1=len(left) + len(gene),
                strand=1)
    m0 = len(left) + len(gene) + len(mid)
    minus = dict(contig="c1", start1=m0 + 1, end1=m0 + len(gene), strand=-1)
    return {"c1": seq}, gene, plus, minus


def test_gene_extraction_both_strands():
    contigs, gene, plus, minus = build_contig()
    check(S.gene_seq(contigs, **plus) == gene, "plus-strand gene extracted")
    check(S.gene_seq(contigs, **minus) == gene, "minus-strand gene extracted")


def test_plus_strand_coordinates():
    contigs, gene, plus, _ = build_contig()
    e = S.gene_edit_to_genome(0, 1, "C", plus["start1"], plus["end1"],
                              1, "c1", "t")
    check(e.pos == plus["start1"] - 1,
          f"plus strand gene offset 0 is genomic {plus['start1'] - 1}, got {e.pos}")
    check(e.alt == "C", "plus strand replacement is not complemented")
    e2 = S.gene_edit_to_genome(7, 1, "T", plus["start1"], plus["end1"],
                               1, "c1", "t")
    check(e2.pos == plus["start1"] - 1 + 7, "plus strand offset 7")


def test_minus_strand_coordinates():
    contigs, gene, _, minus = build_contig()
    s1, e1 = minus["start1"], minus["end1"]

    # gene offset 0 on the minus strand is the genomic base at e1 (1-based),
    # i.e. 0-based e1 - 1, and the base must be complemented.
    ed = S.gene_edit_to_genome(0, 1, "C", s1, e1, -1, "c1", "t")
    check(ed.pos == e1 - 1,
          f"minus strand offset 0 is genomic {e1 - 1}, got {ed.pos}")
    check(ed.alt == "G", f"minus strand 'C' must be stored as 'G', got {ed.alt}")

    # offset 6, deleting 2 bases: gene offsets 6 and 7 are genomic e1-7 and
    # e1-8 (0-based), so the replaced stretch starts at e1-8.
    ed = S.gene_edit_to_genome(6, 2, "", s1, e1, -1, "c1", "t")
    check(ed.pos == e1 - 8,
          f"minus strand 2 bp deletion at offset 6 starts at {e1 - 8}, got {ed.pos}")

    # an insertion before offset 6 goes between genomic e1-7 and e1-6, which is
    # insertion index e1-6.
    ed = S.gene_edit_to_genome(6, 0, "A", s1, e1, -1, "c1", "t")
    check(ed.pos == e1 - 6,
          f"minus strand insertion at offset 6 is index {e1 - 6}, got {ed.pos}")
    check(ed.alt == "T", "minus strand inserted 'A' must be stored as 'T'")


def test_minus_strand_round_trip():
    """Apply a minus-strand lesion and re-extract the gene through the liftover.

    This is the exact path validate() takes, checked here against a gene-space
    expectation computed without reference to any genomic coordinate.
    """
    for gene_pos, ref_len, alt in [(0, 1, "C"), (6, 0, "A"), (9, 1, ""),
                                   (12, 2, ""), (15, 0, "GG"), (3, 3, "")]:
        contigs, gene, _, minus = build_contig()
        s1, e1, cg = minus["start1"], minus["end1"], "c1"
        want = S.apply_edits(gene, [S.Edit("g", gene_pos, ref_len, alt)])
        ed = S.gene_edit_to_genome(gene_pos, ref_len, alt, s1, e1, -1, cg, "t")
        mod = {cg: S.apply_edits(contigs[cg], [ed])}
        lf = S.Liftover([ed])
        delta = len(alt) - ref_len
        # mirror validate(): anchor on the gene's genomic right edge
        try:
            anchor = lf.to_modified(e1 - 1)
        except KeyError:
            anchor = e1 - 1 + lf.delta_before(e1 - 1)
        got = S.revcomp(mod[cg][anchor - (len(gene) + delta) + 1:anchor + 1])
        check(got == want,
              f"minus strand round trip failed for ({gene_pos},{ref_len},{alt!r}): "
              f"expected {want[:18]}..., got {got[:18]}...")


# --------------------------------------------------------------------------
def test_make_stop():
    rng = random.Random(3)
    _, gene, _, _ = build_contig()
    found = 0
    for codon in range(1, len(gene) // 3 - 1):
        hit = S.make_stop(gene, codon, rng)
        if hit is None:
            continue
        found += 1
        gp, b = hit
        new = S.apply_edits(gene, [S.Edit("g", gp, 1, b)])
        check(new[gp // 3 * 3:gp // 3 * 3 + 3] in S.STOPS,
              f"make_stop did not create a stop at codon {codon}")
        check("*" in S.translate(new)[:-1],
              f"make_stop left no premature stop at codon {codon}")
    check(found > 0, "make_stop never found a one-substitution stop to make")


def test_synonymous_is_synonymous():
    contigs, gene, plus, _ = build_contig()
    c = S.Cds("t1", "p1", "c1", plus["start1"], plus["end1"], 1, len(gene) // 3)
    les = S.draw_lesion(c, contigs, "C1_synonymous", random.Random(11))
    if check(les is not None, "C1_synonymous produced no lesion"):
        new = S.apply_edits(gene, [S.Edit("g", gp, rl, a)
                                   for gp, rl, a in les.gene_edits])
        check(S.translate(new) == S.translate(gene),
              "C1_synonymous changed the protein")
        check(new != gene, "C1_synonymous changed nothing at all")
        check(les.tier == "negative", "C1_synonymous must be tier negative")


def test_frameshift_arms_change_frame():
    contigs, gene, plus, _ = build_contig()
    c = S.Cds("t1", "p1", "c1", plus["start1"], plus["end1"], 1, len(gene) // 3)
    for code in ("M2_fs_ins1", "M3_fs_del1", "M4_fs_indel2"):
        ok = 0
        for s in range(40):
            les = S.draw_lesion(c, contigs, code, random.Random(s))
            if les is None:
                continue
            new = S.apply_edits(gene, [S.Edit("g", gp, rl, a)
                                       for gp, rl, a in les.gene_edits])
            if len(new) % 3 == 0:
                check(False, f"{code} left the frame intact")
                break
            ok += 1
        check(ok > 0, f"{code} never produced a lesion on the toy gene")


def test_tiers_never_claim_too_much():
    """A late lesion is medium, never high; an ambiguous start loss is ambiguous."""
    contigs, gene, plus, _ = build_contig()
    c = S.Cds("t1", "p1", "c1", plus["start1"], plus["end1"], 1, len(gene) // 3)
    for s in range(60):
        for code in S.POSITIONED:
            les = S.draw_lesion(c, contigs, code, random.Random(s))
            if les is None:
                continue
            if les.position == "late":
                check(les.tier == "medium",
                      f"{code} late must be tier medium, got {les.tier}")
            else:
                check(les.tier == "high",
                      f"{code} {les.position} must be tier high, got {les.tier}")
    les = S.draw_lesion(c, contigs, "M8_disrupt_int", random.Random(5))
    if les:
        check(les.tier == "medium", "M8 must be tier medium")
    les = S.draw_lesion(c, contigs, "C3_inframe_indel", random.Random(5))
    if les:
        check(les.tier == "ambiguous", "C3 must be tier ambiguous, never negative")


def test_plan_hash_is_enforced():
    check(S.PLAN.exists(), "the frozen plan file is missing")
    import hashlib
    h = hashlib.sha256(S.PLAN.read_bytes()).hexdigest()
    check(h == S.PLAN_SHA,
          f"plan hash drift: simulate.py expects {S.PLAN_SHA[:16]}, file is {h[:16]}")


# --------------------------------------------------------------------------
def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    print(f"running {len(tests)} test groups")
    for t in tests:
        print(f"  {t.__name__}")
        t()
    if FAILED:
        print(f"\n{len(FAILED)} failures")
        sys.exit(1)
    print("\nall simulator tests pass")


if __name__ == "__main__":
    main()
