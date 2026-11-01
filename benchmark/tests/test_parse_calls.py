#!/usr/bin/env python3
"""Fixture tests for the output parser: every class, from hand-written output.

SIMULATION_PLAN.md section 6 requires the parser to be driven by a fixture
covering every class, and its tests to be part of the Stage A gate. The rules
most likely to be got wrong, and therefore tested hardest:

  * DFAST `partial` must come out `partial`, not `explicit_pseudogene`
  * a missing CDS must come out `absent`, not a pseudogene call
  * a pseudogene assertion plus an intact CDS over one locus is `ambiguous`,
    not silently resolved
  * Bakta's two keys, `pseudogene=` and `pseudo=`, must BOTH assert
  * a gene split in two must come out `disrupted`

    python3 benchmark/tests/test_parse_calls.py
"""
import pathlib
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import parse_calls as PC  # noqa: E402

FAILED: list[str] = []


def check(cond, msg):
    if not cond:
        FAILED.append(msg)
        print(f"  FAIL  {msg}")


BAKTA_GFF = """\
##gff-version 3
##sequence-region c1 1 20000
c1\tbakta\tCDS\t1001\t1600\t.\t+\t0\tID=g_intact;product=ordinary protein
c1\tbakta\tCDS\t2001\t2600\t.\t+\t0\tID=g_unitary;product=x;pseudogene=unitary
c1\tbakta\tCDS\t3001\t3600\t.\t+\t0\tID=g_true;product=x;pseudo=True
c1\tbakta\tCDS\t4001\t4280\t.\t+\t0\tID=g_split_a;product=x
c1\tbakta\tCDS\t4320\t4600\t.\t+\t0\tID=g_split_b;product=x
c1\tbakta\tCDS\t6001\t6600\t.\t+\t0\tID=g_amb_p;product=x;pseudogene=unitary
c1\tbakta\tCDS\t6001\t6600\t.\t+\t0\tID=g_amb_c;product=x
c1\tbakta\ttRNA\t7001\t7100\t.\t+\t.\tID=t1;pseudogene=unknown
c1\tbakta\tCDS\t8001\t8600\t.\t+\t0\tID=g_uridine;product=tRNA pseudouridine synthase
##FASTA
>c1
ACGT
"""

DFAST_GFF = """\
##gff-version 3
c1\tDFAST\tCDS\t1001\t1600\t.\t+\t0\tID=MGA_1;product=ordinary protein
c1\tDFAST\tCDS\t2001\t2600\t.\t+\t0\tID=MGA_2;product=x
c1\tDFAST\tCDS\t3001\t3600\t.\t+\t0\tID=MGA_3;product=x
c1\tDFAST\tCDS\t5001\t5600\t.\t+\t0\tID=MGA_5;product=x
##FASTA
>c1
ACGT
"""

DFAST_SUMMARY = """\
# locus_id\tlocus_tag\tlocation\tclassification\tinternal_stop\tindel\tref_id\tdescription
MGA_2\tLOCUS_20\tc1:2001..2600(+)\tframeshift\t\t2100..2101(+)\tWP_1\tx
MGA_3\tLOCUS_30\tc1:3001..3600(+)\tpartial\t\t\tWP_2\tx
MGA_5\tLOCUS_50\tc1:5001..5600(+)\tinternal_stop_codon\t5100..5102(+)\t\tWP_3\tx
"""


def setup(tmp):
    b = tmp / "bakta" / "ASM"
    d = tmp / "dfast" / "ASM"
    b.mkdir(parents=True)
    d.mkdir(parents=True)
    (b / "ASM.gff3").write_text(BAKTA_GFF)
    (d / "genome.gff").write_text(DFAST_GFF)
    (d / "pseudogene_summary.tsv").write_text(DFAST_SUMMARY)
    return b, d


def cls_at(caller, run, asm, start1, end1):
    feats = PC.load_features(caller, run, asm)
    idx = PC.Index(feats)
    s, e = start1 - 1, end1
    return PC.classify(idx.overlapping("c1", s, e), idx.fragments("c1", s, e))


def test_bakta_classes():
    with tempfile.TemporaryDirectory() as td:
        b, _ = setup(pathlib.Path(td))
        check(cls_at("bakta", b, "ASM", 1001, 1600) == "intact_cds",
              "Bakta ordinary CDS should be intact_cds")
        check(cls_at("bakta", b, "ASM", 2001, 2600) == "explicit_pseudogene",
              "Bakta pseudogene=unitary should assert a pseudogene")
        check(cls_at("bakta", b, "ASM", 3001, 3600) == "explicit_pseudogene",
              "Bakta pseudo=True should also assert a pseudogene")
        check(cls_at("bakta", b, "ASM", 4001, 4600) == "disrupted",
              "a gene split into two CDS should be disrupted")
        check(cls_at("bakta", b, "ASM", 6001, 6600) == "ambiguous",
              "a pseudogene assertion plus an intact CDS should be ambiguous")
        check(cls_at("bakta", b, "ASM", 12000, 12600) == "absent",
              "a region with no feature should be absent")
        check(cls_at("bakta", b, "ASM", 8001, 8600) == "intact_cds",
              "a pseudouridine synthase is an intact CDS, not a pseudogene")


def test_dfast_classes():
    with tempfile.TemporaryDirectory() as td:
        _, d = setup(pathlib.Path(td))
        check(cls_at("dfast", d, "ASM", 1001, 1600) == "intact_cds",
              "DFAST ordinary CDS should be intact_cds")
        check(cls_at("dfast", d, "ASM", 2001, 2600) == "explicit_pseudogene",
              "DFAST frameshift should assert a pseudogene")
        check(cls_at("dfast", d, "ASM", 5001, 5600) == "explicit_pseudogene",
              "DFAST internal_stop_codon should assert a pseudogene")
        check(cls_at("dfast", d, "ASM", 3001, 3600) == "partial",
              "DFAST partial must be partial, NOT explicit_pseudogene")
        check(cls_at("dfast", d, "ASM", 12000, 12600) == "absent",
              "a region with no feature should be absent")


def test_tRNA_pseudogene_is_not_a_cds_call():
    """Bakta marks some tRNAs pseudogene=unknown; those are not CDS calls."""
    with tempfile.TemporaryDirectory() as td:
        b, _ = setup(pathlib.Path(td))
        feats = PC.load_features("bakta", b, "ASM")
        t = [f for f in feats if f.ident == "t1"]
        check(len(t) == 1 and t[0].kind == PC.OTHER,
              f"a tRNA row must be kind 'other', got "
              f"{t[0].kind if t else 'missing'}")


def test_endpoints():
    """A and B must not be the same thing, and absent must differ between them."""
    check(PC.endpoint_positive("explicit_pseudogene", "A"), "A: pseudogene")
    check(not PC.endpoint_positive("partial", "A"),
          "A must NOT count partial")
    check(not PC.endpoint_positive("absent", "A"),
          "A must NOT count absent")
    check(not PC.endpoint_positive("disrupted", "A"),
          "A must NOT count disrupted")
    for c in ("explicit_pseudogene", "disrupted", "partial", "absent"):
        check(PC.endpoint_positive(c, "B"), f"B: {c}")
    check(not PC.endpoint_positive("intact_cds", "B"), "B: intact is negative")
    check(PC.endpoint_positive("intact_cds", "C"), "C: intact is recovered")
    check(not PC.endpoint_positive("absent", "C"), "C: absent is not recovered")


def test_reciprocal_overlap_rule():
    """The plan's reciprocal rule, and the fragment rule, kept distinct."""
    feats = [PC.Feat("c1", 1000, 1060, "+", PC.PSEUDO, "tiny")]
    idx = PC.Index(feats)
    check(idx.overlapping("c1", 1000, 1600) == [],
          "a 60 bp feature must not match a 600 bp query at 50% reciprocal")
    check(len(idx.fragments("c1", 1000, 1600)) == 1,
          "but it IS a fragment: all of it lies inside the query")
    check(len(idx.any_overlap("c1", 1000, 1600)) == 1,
          "and it overlaps, which is what endpoint C asks")

    idx = PC.Index([PC.Feat("c1", 1000, 1600, "+", PC.PSEUDO, "same")])
    check(len(idx.overlapping("c1", 1000, 1600)) == 1,
          "an identical interval must match")
    # 300 of 600 each way is exactly at the bar
    idx = PC.Index([PC.Feat("c1", 1300, 1900, "+", PC.PSEUDO, "half")])
    check(len(idx.overlapping("c1", 1000, 1600)) == 1,
          "exactly 50% reciprocal must match")


def test_fragmentation_is_detectable():
    """Two half-length fragments must read as disrupted, not absent.

    This is why Index.fragments exists: under the reciprocal rule alone a gene
    split in half has no qualifying feature at all, so `disrupted` would be
    unreachable and endpoint B would miss every fragmentation event.
    """
    halves = [PC.Feat("c1", 4000, 4280, "+", PC.CDS, "a"),
              PC.Feat("c1", 4319, 4600, "+", PC.CDS, "b")]
    idx = PC.Index(halves)
    check(idx.overlapping("c1", 4000, 4600) == [],
          "neither half clears the reciprocal bar, as expected")
    check(PC.classify(idx.overlapping("c1", 4000, 4600),
                      idx.fragments("c1", 4000, 4600)) == "disrupted",
          "two contained CDS fragments must classify as disrupted")
    # one fragment alone is not disruption by the plan's definition
    idx1 = PC.Index([halves[0]])
    check(PC.classify(idx1.overlapping("c1", 4000, 4600),
                      idx1.fragments("c1", 4000, 4600)) == "absent",
          "a single fragment must stay absent, which endpoint B already counts")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    print(f"running {len(tests)} test groups")
    for t in tests:
        print(f"  {t.__name__}")
        t()
    if FAILED:
        print(f"\n{len(FAILED)} failures")
        sys.exit(1)
    print("\nall parser fixture tests pass")


if __name__ == "__main__":
    main()
