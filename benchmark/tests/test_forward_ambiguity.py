#!/usr/bin/env python3
"""Unit tests for analysis_forward_ambiguity.py.

The defect this analysis fixes was a POOLED rate hiding a sevenfold class
difference, so the tests concentrate on the thing that went wrong: a summary
statistic that is right on average and wrong per class. They also pin the FASTA
reader, since a reader that silently drops or mangles a sequence would make the
paralogy test report whatever it liked.

    python3 benchmark/tests/test_forward_ambiguity.py
"""
import ast
import pathlib
import sys
import tempfile
import types

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]
FAILED: list[str] = []


def check(cond, msg):
    if not cond:
        FAILED.append(msg)


def close(a, b, tol=1e-9):
    return a is not None and abs(a - b) <= tol


def load():
    path = ROOT / "analysis_forward_ambiguity.py"
    tree = ast.parse(path.read_text())
    keep = [n for n in tree.body
            if isinstance(n, (ast.Import, ast.ImportFrom, ast.FunctionDef,
                              ast.ClassDef))
            or (isinstance(n, ast.Assign)
                and all(isinstance(t, ast.Name) and t.id.isupper()
                        for t in n.targets))]
    m = types.ModuleType("fa")
    m.__dict__["__file__"] = str(path)
    exec(compile(ast.Module(body=keep, type_ignores=[]), str(path), "exec"),
         m.__dict__)
    return m


M = load()


def _df(n_dis, amb_dis, n_con, amb_con, acc="A", genus="G"):
    rows = []
    for i in range(n_dis):
        rows.append({"accession": acc, "genus": genus, "kind": "discordant",
                     "forward_ambiguous": 1 if i < amb_dis else 0})
    for i in range(n_con):
        rows.append({"accession": acc, "genus": genus,
                     "kind": "concordant_intact",
                     "forward_ambiguous": 1 if i < amb_con else 0})
    return pd.DataFrame(rows)


# ------------------------------------------------- the defect, reproduced
def test_pooled_rate_hides_the_class_difference():
    """The real shape: 1871 discordant at 13.9%, 4489 baseline at 1.9%.

    Pooled this is about 5.4% and passes a 10% threshold; per class it fails.
    The test asserts both, because the point is that one number was hiding the
    other.
    """
    df = _df(1871, 260, 4489, 85)
    t = M.rate_table(df)
    dis = float(t.rate_discordant.iloc[0])
    con = float(t.rate_concordant_intact.iloc[0])
    pooled = float(pd.to_numeric(df.forward_ambiguous).mean())
    check(abs(dis - 0.139) < 0.002, f"discordant rate should be ~13.9%, {dis}")
    check(abs(con - 0.019) < 0.002, f"baseline rate should be ~1.9%, {con}")
    check(pooled <= M.PROXY_THRESHOLD,
          f"pooled {pooled:.3f} must pass the old 10% threshold")
    check(max(dis, con) > M.PROXY_THRESHOLD,
          "the per-class maximum must FAIL the same threshold")
    check(abs(float(t.ratio.iloc[0]) - dis / con) < 1e-9,
          "the reported ratio must be discordant over baseline")
    check(close(float(t.difference_pp.iloc[0]), 100 * (dis - con), 1e-6),
          "the difference must be in percentage points")


def test_rate_table_groups_without_mixing_classes():
    df = pd.concat([_df(10, 5, 10, 0, acc="A", genus="GA"),
                    _df(10, 0, 10, 5, acc="B", genus="GB")])
    t = M.rate_table(df, "accession").set_index("key")
    check(close(float(t.loc["A"].rate_discordant), 0.5),
          "genome A's discordant rate must be 0.5")
    check(close(float(t.loc["A"].rate_concordant_intact), 0.0),
          "genome A's baseline rate must be 0.0")
    check(close(float(t.loc["B"].rate_discordant), 0.0),
          "genome B must not inherit A's rate")
    check(close(float(t.loc["B"].rate_concordant_intact), 0.5),
          "genome B's baseline rate must be 0.5")
    # pooling these two genomes would give 0.25 in both classes and hide both
    allt = M.rate_table(df)
    check(close(float(allt.rate_discordant.iloc[0]), 0.25),
          "the pooled view genuinely loses the per-genome signal")


def test_rate_table_handles_an_empty_class():
    df = _df(5, 2, 0, 0)
    t = M.rate_table(df)
    check(close(float(t.rate_discordant.iloc[0]), 0.4), "discordant still works")
    check(pd.isna(t.rate_concordant_intact.iloc[0]),
          "an absent class must be None, not zero")
    check(pd.isna(t.difference_pp.iloc[0]),
          "no difference can be computed against an absent class")
    check(pd.isna(t.ratio.iloc[0]), "and no ratio either")


def test_rate_table_tolerates_missing_flags():
    df = _df(4, 2, 4, 0)
    df.loc[0, "forward_ambiguous"] = np.nan
    t = M.rate_table(df)
    r = float(t.rate_discordant.iloc[0])
    check(close(r, 1 / 3), f"a NaN flag must be skipped, not counted: {r}")


# ------------------------------------------------------------ fasta reader
def test_read_fasta():
    with tempfile.TemporaryDirectory() as td:
        p = pathlib.Path(td) / "x.faa"
        p.write_text(">a desc here\nMKV\nLLA\n>b\nMMM\n")
        got = M.read_fasta(p)
    check(got == {"a": "MKVLLA", "b": "MMM"},
          f"multi-line sequences must be joined and ids split on space: {got}")


def test_read_fasta_empty_and_trailing():
    with tempfile.TemporaryDirectory() as td:
        p = pathlib.Path(td) / "x.faa"
        p.write_text(">only\nMK\n")
        check(M.read_fasta(p) == {"only": "MK"},
              "the last record must not be dropped")
        p.write_text("")
        check(M.read_fasta(p) == {}, "an empty file must give an empty dict")


# ---------------------------------------------- the borrowed frozen machinery
def test_frozen_machinery_is_reused_not_reimplemented():
    check(hasattr(M.S, "wmedian") and hasattr(M.S, "ipw_weights")
          and hasattr(M.S, "combine") and hasattr(M.S, "logistic_fit"),
          "the weighting machinery must come from analysis_rbh_stratified.py")
    check(close(M.S.wmedian([1.0, 2.0, 3.0], [1, 1, 1]), 2.0),
          "the borrowed weighted median must behave as its own tests pin it")
    check(M.S.PLAN_SHA ==
          "1fadc6db0ee4310b8004a9c488deb262d9d0815687423200edafa3e7fc77d4ef",
          "the borrowed module must still be the one whose plan is frozen")


def test_declared_constants_match_the_plan():
    check(M.CONSEQUENCE_TOL == 0.005, "the plan declares a 0.005 tolerance")
    check(M.PROXY_THRESHOLD == 0.10, "the old proxy threshold was 10%")
    check(M.PARALOGUE_REL_BITSCORE == 0.20,
          "the plan declares a 20% relative bitscore for homology")
    check(M.IPW_COVS == ["log_span", "forward_bitscore", "log_candidates"],
          f"the plan declares those three covariates, got {M.IPW_COVS}")


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
    print("\nall forward-ambiguity tests pass")


if __name__ == "__main__":
    main()
