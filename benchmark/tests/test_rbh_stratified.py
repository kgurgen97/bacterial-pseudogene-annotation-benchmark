#!/usr/bin/env python3
"""Unit tests for the weighting machinery in analysis_rbh_stratified.py.

An inverse-probability estimator is easy to get silently wrong: a weighted
median that quietly uses a different convention from the unweighted one, a
logistic fit that does not converge, a stabilisation that changes the estimate
rather than just the scale. Each of those would move the headline number without
raising an error, so each is pinned here.

The decisive test is the last one: a simulated dataset where retention depends
on a covariate and the TRUE full-sample gap is known by construction. IPW must
recover it from the retained rows alone, and the naive retained estimate must
not. That is the same check the real analysis performs against S0, run on data
where the answer is unarguable.

    python3 benchmark/tests/test_rbh_stratified.py
"""
import ast
import pathlib
import sys
import types

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
FAILED: list[str] = []


def check(cond, msg):
    if not cond:
        FAILED.append(msg)


def close(a, b, tol=1e-9):
    return a is not None and abs(a - b) <= tol


def load():
    path = ROOT / "analysis_rbh_stratified.py"
    tree = ast.parse(path.read_text())
    keep = [n for n in tree.body
            if isinstance(n, (ast.Import, ast.ImportFrom, ast.FunctionDef,
                              ast.ClassDef))
            or (isinstance(n, ast.Assign)
                and all(isinstance(t, ast.Name) and t.id.isupper()
                        for t in n.targets))]
    m = types.ModuleType("rbh_s")
    m.__dict__["__file__"] = str(path)
    exec(compile(ast.Module(body=keep, type_ignores=[]), str(path), "exec"),
         m.__dict__)
    return m


M = load()


# ------------------------------------------------------------ weighted median
def test_weighted_median_matches_the_unweighted_convention():
    v = [1.0, 2.0, 3.0, 4.0, 5.0]
    check(close(M.wmedian(v, np.ones(5)), 3.0), "odd count, equal weights")
    # even count: the documented 'lower' convention, NOT numpy's midpoint
    check(close(M.wmedian([1.0, 2.0, 3.0, 4.0], np.ones(4)), 2.0),
          "even count must take the lower central value, as documented")
    check(close(M.wmedian([5.0, 1.0, 3.0], [1, 1, 1]), 3.0),
          "input order must not matter")


def test_weights_actually_move_the_median():
    v = [1.0, 2.0, 10.0]
    check(close(M.wmedian(v, [1, 1, 1]), 2.0), "unweighted centre")
    check(close(M.wmedian(v, [1, 1, 100]), 10.0),
          "a dominant weight must pull the median to its value")
    check(close(M.wmedian(v, [100, 1, 1]), 1.0),
          "and the same the other way")


def test_weighted_median_is_scale_invariant():
    v = [1.0, 7.0, 9.0, 11.0]
    w = [2.0, 3.0, 1.0, 4.0]
    check(close(M.wmedian(v, w), M.wmedian(v, np.array(w) * 17.0)),
          "multiplying every weight must change nothing")


def test_weighted_median_handles_degenerate_input():
    check(M.wmedian([], []) is None, "no data must give None")
    check(M.wmedian([1.0, 2.0], [0.0, 0.0]) is None,
          "all-zero weights must give None, not a divide-by-zero")
    # dropping the NaN leaves [1.0, 3.0], and the documented lower convention
    # on an even count takes the LOWER of the two. The first version of this
    # test asserted 3.0, which was the test being wrong about the convention
    # it had just pinned two cases earlier, not the code.
    check(close(M.wmedian([1.0, np.nan, 3.0], [1, 1, 1]), 1.0),
          "a NaN value must be dropped, then the lower convention applies")
    check(close(M.wmedian([1.0, 2.0, np.nan, 3.0, 4.0, 5.0], np.ones(6)), 3.0),
          "dropping a NaN must leave an odd count behaving normally")


# ---------------------------------------------------------- logistic fitting
def test_logistic_recovers_a_known_slope():
    g = np.random.default_rng(7)
    n = 4000
    x = g.normal(size=n)
    eta = -0.5 + 1.5 * x
    y = (g.random(n) < 1 / (1 + np.exp(-eta))).astype(float)
    b = M.logistic_fit(np.column_stack([np.ones(n), x]), y)
    check(b is not None, "the fit must converge")
    check(abs(b[0] + 0.5) < 0.15, f"intercept should be near -0.5, got {b[0]:.3f}")
    check(abs(b[1] - 1.5) < 0.15, f"slope should be near 1.5, got {b[1]:.3f}")


def test_logistic_survives_separable_data():
    """Perfect separation makes the likelihood unbounded; it must not crash."""
    x = np.array([-3.0, -2.0, -1.0, 1.0, 2.0, 3.0])
    y = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
    b = M.logistic_fit(np.column_stack([np.ones(6), x]), y)
    check(b is None or np.all(np.isfinite(b)),
          "separable data must give None or finite coefficients, never NaN")


# ------------------------------------------------------------------ weights
def _cell(n_ret, n_exc, cov_ret, cov_exc):
    import pandas as pd
    rows = [{"hiconf": 1, "log_span": 3.0, "forward_coverage": cov_ret,
             "omega": 0.2} for _ in range(n_ret)]
    rows += [{"hiconf": 0, "log_span": 3.0, "forward_coverage": cov_exc,
              "omega": 0.1} for _ in range(n_exc)]
    return pd.DataFrame(rows)


def test_weights_are_stabilised_to_mean_one():
    g = np.random.default_rng(11)
    import pandas as pd
    n = 400
    cov = g.random(n)
    keep = (g.random(n) < 0.3 + 0.5 * cov).astype(int)
    cell = pd.DataFrame({"hiconf": keep, "log_span": 3.0,
                         "forward_coverage": cov, "omega": 0.1})
    w, d = M.ipw_weights(cell)
    check(w is not None and d["status"] == "ok", f"fit should work: {d}")
    check(close(float(np.mean(w)), 1.0, 1e-9),
          f"stabilised weights must average 1, got {np.mean(w):.6f}")
    check(d["ess"] <= d["n_retained"] + 1e-9,
          "effective sample size cannot exceed the retained count")


def test_weight_cap_is_applied_and_counted():
    check(M.WEIGHT_CAP == 10.0, "the cap must be the declared 10.0")
    g = np.random.default_rng(3)
    import pandas as pd
    n = 600
    cov = g.random(n)
    keep = (g.random(n) < np.clip(0.01 + 0.95 * cov ** 6, 0, 1)).astype(int)
    cell = pd.DataFrame({"hiconf": keep, "log_span": 3.0,
                         "forward_coverage": cov, "omega": 0.1})
    w, d = M.ipw_weights(cell)
    if w is not None:
        check(float(np.max(w)) <= M.WEIGHT_CAP + 1e-9,
              "no weight may exceed the cap")
        check("n_clipped" in d, "the clipped count must be reported")


def test_cells_too_small_to_model_are_refused_not_guessed():
    w, d = M.ipw_weights(_cell(5, 20, 1.0, 1.0))
    check(w is None and d["status"] == "too_few_to_model",
          "too few retained rows must refuse")
    w, d = M.ipw_weights(_cell(50, 1, 1.0, 1.0))
    check(w is None and d["status"] == "too_few_to_model",
          "too few excluded rows must refuse")


# ----------------------------------------------------------------- balance
def test_smd_weighted():
    check(close(M.smd_w([1, 2, 3], [1, 1, 1], [1, 2, 3], [1, 1, 1]), 0.0),
          "identical groups are balanced")
    check(M.smd_w([1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]) == 0.0,
          "no variance must not divide by zero")
    v = M.smd_w([10, 11, 12], [1, 1, 1], [1, 2, 3], [1, 1, 1])
    check(v is not None and v > 5, f"a large shift must show, got {v}")


# ------------------------------------------- the decisive end-to-end check
def test_ipw_recovers_a_known_gap_that_the_naive_estimate_misses():
    """Retention depends on a covariate that also tracks omega.

    Built so the naive retained gap is wrong by a wide margin and IPW, using
    only the retained rows, gets back to the truth. This is the same logic the
    real analysis applies against S0, on data where the answer is not in doubt.
    """
    import pandas as pd
    g = np.random.default_rng(42)
    n = 6000
    rows = []
    for kind, base in (("discordant", 0.14), ("concordant_intact", 0.10)):
        cov = g.random(n)                       # the covariate
        omega = base + 0.10 * cov               # omega rises with it
        # retention rises steeply with the covariate, and much harder on the
        # discordant side: exactly the real pattern
        pr = np.clip((0.05 if kind == "discordant" else 0.60) + 0.9 * cov,
                     0.02, 0.99)
        keep = (g.random(n) < pr).astype(int)
        for c, om, k in zip(cov, omega, keep):
            rows.append({"accession": "A", "kind": kind, "hiconf": int(k),
                         "log_span": 3.0, "forward_coverage": float(c),
                         "omega": float(om)})
    df = pd.DataFrame(rows)

    truth = M.gap_unweighted(df)
    naive = M.gap_unweighted(df, mask=(df.hiconf == 1))
    wd, dd = M.ipw_weights(df[df.kind == "discordant"])
    wc, dc = M.ipw_weights(df[df.kind == "concordant_intact"])
    check(wd is not None and wc is not None, f"both fits must work: {dd} {dc}")
    rd = df[(df.kind == "discordant") & (df.hiconf == 1)]
    rc = df[(df.kind == "concordant_intact") & (df.hiconf == 1)]
    ipw = M.wmedian(rd.omega, wd) - M.wmedian(rc.omega, wc)

    check(close(truth, 0.04, 0.004), f"truth should be ~0.040, got {truth:.4f}")
    check(abs(naive - truth) > 0.010,
          f"the naive retained gap should be badly biased here, "
          f"got {naive:.4f} against {truth:.4f}")
    check(abs(ipw - truth) < abs(naive - truth),
          f"IPW must beat the naive estimate: |{ipw:.4f}-{truth:.4f}| vs "
          f"|{naive:.4f}-{truth:.4f}|")
    check(abs(ipw - truth) < 0.006,
          f"IPW should land within 0.006 of the truth, got {ipw:.4f} "
          f"against {truth:.4f}")


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
    print("\nall stratified-estimator tests pass")


if __name__ == "__main__":
    main()
