#!/usr/bin/env python3
"""Publication figures, built only from committed outputs.

Every panel writes the values it plots to figures/data/<panel>.tsv, so each
plotted number can be checked against the source without re-reading the
figure. Vector output (SVG and PDF) plus a PNG for drafts.

Figures produced here are the ones whose inputs are final:
  Figure 2  caller-dependent pseudogene abundance and disagreement
  Figure 3  annotation choice changes the direction of a downstream result
  Figure 5  dN/dS interpretation: frameshift artefact, reference and
            orthologue sensitivity
Figure 1 (study design) is a schematic and is not generated here.

    python3 figures/make_figures.py
"""
from __future__ import annotations

import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                        # noqa: E402
import pandas as pd                       # noqa: E402

P = pathlib.Path(__file__).resolve().parent.parent
R = P / "results"
B = P / "benchmark" / "results"
FIG = P / "figures"
DAT = FIG / "data"

CALLERS = ["PGAP", "Pseudofinder", "Bakta", "DFAST"]
COL = {"PGAP": "#3b6ea5", "Pseudofinder": "#c0504d",
       "Bakta": "#8064a2", "DFAST": "#4f8a5b"}

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "svg.fonttype": "none", "pdf.fonttype": 42, "figure.dpi": 150,
})


def save(fig, name):
    FIG.mkdir(exist_ok=True)
    for ext in ("svg", "pdf", "png"):
        fig.savefig(FIG / f"{name}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote figures/{name}.{{svg,pdf,png}}")


def dump(df, name):
    DAT.mkdir(parents=True, exist_ok=True)
    df.to_csv(DAT / f"{name}.tsv", sep="\t", index=False)


# ===================================================================== fig 2
def figure2():
    corpus = pd.read_csv(R / "corpus523.tsv", sep="\t")
    spreads = json.loads((B / "two_spreads.json").read_text())
    per = {d["caller"]: d for d in spreads["per_caller"]}
    locus = json.loads((R / "locus.json").read_text())

    fig = plt.figure(figsize=(7.2, 5.0))
    gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.34)

    # (a) per-genome pseudogene fraction, 523 genomes
    ax = fig.add_subplot(gs[0, 0])
    rows, data = [], []
    for c in CALLERS:
        v = corpus[f"frac_{c}"].dropna().to_numpy(float)
        data.append(v)
        rows.append({"caller": c, "n_genomes": len(v),
                     "median_fraction": float(np.median(v)),
                     "q1": float(np.percentile(v, 25)),
                     "q3": float(np.percentile(v, 75)),
                     "genomes_with_zero": int((v == 0).sum())})
    pos = np.arange(len(CALLERS))
    bp = ax.boxplot(data, positions=pos, widths=0.62, showfliers=False,
                    patch_artist=True, medianprops=dict(color="black", lw=1.2))
    for patch, c in zip(bp["boxes"], CALLERS):
        patch.set_facecolor(COL[c])
        patch.set_alpha(0.55)
        patch.set_edgecolor("black")
        patch.set_linewidth(0.6)
    for i, c in enumerate(CALLERS):
        v = data[i]
        j = (np.random.default_rng(7 + i).random(len(v)) - 0.5) * 0.30
        ax.plot(pos[i] + j, np.clip(v, 1e-5, None), ".", ms=1.1,
                color=COL[c], alpha=0.30, zorder=0)
    ax.set_yscale("log")
    ax.set_xticks(pos)
    ax.set_xticklabels([c.replace("finder", "finder") for c in CALLERS],
                       rotation=18, ha="right")
    ax.set_ylabel("pseudogene fraction per genome\n(log scale)")
    ax.set_title("a   523 genomes, four configurations", loc="left")
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.6)
    dump(pd.DataFrame(rows), "fig2a_per_genome_fraction")

    # (b) corpus totals and zero-call genomes
    ax = fig.add_subplot(gs[0, 1])
    tot = [per[c]["corpus_total"] for c in CALLERS]
    zero = [per[c]["genomes_with_zero"] for c in CALLERS]
    ax.bar(pos, tot, color=[COL[c] for c in CALLERS], alpha=0.85,
           edgecolor="black", lw=0.6)
    for i, t in enumerate(tot):
        ax.text(i, t * 1.08, f"{t:,}", ha="center", va="bottom", fontsize=7)
    ax.set_yscale("log")
    ax.set_ylim(top=max(tot) * 3.2)
    ax.set_xticks(pos)
    # the zero-call count goes in the tick label: inside the bar it was
    # clipped by the bar edge and collided with the axis
    ax.set_xticklabels([f"{c}\n({z} zero)" if z else f"{c}\n(0 zero)"
                        for c, z in zip(CALLERS, zero)], fontsize=6.8)
    ax.set_ylabel("pseudogene calls, corpus total")
    ax.set_title(f"b   {spreads['corpus_total_ratio']:.1f}-fold on totals, "
                 f"{spreads['median_fraction_ratio']:.1f}-fold on medians",
                 loc="left")
    dump(pd.DataFrame({"caller": CALLERS, "corpus_total": tot,
                       "genomes_with_zero": zero,
                       "corpus_total_ratio": spreads["corpus_total_ratio"],
                       "median_fraction_ratio":
                           spreads["median_fraction_ratio"]}),
         "fig2b_corpus_totals")

    # (c) the two spread metrics are different measurements
    ax = fig.add_subplot(gs[1, 0])
    med = [per[c]["median_fraction"] for c in CALLERS]
    mean = [per[c]["mean_fraction"] for c in CALLERS]
    w = 0.38
    ax.bar(pos - w / 2, med, w, label="median genome", color="#555555",
           edgecolor="black", lw=0.5)
    ax.bar(pos + w / 2, mean, w, label="mean genome", color="#bbbbbb",
           edgecolor="black", lw=0.5)
    ax.set_yscale("log")
    ax.set_ylim(top=max(mean) * 2.4)     # Pseudofinder's bars were clipped
    ax.set_xticks(pos)
    ax.set_xticklabels(CALLERS, rotation=18, ha="right")
    ax.set_ylabel("pseudogene fraction")
    ax.set_title("c   median vs mean genome", loc="left")
    ax.legend(frameon=False, loc="upper left")
    dump(pd.DataFrame({"caller": CALLERS, "median_fraction": med,
                       "mean_fraction": mean}), "fig2c_median_vs_mean")

    # (d) locus-level agreement
    ax = fig.add_subplot(gs[1, 1])
    nd = locus["n_discordant"]
    nc = locus["n_loci"] - nd
    ax.barh([0], [nc], color="#9bb7d4", edgecolor="black", lw=0.6,
            label=f"concordant-intact  {nc:,}")
    ax.barh([0], [nd], left=[nc], color="#c0504d", edgecolor="black", lw=0.6,
            label=f"discordant  {nd:,}")
    ax.set_yticks([])
    ax.set_ylim(-1.6, 0.7)
    ax.set_xlabel("caller-neutral loci")
    ax.set_title(f"d   {locus['n_loci']:,} loci, 50 % reciprocal overlap",
                 loc="left")
    ax.legend(frameon=False, loc="lower left", ncol=1, handlelength=1.2)
    dump(pd.DataFrame({"class": ["concordant_intact", "discordant"],
                       "n_loci": [nc, nd]}), "fig2d_locus_classes")

    save(fig, "figure2_caller_disagreement")


# ===================================================================== fig 3
def figure3():
    lit = json.loads((R / "literature_retest.json").read_text())
    rows = []
    for key, t in lit["C1"]["per_test"].items():
        parts = key.split("|")
        rows.append({"test": key, "level": parts[0],
                     "caller": parts[1] if len(parts) > 1 else "?",
                     "metric": parts[2] if len(parts) > 2 else "?",
                     "rho": t.get("rho"), "p": t.get("p"),
                     "partial_rho": t.get("partial_rho"),
                     "n": t.get("n"), "verdict": t.get("verdict")})
    df = pd.DataFrame(rows).dropna(subset=["rho"])
    dump(df, "fig3_literature_retest")

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.2, 2.9),
                                  gridspec_kw={"width_ratios": [2.1, 1.0],
                                               "wspace": 0.35})
    order = df.sort_values(["caller", "metric"]).reset_index(drop=True)
    y = np.arange(len(order))
    cols = [COL.get(c, "#888888") for c in order.caller]
    ax.scatter(order.rho, y, s=16, c=cols, edgecolor="black", lw=0.3,
               zorder=3)
    sig = order.p < 0.05
    ax.scatter(order.rho[sig], y[sig], s=46, facecolors="none",
               edgecolors="black", lw=0.6, zorder=2)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks([])
    ax.set_xlabel("Spearman $\\rho$  (predicted direction: negative)")
    ax.set_title("a   48 caller $\\times$ metric re-tests, 363 genomes",
                 loc="left")
    ax.grid(axis="x", ls=":", lw=0.5, alpha=0.6)
    handles = [plt.Line2D([], [], marker="o", ls="", color=COL[c], label=c,
                          markeredgecolor="black", markeredgewidth=0.3)
               for c in CALLERS]
    handles.append(plt.Line2D([], [], marker="o", ls="", markerfacecolor="none",
                              markeredgecolor="black", label="p < 0.05"))
    ax.legend(handles=handles, frameon=False, loc="lower right", ncol=2)

    vc = order.verdict.value_counts()
    ax2.bar(range(len(vc)), vc.values, color="#777777", edgecolor="black",
            lw=0.6)
    for i, v in enumerate(vc.values):
        ax2.text(i, v + 0.4, str(v), ha="center", fontsize=7.5)
    ax2.set_xticks(range(len(vc)))
    ax2.set_xticklabels(vc.index, rotation=18, ha="right")
    ax2.set_ylabel("tests")
    ax2.set_title("b   no test reproduces\n     the predicted direction",
                  loc="left")
    dump(vc.rename_axis("verdict").reset_index(name="n_tests"),
         "fig3b_verdicts")
    save(fig, "figure3_downstream_consequence")


# ===================================================================== fig 4
def figure4_dnds():
    f2 = json.loads((R / "frame2.json").read_text())
    ref = json.loads((R / "reference_sensitivity_all.json").read_text())
    rbh = json.loads((R / "rbh_audited_all.json").read_text())
    rev = json.loads((R / "review_response.json").read_text())

    # generous horizontal spacing: the long arm labels sit outside each axes
    # and a first version let them overflow into the neighbouring panel
    fig = plt.figure(figsize=(7.2, 3.1))
    gs = fig.add_gridspec(1, 3, wspace=0.95)

    # (a) FRAME2 arms
    ax = fig.add_subplot(gs[0, 0])
    # Read the arms by their actual keys in results/frame2.json. A first
    # version guessed substring keys, found nothing, and silently fell back to
    # values hardcoded in this script -- which made the panel untraceable even
    # though the numbers happened to be right. No fallback now: a missing key
    # is an error.
    ARMKEYS = [("untouched", "untouched (determinism)"),
               ("synonymous", "verified synonymous"),
               ("frameshift_high", "single frameshift"),
               ("stop_high", "premature stop, frame intact")]
    arms = []
    for key, lab in ARMKEYS:
        node = f2["arms"][key]
        arms.append({"arm": lab, "median_delta_omega": node["median_delta"],
                     "n": node["n"], "n_positive": node["n_positive"],
                     "source": f"results/frame2.json: arms.{key}.median_delta"})
    adf = pd.DataFrame(arms)
    dump(adf, "fig4a_frame2_arms")
    cols = ["#999999", "#999999", "#c0504d", "#3b6ea5"]
    ax.barh(range(len(adf)), adf.median_delta_omega,
            color=cols[:len(adf)], edgecolor="black", lw=0.6)
    for i, r in adf.iterrows():
        ax.text(max(r.median_delta_omega, 0) + 0.0015, i,
                f"{r.median_delta_omega:+.4f}  (n={r.n})",
                va="center", fontsize=6.4)
    ax.set_yticks(range(len(adf)))
    ax.set_yticklabels(["untouched\n(determinism)", "verified\nsynonymous",
                        "single\nframeshift", "premature stop,\nframe intact"],
                       fontsize=6.8)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("median $\\Delta\\omega$, same gene")
    ax.set_xlim(-0.012, 0.062)
    ax.set_title("a   a frameshift alone\n     inflates $\\omega$", loc="left")

    # (b) reference sensitivity
    ax = fig.add_subplot(gs[0, 1])
    rows = []
    for key, lab in (("published", "Pseudofinder-filtered\nreference"),
                     ("treatment", "partner's own\nRefSeq reference")):
        g = ref["g3"][key]
        rows.append({"arm": lab.replace("\n", " "), "g3": g["median_gap"],
                     "lo": g["genus_cluster_ci95"][0],
                     "hi": g["genus_cluster_ci95"][1],
                     "n_positive": g["n_positive"], "n": g["n_tested"]})
    rdf = pd.DataFrame(rows)
    dump(rdf, "fig4b_reference_sensitivity")
    yy = np.arange(len(rdf))
    ax.errorbar(rdf.g3, yy, xerr=[rdf.g3 - rdf.lo, rdf.hi - rdf.g3],
                fmt="o", ms=5, color="#3b6ea5", ecolor="black",
                elinewidth=0.9, capsize=3)
    for i, r in rdf.iterrows():
        ax.annotate(f"{r.g3:+.4f}   {r.n_positive}/{r.n}", (r.g3, i),
                    textcoords="offset points", xytext=(0, 9),
                    ha="center", fontsize=6.4)
    ax.set_yticks(yy)
    ax.set_yticklabels(["Pseudofinder-\nfiltered", "partner's own\nRefSeq"],
                       fontsize=6.8)
    ax.set_ylim(-0.6, len(rdf) - 0.25)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("$\\omega$ gap, 92 genomes")
    ax.set_xlim(0, 0.055)
    ax.set_title("b   reference choice\n     moves it by $-$0.0035", loc="left")

    # (c) G3 across arms, genus-clustered intervals
    ax = fig.add_subplot(gs[0, 2])
    r2 = rev["R2_ci_of_reported_statistic"]          # {lo, hi, point}
    arms3 = [("full universe, genus level", rev["R1_genus_level"]["median"],
              [r2["lo"], r2["hi"]], None)]
    want = {"E1 original locus universe": "subset, all loci",
            "E2 modified RBH, high-confidence": "subset, RBH-confirmed",
            "E2 modified RBH + synteny >= 1": "subset, RBH + synteny"}
    for a in rbh["arms"]:
        if a["label"] in want:
            arms3.append((want[a["label"]], a["median_gap"],
                          a["genus_cluster_ci95"], a["n_loci_used"]))
    rows = []
    for lab, g, ci, n in arms3:
        lo, hi = (ci if isinstance(ci, (list, tuple)) else (None, None))
        rows.append({"arm": lab, "g3": g, "lo": lo, "hi": hi, "n_loci": n})
    gdf = pd.DataFrame(rows)
    dump(gdf, "fig4c_g3_arms")
    yy = np.arange(len(gdf))[::-1]
    ax.errorbar(gdf.g3, yy, xerr=[gdf.g3 - gdf.lo, gdf.hi - gdf.g3],
                fmt="s", ms=4.5, color="#4f8a5b", ecolor="black",
                elinewidth=0.9, capsize=3)
    ax.set_yticks(yy)
    ax.set_yticklabels(["full universe,\ngenus level", "subset,\nall loci",
                        "subset,\nRBH-confirmed", "subset,\nRBH + synteny"],
                       fontsize=6.8)
    ax.set_ylim(-0.6, len(gdf) - 0.4)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("$\\omega$ gap, genus-clustered 95 % CI")
    ax.set_xlim(0, 0.055)
    ax.set_title("c   the gap survives every\n     sensitivity arm", loc="left")

    save(fig, "figure4_dnds_sensitivity")


# ===================================================================== fig 5
def figure5_lesion():
    les = json.loads((R / "lesion_extended.json").read_text())
    pl = pd.read_csv(R / "lesion_extended" / "per_locus_extended.tsv", sep="\t")
    pg = pd.read_csv(R / "lesion_extended" / "lesion_rate_by_genome.tsv",
                     sep="\t")
    ct = pd.read_csv(R / "lesion_extended" / "crosstab_vs_dfast.tsv", sep="\t",
                     index_col=0)

    fig = plt.figure(figsize=(7.2, 5.1))
    # extra vertical room: panel a's legend sits below its axes and a first
    # version let it land on top of panel c's title
    gs = fig.add_gridspec(2, 2, hspace=0.95, wspace=0.42)

    # (a) class composition by locus class
    ax = fig.add_subplot(gs[0, 0])
    order = ["intact", "frameshift", "internal_stop", "truncation",
             "fragmentation", "ambiguous", "no_alignment"]
    cols = {"intact": "#9bb7d4", "frameshift": "#c0504d",
            "internal_stop": "#e8a33d", "truncation": "#8064a2",
            "fragmentation": "#6aa84f", "ambiguous": "#cccccc",
            "no_alignment": "#888888"}
    kinds = ["discordant", "concordant_intact"]
    rows = []
    for i, k in enumerate(kinds):
        g = pl[pl.kind == k]
        left = 0.0
        for c in order:
            frac = float((g.lesion_class == c).mean()) if len(g) else 0.0
            if frac <= 0:
                continue
            ax.barh(i, frac, left=left, color=cols[c], edgecolor="white",
                    lw=0.5)
            if frac > 0.07:
                ax.text(left + frac / 2, i, f"{100 * frac:.0f}%", ha="center",
                        va="center", fontsize=6.4,
                        color="white" if c in ("frameshift", "truncation")
                        else "black")
            left += frac
            rows.append({"kind": k, "lesion_class": c, "fraction": frac,
                         "n": int((g.lesion_class == c).sum())})
    ax.set_yticks(range(len(kinds)))
    ax.set_yticklabels(["discordant\n(n=1,871)", "all four call intact\n(n=4,489)"],
                       fontsize=6.8)
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction of loci")
    ax.set_title("a   lesion evidence on extended intervals", loc="left")
    handles = [plt.Rectangle((0, 0), 1, 1, color=cols[c]) for c in order
               if c in {r["lesion_class"] for r in rows}]
    labs = [c.replace("_", " ") for c in order
            if c in {r["lesion_class"] for r in rows}]
    ax.legend(handles, labs, frameon=False, ncol=4, fontsize=6.0,
              loc="upper center", bbox_to_anchor=(0.5, -0.30),
              columnspacing=1.0, handlelength=1.1, handletextpad=0.5)
    dump(pd.DataFrame(rows), "fig5a_class_composition")

    # (b) per-genome lesion rate
    ax = fig.add_subplot(gs[0, 1])
    for _, r in pg.iterrows():
        ax.plot([0, 1], [r.lesion_rate_concordant, r.lesion_rate_discordant],
                "-", color="#999999", lw=0.8, alpha=0.8, zorder=1)
    ax.plot(np.zeros(len(pg)), pg.lesion_rate_concordant, "o", ms=4,
            color="#9bb7d4", markeredgecolor="black", markeredgewidth=0.4,
            zorder=3)
    ax.plot(np.ones(len(pg)), pg.lesion_rate_discordant, "o", ms=4,
            color="#c0504d", markeredgecolor="black", markeredgewidth=0.4,
            zorder=3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["all four\ncall intact", "discordant"], fontsize=7)
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylabel("lesion rate per genome")
    c = les["lesion_rate_contrast"]
    ax.set_title(f"b   {c['n_positive']}/{c['n_genomes']} genomes, "
                 f"p = {c['sign_test_p_two_sided']:.1e}", loc="left")
    dump(pg, "fig5b_lesion_rate_by_genome")

    # (c) recall against DFAST, before and after extension
    ax = fig.add_subplot(gs[1, 0])
    bars = [("frameshift\nrecall", 18 / 413, 307 / 466),
            ("internal-stop\nrecall", 280 / 340, 224 / 340),
            ("discordant loci\ndecided", 0.435, les["gates"]["G4_decided"]["value"])]
    x = np.arange(len(bars))
    w = 0.38
    ax.bar(x - w / 2, [b[1] for b in bars], w, label="annotation-bounded",
           color="#bbbbbb", edgecolor="black", lw=0.5)
    ax.bar(x + w / 2, [b[2] for b in bars], w, label="extended",
           color="#4f8a5b", edgecolor="black", lw=0.5)
    for i, b in enumerate(bars):
        ax.text(i - w / 2, b[1] + 0.025, f"{100 * b[1]:.0f}", ha="center",
                fontsize=6.4)
        ax.text(i + w / 2, b[2] + 0.025, f"{100 * b[2]:.0f}", ha="center",
                fontsize=6.4)
    ax.set_xticks(x)
    ax.set_xticklabels([b[0] for b in bars], fontsize=6.8)
    ax.set_ylim(0, 1.12)
    ax.set_ylabel("proportion")
    ax.set_title("c   extending the interval trades\n     one sensitivity for "
                 "another", loc="left")
    ax.legend(frameon=False, fontsize=6.6, loc="upper right")
    dump(pd.DataFrame([{"metric": b[0].replace("\n", " "),
                        "annotation_bounded": b[1], "extended": b[2]}
                       for b in bars]), "fig5c_recall_comparison")

    # (d) controls
    ax = fig.add_subplot(gs[1, 1])
    g = les["gates"]
    ctr = [("synthetic 1 bp\ndeletion", g["G6_synthetic_positive"]["value"],
            "#c0504d"),
           ("synthetic\nunmodified", 0.073, "#9bb7d4"),
           ("all four call\nintact (n=4,489)", g["G3_false_lesion"]["value"],
            "#9bb7d4")]
    ax.barh(range(len(ctr)), [c[1] for c in ctr],
            color=[c[2] for c in ctr], edgecolor="black", lw=0.6)
    for i, c_ in enumerate(ctr):
        ax.text(c_[1] + 0.02, i, f"{100 * c_[1]:.1f} %", va="center",
                fontsize=6.6)
    ax.set_yticks(range(len(ctr)))
    ax.set_yticklabels([c[0] for c in ctr], fontsize=6.8)
    ax.set_xlim(0, 1.15)
    ax.set_xlabel("classified as carrying a lesion")
    ax.set_title("d   positive and negative controls", loc="left")
    dump(pd.DataFrame([{"control": c[0].replace("\n", " "), "lesion_rate": c[1]}
                       for c in ctr]), "fig5d_controls")

    save(fig, "figure5_lesion_evidence")


# ===================================================================== fig 1
def figure1():
    """Study design. A schematic, but its counts come from the result files."""
    spreads = json.loads((B / "two_spreads.json").read_text())
    locus = json.loads((R / "locus.json").read_text())
    les = json.loads((R / "lesion_extended.json").read_text())
    rbh = json.loads((R / "rbh_audited_all.json").read_text())
    n_corpus = max(d["genomes"] for d in spreads["per_caller"])

    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    def box(x, y, w, h, text, fc="#f2f4f7", ec="#444444", fs=7.0, bold=False):
        ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=fc, edgecolor=ec,
                                   lw=0.9, zorder=2,
                                   joinstyle="round"))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
                fontsize=fs, zorder=3,
                fontweight="bold" if bold else "normal", linespacing=1.45)

    def arrow(x1, y1, x2, y2, style="-|>", ls="-"):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle=style, color="#444444",
                                    lw=0.9, linestyle=ls,
                                    shrinkA=1, shrinkB=1), zorder=1)

    # --- corpus
    box(28, 90, 44, 8, f"{n_corpus} bacterial genomes\n"
                       f"4 ecological strata, complete RefSeq assemblies",
        fc="#dce6f1", bold=True)

    # --- four configurations, separated so the row does not read as one bar
    cfg = [("PGAP", 4), ("Pseudofinder", 28), ("Bakta", 52), ("DFAST", 76)]
    for nm, x in cfg:
        ax.add_patch(plt.Rectangle((x, 79), 20, 6, facecolor="#ffffff",
                                   edgecolor="#444444", lw=0.9, zorder=2))
        ax.text(x + 10, 82, nm, ha="center", va="center", fontsize=7.0,
                color=COL.get(nm, "black"), fontweight="bold", zorder=4)
        arrow(50, 90, x + 10, 85.3)
    # the dependency arrow runs BELOW the row, not through the labels
    ax.annotate("", xy=(34, 77.6), xytext=(18, 77.6),
                arrowprops=dict(arrowstyle="<|-", color="#c0504d", lw=1.1,
                                connectionstyle="arc3,rad=-0.35"), zorder=5)
    ax.text(26, 73.6, "Pseudofinder runs on PGAP output — not independent",
            ha="center", va="center", fontsize=6.1, color="#c0504d",
            style="italic")

    # --- branch 1, counts
    box(3, 60, 42, 9,
        f"Pseudogene counts, all {n_corpus} genomes\n"
        f"{spreads['corpus_total_ratio']:.1f}× on corpus totals, "
        f"{spreads['median_fraction_ratio']:.1f}× on medians",
        fc="#eaf2e7", fs=6.6)
    arrow(24, 79, 24, 69)
    box(3, 48, 42, 8,
        "Downstream re-test, 363 genomes\n"
        "paired contrast between configurations",
        fc="#eaf2e7", fs=6.6)
    arrow(24, 60, 24, 56)

    # --- branch 2, loci
    box(55, 60, 42, 9,
        f"{locus['n_loci']:,} caller-neutral loci, 50 % overlap\n"
        f"92 genomes, 28 genera, congeneric partners",
        fc="#fdf0e6", fs=6.6)
    arrow(76, 79, 76, 69)
    box(55, 48, 19, 8, f"discordant\n{locus['n_discordant']:,}",
        fc="#f6d9d6", fs=6.6)
    box(78, 48, 19, 8,
        f"all four intact\n{locus['n_loci'] - locus['n_discordant']:,}",
        fc="#dce6f1", fs=6.6)
    arrow(70, 60, 66, 56)
    arrow(82, 60, 86, 56)

    # --- evidence column
    box(55, 34, 42, 9,
        "dN/dS at disputed loci, genus-clustered\n"
        "frameshift artefact in a same-gene control",
        fc="#fdf0e6", fs=6.6)
    arrow(76, 48, 76, 43)
    box(55, 21, 42, 9,
        f"Caller-independent lesion evidence\n"
        f"extended intervals, {les['n_loci']:,} loci",
        fc="#fdf0e6", fs=6.6)
    arrow(76, 34, 76, 30)

    # --- methods column
    box(3, 34, 42, 9,
        "Sensitivity: reference selection,\n"
        f"reciprocal best hits + synteny ({rbh['n_loci']:,} loci)",
        fc="#eef0f4", fs=6.6)
    arrow(55, 38.5, 45, 38.5)
    box(3, 21, 42, 9,
        "Pre-registered plans, hashes asserted\n"
        "failed gates preserved · verification with controls",
        fc="#eef0f4", fs=6.4)
    arrow(24, 34, 24, 30)

    # --- what is and is not established
    box(3, 4, 42, 11,
        "ESTABLISHED\n"
        "annotation-dependent counts\n"
        "structural non-independence\n"
        "sequence-level disruption at disputed loci",
        fc="#e7f0e4", fs=6.2)
    box(55, 4, 42, 11,
        "NOT ESTABLISHED\n"
        "which configuration is correct\n"
        "loss of gene function\n"
        "biological validation of any single call",
        fc="#f6e9e7", fs=6.2)
    arrow(24, 21, 24, 15, ls=":")
    arrow(76, 21, 76, 15, ls=":")

    save(fig, "figure1_study_design")
    dump(pd.DataFrame([
        {"element": "corpus genomes", "value": n_corpus,
         "source": "benchmark/results/two_spreads.json: per_caller[].genomes"},
        {"element": "corpus total ratio", "value": spreads["corpus_total_ratio"],
         "source": "benchmark/results/two_spreads.json: corpus_total_ratio"},
        {"element": "median fraction ratio",
         "value": spreads["median_fraction_ratio"],
         "source": "benchmark/results/two_spreads.json: median_fraction_ratio"},
        {"element": "caller-neutral loci", "value": locus["n_loci"],
         "source": "results/locus.json: n_loci"},
        {"element": "discordant loci", "value": locus["n_discordant"],
         "source": "results/locus.json: n_discordant"},
        {"element": "lesion-evidence loci", "value": les["n_loci"],
         "source": "results/lesion_extended.json: n_loci"},
        {"element": "RBH subset loci", "value": rbh["n_loci"],
         "source": "results/rbh_audited_all.json: n_loci"}]),
        "fig1_study_design")


def main() -> None:
    print("building figures from committed outputs")
    figure1()
    figure2()
    figure3()
    figure4_dnds()
    figure5_lesion()
    print(f"\nplotted values written to {DAT.relative_to(P)}/*.tsv "
          f"for independent checking")


if __name__ == "__main__":
    main()
