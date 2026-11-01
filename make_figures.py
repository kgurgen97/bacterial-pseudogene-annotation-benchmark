#!/usr/bin/env python3
"""The four manuscript figures, in journal style. Vector PDF plus 600 dpi PNG.

Every value is read from a file in results/. Style matches the companion FunAqua
paper so one author's figures look like one author's figures: Liberation Sans at
6-7 pt, no titles or explanatory text inside the panels, panel letters bold
outside the axes, left and bottom spines only, thin marks.

Palette: blue #2a78d6 and orange #eb6834, the two categorical slots validated for
that paper, plus neutral greys.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

P = Path.home() / "research/projects/n2-pseudogene-caller-bias"
R, F = P / "results", P / "figures"
F.mkdir(exist_ok=True)

BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, MUTED = "#111111", "#4a4a4a", "#8a8a8a"
GRID, SPINE, SURF = "#e8e8e4", "#4a4a4a", "#ffffff"
MM = 1 / 25.4
W1, W2 = 89 * MM, 183 * MM

FAM = "Liberation Sans" if any(f.name == "Liberation Sans" for f in fm.fontManager.ttflist) \
      else "DejaVu Sans"

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": [FAM],
    "font.size": 6.5, "axes.labelsize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
    "axes.edgecolor": SPINE, "axes.linewidth": 0.6,
    "axes.labelcolor": INK, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.color": SPINE, "ytick.color": SPINE,
    "xtick.labelcolor": INK, "ytick.labelcolor": INK,
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.size": 2.4, "ytick.major.size": 2.4,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "axes.grid": False, "grid.color": GRID, "grid.linewidth": 0.5,
    "axes.axisbelow": True,
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.015,
    "lines.linewidth": 0.9, "lines.markersize": 2.8,
    "legend.frameon": False, "legend.handlelength": 1.2,
    "legend.handletextpad": 0.45, "legend.borderpad": 0.0,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def panel(ax, letter, dx=-0.17, dy=1.03):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=8,
            fontweight="bold", va="bottom", ha="left", color=INK)


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(F / f"{name}.{ext}", dpi=600 if ext == "png" else None)
    plt.close(fig)
    print(f"  {name}.pdf / .png")


d = pd.read_csv(R / "corpus523.tsv", sep="\t")
d["PGAP"] = d.pgap_pseudo / (d.pgap_pseudo + d.pgap_cds)
d["Pseudofinder"] = d.pf_pseudo / (d.pf_pseudo + d.pf_intact)
d["DFAST"] = d.dfast_disrupted / (d.dfast_disrupted + d.dfast_cds)
d["Bakta"] = d.bakta_pseudo / (d.bakta_pseudo + d.bakta_cds)
CALLERS = ["Bakta", "PGAP", "DFAST", "Pseudofinder"]        # by median, ascending
cfg = json.load((R / "pf_config.json").open())
c3 = json.load((R / "c3_strata.json").open())

meds = {c: float(np.nanmedian(d[c])) for c in CALLERS}
fold = max(meds.values()) / min(meds.values())
assert abs(fold - cfg["between_caller_fold"]) < 0.5, \
    f"figure computes {fold:.1f}x, pf_config.json says {cfg['between_caller_fold']:.1f}x"
print(f"  between-caller fold {fold:.1f}x, agrees with pf_config.json")

# ---------------------------------------------------------------- Figure 1
fig, ax = plt.subplots(figsize=(W1, W1 * 0.78))
rng = np.random.default_rng(7)
for i, c in enumerate(CALLERS):
    v = d[c].to_numpy()
    v = v[np.isfinite(v) & (v > 0)]
    ax.scatter(np.full(len(v), i) + rng.uniform(-.16, .16, len(v)), v,
               s=1.4, color=MUTED, alpha=.35, linewidths=0, zorder=1)
    q1, q3 = np.percentile(v, [25, 75])
    ax.plot([i, i], [q1, q3], color=INK, linewidth=1.2, zorder=2,
            solid_capstyle="butt")
    ax.plot([i - .26, i + .26], [meds[c]] * 2, color=BLUE, linewidth=1.8, zorder=3)
ax.set_yscale("log")
ax.set_xticks(range(len(CALLERS)))
ax.set_xticklabels(CALLERS)
ax.set_ylabel("pseudogene fraction")
ax.set_xlim(-.6, len(CALLERS) - .4)
save(fig, "fig1_caller_spread")

# ---------------------------------------------------------------- Figure 2
fig, axes = plt.subplots(1, 3, figsize=(W2, W2 * 0.33), sharex=True, sharey=True)
for ax, c in zip(axes, ["Pseudofinder", "DFAST", "Bakta"]):
    m = np.isfinite(d.PGAP) & np.isfinite(d[c]) & (d.PGAP > 0) & (d[c] > 0)
    ax.scatter(d.PGAP[m], d[c][m], s=2.2, color=BLUE, alpha=.45, linewidths=0)
    lim = [5e-5, 0.6]
    ax.plot(lim, lim, color=INK2, linewidth=0.7, linestyle=(0, (3, 1.6)))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(*lim); ax.set_ylim(*lim)
    ax.set_xlabel("PGAP fraction")
    ax.set_ylabel(f"{c} fraction" if c == "Pseudofinder" else "")
    n_above = int((d[c][m] > d.PGAP[m]).sum())
    ax.text(.04, .93, f"{c}\n{n_above} of {int(m.sum())} above the line",
            transform=ax.transAxes, fontsize=5.8, va="top", color=INK2)
for a, L in zip(axes, "abc"):
    panel(a, L, dx=-0.20)
fig.subplots_adjust(wspace=0.16)
save(fig, "fig2_pairwise")

# ---------------------------------------------------------------- Figure 3
S = ["S4_free_living_breadth", "S1_extreme_obligate",
     "S2_obligate_intracellular", "S3_recent_restriction"]
LAB = ["S4\nfree-living", "S1\nextreme\nobligate", "S2\nobligate\nintracellular",
       "S3\nrecent\nrestriction"]
fig, (axa, axb) = plt.subplots(1, 2, figsize=(W2, W2 * 0.40))

STY = {"PGAP": (BLUE, "o", "solid"), "Pseudofinder": (ORANGE, "s", "solid"),
       "DFAST": (BLUE, "^", (0, (3.4, 1.4))), "Bakta": (ORANGE, "D", (0, (1.2, 1.2)))}
for nm, (col, mk, ls) in STY.items():
    med = c3["per_caller"][nm]["medians"]
    base = med["S4_free_living_breadth"]
    y = [med[s] / base for s in S]
    axa.plot(range(4), y, color=col, marker=mk, linestyle=ls, markersize=3.2,
             markerfacecolor="none", markeredgewidth=0.8, label=nm)
axa.axhline(1.0, color=MUTED, linewidth=0.6)
axa.set_xticks(range(4)); axa.set_xticklabels(LAB)
axa.set_ylabel("median fraction, relative to S4")
axa.legend(loc="upper left", ncol=2)

rbc, lo, hi, names = [], [], [], []
for nm in ["Pseudofinder", "Bakta", "DFAST", "PGAP"]:
    v = c3["per_caller"][nm].get("S1_vs_S4")
    if not v:
        continue
    names.append(nm); rbc.append(v["rbc"])
    g1 = d[d.stratum == "S1_extreme_obligate"][nm].dropna().to_numpy()
    g4 = d[d.stratum == "S4_free_living_breadth"][nm].dropna().to_numpy()
    # Percentile bootstrap on the rank-biserial. The stored json has the point
    # estimate and the p value but no interval, and an interval is what shows
    # that these callers do not merely differ, they exclude each other.
    bs = []
    for _ in range(2000):
        a = rng.choice(g1, len(g1), replace=True)
        b = rng.choice(g4, len(g4), replace=True)
        u = float((a[:, None] > b[None, :]).sum())
        bs.append(2 * u / (len(a) * len(b)) - 1)
    lo.append(np.percentile(bs, 2.5)); hi.append(np.percentile(bs, 97.5))
yy = range(len(names))
for i, nm in enumerate(names):
    col = ORANGE if rbc[i] > 0 else BLUE
    axb.plot([lo[i], hi[i]], [i, i], color=col, linewidth=1.1)
    axb.plot([rbc[i]], [i], marker="o", markersize=3.6, color=col,
             markerfacecolor="none", markeredgewidth=1.1)
axb.axvline(0, color=INK, linewidth=0.8)
axb.set_yticks(list(yy)); axb.set_yticklabels(names)
axb.set_xlabel("rank-biserial, S1 against S4")
axb.set_ylim(-0.6, len(names) - 0.4)
axb.tick_params(axis="y", length=0)
for a, L in zip((axa, axb), "ab"):
    panel(a, L, dx=-0.18)
fig.subplots_adjust(wspace=0.34)
save(fig, "fig3_gradient_and_sign")

# ---------------------------------------------------------------- Figure 4
fig, ax = plt.subplots(figsize=(W1, W1 * 0.70))
order = CALLERS + ["Pseudofinder\nconfig B"]
vals = [meds[c] for c in CALLERS] + [cfg["medians"]["frac_PF_B"]]
cols = [MUTED] * len(CALLERS) + [ORANGE]
cols[CALLERS.index("Pseudofinder")] = ORANGE
ax.bar(range(len(order)), vals, color=cols, width=.62, edgecolor="none")
ax.set_yscale("log")
ax.set_xticks(range(len(order)))
ax.set_xticklabels(order, fontsize=5.6)
ax.set_ylabel("median pseudogene fraction")
i_pf, i_b = CALLERS.index("Pseudofinder"), len(order) - 1
top = max(vals) * 2.4
ax.plot([i_pf, i_pf, i_b, i_b],
        [vals[i_pf] * 1.3, top, top, vals[i_b] * 1.3], color=INK, linewidth=0.7)
ax.text((i_pf + i_b) / 2, top * 1.12, f"{cfg['within_caller_fold']:.1f}x",
        ha="center", fontsize=6, color=INK)
lo_i = int(np.argmin(vals[:len(CALLERS)]))
top2 = top * 3.4
ax.plot([lo_i, lo_i, i_pf, i_pf],
        [vals[lo_i] * 1.3, top2, top2, top * 1.6], color=INK, linewidth=0.7)
ax.text((lo_i + i_pf) / 2, top2 * 1.12, f"{fold:.0f}x",
        ha="center", fontsize=6, color=INK)
ax.set_ylim(min(vals) / 3, top2 * 4)
save(fig, "fig4_caller_vs_config")

# ---------------------------------------------------------------- figure 5
# Does selection back the calls? Panel a is the per-caller separation, panel b
# the sharpest contrast in it: Pseudofinder's two distributions come apart and
# DFAST's do not.
dn = pd.read_csv(R / "dnds_pairs.tsv", sep="\t")
dj = json.load((R / "dnds.json").open())
pc = dj["per_caller"]
for c in CALLERS:
    got = float(np.median(dn[(dn.caller == c) & (dn.call == "pseudo")].omega))
    assert abs(got - pc[c]["median_pseudo"]) < 1e-9, \
        f"figure computes {got:.4f} for {c}, dnds.json says {pc[c]['median_pseudo']:.4f}"
print(f"  dN/dS medians agree with dnds.json for all {len(CALLERS)} callers")

fig, (a, b) = plt.subplots(1, 2, figsize=(W2, 62 * MM))

panel(a, "a", dx=-0.22)
y = np.arange(len(CALLERS))
for i, c in enumerate(CALLERS):
    lo, hi = pc[c]["median_intact"], pc[c]["median_pseudo"]
    a.plot([lo, hi], [i, i], color=MUTED, linewidth=1.0, zorder=1,
           solid_capstyle="round")
    a.scatter([lo], [i], s=22, color="white", edgecolor=SPINE,
              linewidth=0.7, zorder=2)
    a.scatter([hi], [i], s=22, color=ORANGE if pc[c]["p"] < 0.05 else "white",
              edgecolor=ORANGE if pc[c]["p"] < 0.05 else SPINE,
              linewidth=0.7, zorder=2)
    lab = f"+{pc[c]['gap']:.3f}"
    if pc[c]["p"] >= 0.05:
        lab += f", p = {pc[c]['p']:.2f}"
    a.text(hi + 0.004, i, lab, va="center", fontsize=5.8,
           color=INK if pc[c]["p"] < 0.05 else INK2)
a.set_yticks(y)
a.set_yticklabels(CALLERS, fontsize=6)
a.set_ylim(-0.6, len(CALLERS) - 0.4)
a.set_xlim(0.085, 0.152)
a.set_xlabel("median dN/dS")
a.scatter([], [], s=22, color="white", edgecolor=SPINE, linewidth=0.7,
          label="intact controls")
a.scatter([], [], s=22, color=ORANGE, edgecolor=ORANGE, linewidth=0.7,
          label="pseudogene calls")
a.legend(frameon=False, loc="lower right", handletextpad=0.3, labelspacing=0.3)

panel(b, "b", dx=-0.17)
for c, col in (("Pseudofinder", ORANGE), ("DFAST", BLUE)):
    for call, style in (("pseudo", "-"), ("intact", (0, (2.2, 1.4)))):
        v = np.sort(dn[(dn.caller == c) & (dn.call == call)].omega.values)
        b.plot(v, np.arange(1, len(v) + 1) / len(v), linestyle=style,
               color=col, linewidth=0.9)
b.set_xlim(0, 0.6)
b.set_xlabel("dN/dS")
b.set_ylabel("cumulative fraction of sequences")
b.plot([], [], linestyle="-", color=INK2, linewidth=0.9, label="pseudogene calls")
b.plot([], [], linestyle=(0, (2.2, 1.4)), color=INK2, linewidth=0.9,
       label="intact controls")
b.legend(frameon=False, loc="lower right", handletextpad=0.6, labelspacing=0.3)
b.text(0.41, 0.42, "Pseudofinder", fontsize=6, color=ORANGE)
b.text(0.41, 0.33, "DFAST", fontsize=6, color=BLUE)

save(fig, "fig5_dnds")

print(f"\n{len(list(F.glob('fig*.png')))} figures in figures/, set in {FAM}")
