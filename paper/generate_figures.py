"""Generate all 5 figures for combined_paper_v3a at NeurIPS quality."""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
DRUG_CSV = ROOT.parent / "drug-perturbation-geometry" / "zenodo_v1" / "drug_instability_8949.csv"
FIG5_SOURCE = ROOT / "results" / "fig5_hdac_source" / "fig5_hdac_source.json"
OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(exist_ok=True)

# ── NeurIPS style ──────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8,
    "axes.labelsize": 9,
    "axes.titlesize": 10,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "lines.linewidth": 1.2,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# ── Palette ────────────────────────────────────────────────────────
C_DI = "#2166AC"       # direction instability — strong blue
C_FRECH = "#B2182B"    # Fréchet variance — brick red
C_MAGCV = "#878787"    # magnitude CV — neutral gray
C_JACC = "#4DAF4A"     # Jaccard — green
C_NHR = "#E66101"      # nuclear hormone receptor — warm orange
C_MACH = "#5E3C99"     # constitutive machinery — deep purple
C_LIGHT = "#F7F7F7"    # background accents

K_COLORS = {5: "#FDB863", 10: "#E66101", 13: "#B2182B", 20: "#762A83", 40: "#2166AC"}


# ═══════════════════════════════════════════════════════════════════
# FIGURE 1 — Conceptual schematic: the magnitude confound
# ═══════════════════════════════════════════════════════════════════
def fig1_schematic():
    fig, ax = plt.subplots(figsize=(5.5, 2.8))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5)
    ax.axis("off")

    # --- Boxes ---
    box_kw = dict(boxstyle="round,pad=0.35", linewidth=0.8)

    def draw_box(x, y, w, h, text, color, fc=None, fontsize=8, bold=False):
        fc = fc or (color + "18")
        rect = mpatches.FancyBboxPatch(
            (x - w/2, y - h/2), w, h,
            boxstyle="round,pad=0.2", facecolor=fc,
            edgecolor=color, linewidth=0.9,
        )
        ax.add_patch(rect)
        weight = "bold" if bold else "normal"
        ax.text(x, y, text, ha="center", va="center", fontsize=fontsize,
                fontweight=weight, color=color,
                path_effects=[pe.withStroke(linewidth=2, foreground="white")])
        return (x, y)

    def arrow(x1, y1, x2, y2, color="0.35", style="-|>", lw=1.0, ls="-"):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                     arrowprops=dict(arrowstyle=style, color=color, lw=lw, ls=ls))

    def label_arrow(x1, y1, x2, y2, text, color="0.35", above=True, lw=1.0, ls="-"):
        arrow(x1, y1, x2, y2, color=color, lw=lw, ls=ls)
        mx, my = (x1+x2)/2, (y1+y2)/2
        offset = 0.18 if above else -0.18
        ax.text(mx, my + offset, text, ha="center", va="center",
                fontsize=6.5, color=color, style="italic")

    # Row 1: the confound path (top)
    draw_box(1.5, 4.0, 2.4, 0.7, "Effect\nmagnitude", "#333333", bold=True)
    draw_box(5.0, 4.0, 2.4, 0.7, "Variance /\nsubspace size", C_FRECH)
    draw_box(8.5, 4.0, 2.2, 0.7, "Spurious\n'divergence'", C_FRECH)

    label_arrow(2.75, 4.0, 3.75, 4.0, "inflates", color=C_FRECH, lw=1.2)
    label_arrow(6.25, 4.0, 7.35, 4.0, "reports", color=C_FRECH, lw=1.2)

    # Row 2: the DI path (bottom)
    draw_box(1.5, 1.5, 2.4, 0.7, "Effect\nmagnitude", "#333333", bold=True)
    draw_box(5.0, 1.5, 2.4, 0.7, "Cosine\nnormalization", C_DI)
    draw_box(8.5, 1.5, 2.2, 0.7, "Direction\ninstability", C_DI)

    # Blocked arrow
    arrow(2.75, 1.5, 3.55, 1.5, color=C_DI, lw=1.2)
    ax.plot([3.6, 3.6], [1.2, 1.8], color=C_DI, lw=2.5, zorder=5)
    ax.text(3.35, 1.85, "blocked", ha="center", va="bottom", fontsize=6.5,
            color=C_DI, style="italic")

    label_arrow(6.25, 1.5, 7.35, 1.5, "measures\ndirection only", color=C_DI, lw=1.2)

    # Biology feeding into both
    draw_box(5.0, 2.75, 2.0, 0.55, "Biology", "#2E7D32", bold=True, fontsize=8)
    arrow(5.0, 3.05, 5.0, 3.6, color="#2E7D32", lw=0.8, ls="--")
    arrow(5.0, 2.45, 5.0, 1.9, color="#2E7D32", lw=0.8, ls="--")

    # Labels
    ax.text(0.15, 4.0, "A", fontsize=11, fontweight="bold", va="center",
            color=C_FRECH)
    ax.text(0.15, 1.5, "B", fontsize=11, fontweight="bold", va="center",
            color=C_DI)
    ax.text(0.15, 4.45, "Confounded path", fontsize=7, color=C_FRECH, va="bottom")
    ax.text(0.15, 1.95, "Direction instability", fontsize=7, color=C_DI, va="bottom")

    fig.savefig(OUT / "fig1_schematic.pdf")
    fig.savefig(OUT / "fig1_schematic.png")
    plt.close(fig)
    print("  fig1_schematic done")


# ═══════════════════════════════════════════════════════════════════
# FIGURE 2 — LOO benchmark: AUROC comparison across metrics
# ═══════════════════════════════════════════════════════════════════
def fig2_benchmark():
    with open(RESULTS / "13_lincs_benchmark" / "bootstrap_auroc_cis.json") as f:
        cis = json.load(f)

    jacc = cis["jaccard_aurocs"]

    metrics = ["direction_instability", "frechet_variance", "magnitude_cv"]
    labels = ["Direction\ninstability", "Frechet\nvariance", "Magnitude\nCV"]
    colors = [C_DI, C_FRECH, C_MAGCV]

    points = [jacc[m]["point"] for m in metrics]
    ci_lo = [jacc[m]["ci_lo"] for m in metrics]
    ci_hi = [jacc[m]["ci_hi"] for m in metrics]
    errs_lo = [p - lo for p, lo in zip(points, ci_lo)]
    errs_hi = [hi - p for p, hi in zip(points, ci_hi)]

    fig, ax = plt.subplots(figsize=(3.2, 2.8))

    x = np.arange(len(metrics))
    for i in range(len(metrics)):
        ax.bar(x[i], points[i], width=0.55, color=colors[i], alpha=0.85,
               edgecolor="white", linewidth=0.5, zorder=3)
        ax.errorbar(x[i], points[i], yerr=[[errs_lo[i]], [errs_hi[i]]],
                    fmt="none", ecolor="0.2", capsize=3, capthick=0.8,
                    elinewidth=0.8, zorder=4)
        ax.text(x[i], points[i] + errs_hi[i] + 0.015,
                f"{points[i]:.3f}", ha="center", va="bottom", fontsize=7,
                fontweight="bold", color=colors[i])

    ax.axhline(0.5, color="0.7", ls="--", lw=0.5, zorder=1)

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("AUROC (Jaccard outcome)")
    ax.set_ylim(0.25, 1.02)

    fig.savefig(OUT / "fig2_benchmark.pdf")
    fig.savefig(OUT / "fig2_benchmark.png")
    plt.close(fig)
    print("  fig2_benchmark done")


# ═══════════════════════════════════════════════════════════════════
# FIGURE 3 — Positive control: NHR vs machinery DI distributions
# ═══════════════════════════════════════════════════════════════════
def fig3_positive_control():
    sys.path.insert(0, str(ROOT / "experiments"))
    from positive_control_validation import (
        is_nuclear_receptor, is_constitutive_machinery,
    )

    df = pd.read_csv(DRUG_CSV)
    nhr_mask = df["moa"].apply(is_nuclear_receptor)
    mach_mask = df["moa"].apply(is_constitutive_machinery)

    di_nhr = df.loc[nhr_mask, "direction_instability"].values
    di_mach = df.loc[mach_mask, "direction_instability"].values

    fig, ax = plt.subplots(figsize=(3.2, 2.8))

    rng = np.random.default_rng(0)
    jitter_mach = rng.uniform(-0.15, 0.15, len(di_mach))
    jitter_nhr = rng.uniform(-0.15, 0.15, len(di_nhr))

    ax.scatter(jitter_mach, di_mach, s=12, color=C_MACH, alpha=0.5,
               edgecolors="none", zorder=3)
    ax.scatter(1 + jitter_nhr, di_nhr, s=12, color=C_NHR, alpha=0.5,
               edgecolors="none", zorder=3)

    bp_kw = dict(widths=0.35, showfliers=False, zorder=4,
                 medianprops=dict(color="white", lw=1.5),
                 boxprops=dict(lw=0.8), whiskerprops=dict(lw=0.6),
                 capprops=dict(lw=0.6))
    bp1 = ax.boxplot([di_mach], positions=[0], patch_artist=True, **bp_kw)
    bp2 = ax.boxplot([di_nhr], positions=[1], patch_artist=True, **bp_kw)
    bp1["boxes"][0].set_facecolor(C_MACH)
    bp1["boxes"][0].set_alpha(0.4)
    bp1["boxes"][0].set_edgecolor(C_MACH)
    bp2["boxes"][0].set_facecolor(C_NHR)
    bp2["boxes"][0].set_alpha(0.4)
    bp2["boxes"][0].set_edgecolor(C_NHR)

    mean_diff = np.mean(di_nhr) - np.mean(di_mach)
    ybar = max(np.max(di_nhr), np.max(di_mach)) + 0.04
    ax.plot([0, 1], [ybar, ybar], color="0.3", lw=0.6)
    ax.plot([0, 0], [ybar - 0.01, ybar], color="0.3", lw=0.6)
    ax.plot([1, 1], [ybar - 0.01, ybar], color="0.3", lw=0.6)
    ax.text(0.5, ybar + 0.015,
            f"+0.139 [0.107, 0.173]",
            ha="center", va="bottom", fontsize=7, fontweight="bold")

    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"Constitutive\nmachinery\n(n = {len(di_mach)})",
                        f"Nuclear hormone\nreceptor\n(n = {len(di_nhr)})"],
                        fontsize=7.5)
    ax.set_ylabel("Direction instability (D)")
    ax.set_ylim(0.45, ybar + 0.1)

    fig.savefig(OUT / "fig3_positive_control.pdf")
    fig.savefig(OUT / "fig3_positive_control.png")
    plt.close(fig)
    print("  fig3_positive_control done")


# ═══════════════════════════════════════════════════════════════════
# FIGURE 4 — LOD curve: DI vs sigma at multiple K values
# ═══════════════════════════════════════════════════════════════════
def fig4_lod_curve():
    with open(RESULTS / "15_positive_control" / "positive_control_validation.json") as f:
        pc = json.load(f)

    sb = pc["section_b_synthetic"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5.5, 2.6), gridspec_kw={"wspace": 0.35})

    # Left panel: mean DI vs sigma for each K
    for k_str in ["5", "10", "13", "20", "40"]:
        kdata = sb["lod_by_K"][k_str]
        sigmas = [pt["sigma"] for pt in kdata["sweep"]]
        di_means = [pt["di_mean"] for pt in kdata["sweep"]]
        k_int = int(k_str)
        ax1.plot(sigmas, di_means, color=K_COLORS[k_int], lw=1.3,
                 label=f"K = {k_str}", marker="o", markersize=2.5, zorder=3)

    ax1.set_xlabel("Angular spread σ (radians)")
    ax1.set_ylabel("Mean DI")
    ax1.legend(fontsize=6.5, frameon=False, loc="upper left")

    # Right panel: detection power (fraction above threshold) vs sigma
    for k_str in ["5", "10", "13", "20", "40"]:
        kdata = sb["lod_by_K"][k_str]
        sigmas = [pt["sigma"] for pt in kdata["sweep"]]
        fracs = [pt["frac_above_threshold"] for pt in kdata["sweep"]]
        k_int = int(k_str)
        ax2.plot(sigmas, fracs, color=K_COLORS[k_int], lw=1.3,
                 label=f"K = {k_str}", marker="o", markersize=2.5, zorder=3)

        # Mark LOD point
        lod_sigma = kdata["lod_sigma"]
        idx = next(i for i, s in enumerate(sigmas) if abs(s - lod_sigma) < 0.001)
        if idx < len(fracs):
            ax2.plot(lod_sigma, fracs[idx], "D", color=K_COLORS[k_int],
                     markersize=5, zorder=5, markeredgecolor="white",
                     markeredgewidth=0.6)

    ax2.axhline(0.95, color="0.5", ls="--", lw=0.5, zorder=1)
    ax2.text(1.35, 0.96, "95% power", fontsize=6, color="0.5", va="bottom", ha="right")
    ax2.set_xlabel("Angular spread σ (radians)")
    ax2.set_ylabel("Detection rate")
    ax2.set_ylim(-0.02, 1.05)
    ax2.legend(fontsize=6.5, frameon=False, loc="center right")

    # Panel labels
    ax1.text(-0.15, 1.08, "a", transform=ax1.transAxes, fontsize=11,
             fontweight="bold", va="top")
    ax2.text(-0.15, 1.08, "b", transform=ax2.transAxes, fontsize=11,
             fontweight="bold", va="top")

    fig.savefig(OUT / "fig4_lod_curve.pdf")
    fig.savefig(OUT / "fig4_lod_curve.png")
    plt.close(fig)
    print("  fig4_lod_curve done")


# ═══════════════════════════════════════════════════════════════════
# FIGURE 5 — HDAC selectivity gradient
# ═══════════════════════════════════════════════════════════════════
def fig5_hdac():
    """Read the committed source table; never hardcode plotted values."""
    src = json.loads((FIG5_SOURCE).read_text())
    rows = src["rows"]
    stats = src["statistics"]

    # ordered by selectivity category, then by plotted value within category
    rows = sorted(rows, key=lambda r: (r["category_rank"], r["corrected_instability"]))

    sel_colors = {
        "Pan-HDAC": "#5E3C99",
        "Class I-selective": "#B2ABD2",
        "Isoform-selective": "#E66101",
    }

    fig, ax = plt.subplots(figsize=(3.6, 2.6))

    y = np.arange(len(rows))[::-1]
    for i, r in enumerate(rows):
        di = r["corrected_instability"]
        ax.barh(y[i], di, height=0.55, color=sel_colors[r["selectivity_category"]],
                alpha=0.85, edgecolor="white", linewidth=0.5, zorder=3)
        ax.text(di + 0.015, y[i], f"{di:.2f} (n={r['n_celllines']})",
                va="center", ha="left", fontsize=7, fontweight="bold", color="0.25")

    ax.set_yticks(y)
    ax.set_yticklabels([r["drug"] for r in rows], fontsize=7.5)
    ax.set_xlabel("Direction instability (D), toxicity-corrected")
    ax.set_xlim(0, 1.12)

    handles = [mpatches.Patch(color=sel_colors[s], label=s, alpha=0.85)
               for s in ["Pan-HDAC", "Class I-selective", "Isoform-selective"]]
    ax.legend(handles=handles, fontsize=6, frameon=False, loc="upper right",
              title="Selectivity", title_fontsize=6.5)

    # every plotted value must equal its source row
    for i, r in enumerate(rows):
        plotted = ax.patches[i].get_width()
        assert abs(plotted - r["corrected_instability"]) < 1e-9, (
            f"Figure 5 plots {plotted} for {r['drug']}, source says "
            f"{r['corrected_instability']}")

    fig.savefig(OUT / "fig5_hdac.pdf")
    fig.savefig(OUT / "fig5_hdac.png")
    plt.close(fig)
    print(f"  fig5_hdac done (tau_b={stats['tau_b']}, "
          f"p={stats['exact_two_sided_p']}, values from source table)")


# ═══════════════════════════════════════════════════════════════════
# FIGURE 6 — Causal DAG: the confound structure (publication quality)
# ═══════════════════════════════════════════════════════════════════
def fig6_causal_dag():
    fig = plt.figure(figsize=(6.8, 2.7))

    # Three panels with controlled spacing
    gs = fig.add_gridspec(1, 3, wspace=0.45, left=0.04, right=0.98,
                          top=0.82, bottom=0.06)
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]

    # Node style
    NODE_RAD = 0.17
    CONFOUND_FC = "#FDDBC7"    # warm blush fill for confounded
    CONFOUND_EC = "#B2182B"    # brick edge
    NEUTRAL_FC = "#F0F0F0"
    NEUTRAL_EC = "#666666"
    DI_FC = "#D1E5F0"          # cool blue fill
    DI_EC = "#2166AC"          # strong blue edge

    def node(ax, x, y, text, fc, ec, fontsize=7, bold=False):
        box = mpatches.FancyBboxPatch(
            (x - 0.22, y - 0.11), 0.44, 0.22,
            boxstyle="round,pad=0.06", facecolor=fc,
            edgecolor=ec, linewidth=0.7, zorder=4, clip_on=False,
        )
        ax.add_patch(box)
        weight = "semibold" if bold else "normal"
        ax.text(x, y, text, ha="center", va="center", fontsize=fontsize,
                fontweight=weight, color=ec if ec != NEUTRAL_EC else "#333",
                zorder=5, clip_on=False)

    def edge(ax, x1, y1, x2, y2, color="0.45", lw=0.8, head="->",
             rad=0.0, ls="-"):
        ax.annotate(
            "", xy=(x2, y2), xytext=(x1, y1),
            arrowprops=dict(
                arrowstyle=head, color=color, lw=lw,
                connectionstyle=f"arc3,rad={rad}",
                shrinkA=0, shrinkB=0, ls=ls,
            ), zorder=3, clip_on=False,
        )

    def blocked_edge(ax, x1, y1, x2, y2, color=DI_EC):
        """Dashed arrow with a perpendicular barrier bar at the midpoint."""
        edge(ax, x1, y1, x2, y2, color=color, lw=0.7, ls="--")
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        dx, dy = x2 - x1, y2 - y1
        length = (dx**2 + dy**2) ** 0.5
        # perpendicular direction
        px, py = -dy / length, dx / length
        bar = 0.055
        ax.plot([mx - px*bar, mx + px*bar],
                [my - py*bar, my + py*bar],
                color=color, lw=2.2, solid_capstyle="round", zorder=6,
                clip_on=False)

    # Layout coordinates (shared across panels)
    TOP = 0.88
    MID = 0.44
    BOT = 0.0
    LEFT = 0.15
    RIGHT = 0.85
    CEN = 0.50

    for ax in axes:
        ax.set_xlim(-0.08, 1.08)
        ax.set_ylim(-0.18, 1.08)
        ax.axis("off")

    # ─── Panel a: DepMap ───────────────────────────────────────────
    ax = axes[0]
    node(ax, CEN, TOP, "Magnitude", NEUTRAL_FC, "#333333", bold=True)
    node(ax, LEFT, MID, "Essentiality\nlabel", CONFOUND_FC, CONFOUND_EC)
    node(ax, RIGHT, MID, "Score\nvariance", CONFOUND_FC, CONFOUND_EC)
    node(ax, CEN, BOT, "Observed\n'divergence'", CONFOUND_FC, CONFOUND_EC)

    edge(ax, CEN - 0.08, TOP - 0.12, LEFT + 0.10, MID + 0.12,
         color=CONFOUND_EC, lw=1.0)
    edge(ax, CEN + 0.08, TOP - 0.12, RIGHT - 0.10, MID + 0.12,
         color=CONFOUND_EC, lw=1.0)
    edge(ax, LEFT + 0.10, MID - 0.12, CEN - 0.08, BOT + 0.12,
         color="#888888", lw=0.7)
    edge(ax, RIGHT - 0.10, MID - 0.12, CEN + 0.08, BOT + 0.12,
         color=CONFOUND_EC, lw=1.0)

    ax.text(CEN, 1.06, "DepMap (variance)", ha="center", fontsize=8,
            fontweight="semibold", color="#333", clip_on=False)

    # ─── Panel b: Perturb-seq ──────────────────────────────────────
    ax = axes[1]
    node(ax, CEN, TOP, "Magnitude", NEUTRAL_FC, "#333333", bold=True)
    node(ax, LEFT, MID, "Essentiality\nstatus", CONFOUND_FC, CONFOUND_EC)
    node(ax, RIGHT, MID, "Estimation\nquality", CONFOUND_FC, CONFOUND_EC)
    node(ax, CEN, BOT, "Geodesic\ndistance", CONFOUND_FC, CONFOUND_EC)

    edge(ax, CEN - 0.08, TOP - 0.12, LEFT + 0.10, MID + 0.12,
         color=CONFOUND_EC, lw=1.0)
    edge(ax, CEN + 0.08, TOP - 0.12, RIGHT - 0.10, MID + 0.12,
         color=CONFOUND_EC, lw=1.0)
    edge(ax, LEFT + 0.10, MID - 0.12, CEN - 0.08, BOT + 0.12,
         color="#888888", lw=0.7)
    edge(ax, RIGHT - 0.10, MID - 0.12, CEN + 0.08, BOT + 0.12,
         color=CONFOUND_EC, lw=1.0)

    ax.text(CEN, 1.06, "Perturb-seq (estimation)", ha="center", fontsize=8,
            fontweight="semibold", color="#333", clip_on=False)

    # ─── Panel c: Direction instability (blocked) ──────────────────
    ax = axes[2]
    node(ax, CEN, TOP, "Magnitude", NEUTRAL_FC, "#333333", bold=True)
    node(ax, LEFT, MID, "Group\nlabel", NEUTRAL_FC, "#999999")
    node(ax, RIGHT, MID, "Cosine\nnormalization", DI_FC, DI_EC)
    node(ax, CEN, BOT, "Direction\ninstability", DI_FC, DI_EC)

    edge(ax, CEN - 0.08, TOP - 0.12, LEFT + 0.10, MID + 0.12,
         color="#999999", lw=0.7)
    blocked_edge(ax, CEN + 0.08, TOP - 0.12, RIGHT - 0.10, MID + 0.12,
                 color=DI_EC)
    edge(ax, LEFT + 0.10, MID - 0.12, CEN - 0.08, BOT + 0.12,
         color="#999999", lw=0.7)
    edge(ax, RIGHT - 0.10, MID - 0.12, CEN + 0.08, BOT + 0.12,
         color=DI_EC, lw=1.0)

    ax.text(CEN, 1.06, "Direction instability (blocked)", ha="center",
            fontsize=8, fontweight="semibold", color="#333", clip_on=False)

    fig.savefig(OUT / "fig6_causal_dag.pdf")
    fig.savefig(OUT / "fig6_causal_dag.png")
    plt.close(fig)
    print("  fig6_causal_dag done")


# ═══════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("Generating figures...")
    fig1_schematic()
    fig2_benchmark()
    fig3_positive_control()
    fig4_lod_curve()
    fig5_hdac()
    fig6_causal_dag()
    print(f"\nAll figures saved to {OUT}/")
