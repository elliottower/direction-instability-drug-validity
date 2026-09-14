"""Regenerate the 4 figures used in the BMC MRM v2 draft.

Fig1/Fig6: Nature-style DAGs — identical node structure in both panels,
color shows which paths are active vs blocked.
Fig5/FigH5: NeurIPS-quality data figures.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(exist_ok=True)

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
    "text.usetex": False,
    "mathtext.fontset": "dejavusans",
})

C_DI = "#2166AC"
C_FRECH = "#B2182B"
C_BIO = "#2E7D32"
C_DARK = "#2C3E50"
C_MUTED = "#BBBBBB"


# ── DAG drawing helpers (Nature-style rounded rectangles) ─────────

def _rect_node(ax, x, y, text, w=1.8, h=0.7, fc="white", ec=C_DARK,
               lw=1.0, fontsize=10, fontcolor=None, bold=False, ls="-"):
    rect = mpatches.FancyBboxPatch(
        (x - w/2, y - h/2), w, h,
        boxstyle="round,pad=0.15", facecolor=fc,
        edgecolor=ec, linewidth=lw, linestyle=ls, zorder=10, clip_on=False)
    ax.add_patch(rect)
    weight = "bold" if bold else "normal"
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize,
            fontweight=weight, color=fontcolor or ec, zorder=11,
            clip_on=False)


def _arrow(ax, x1, y1, x2, y2, color=C_DARK, lw=1.2,
           shrinkA=14, shrinkB=14):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="->,head_length=0.35,head_width=0.2",
                                color=color, lw=lw,
                                shrinkA=shrinkA, shrinkB=shrinkB),
                zorder=3, clip_on=False)


def _gapped_arrow(ax, x1, y1, x2, y2, color=C_MUTED, lw=0.8,
                   shrinkA=12, shrinkB=12):
    """Arrow with solid head but two gaps in the shaft.

    Uses annotate with a custom dash pattern: solid at both ends
    with two gaps in the middle third and two-thirds of the shaft.
    """
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>,head_length=0.35,head_width=0.2",
                                color=color, lw=lw,
                                shrinkA=shrinkA, shrinkB=shrinkB,
                                ls=(0, (6, 2, 6, 2))),
                zorder=3, clip_on=False)


def _blocked_arrow(ax, x1, y1, x2, y2, color=C_FRECH, lw=1.2):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="->,head_length=0.35,head_width=0.2",
                                color=color, lw=lw,
                                shrinkA=14, shrinkB=14, ls=(0, (4, 3))),
                zorder=3, clip_on=False)
    mx, my = (x1 + x2) / 2, (y1 + y2) / 2
    size = 0.22
    ax.plot([mx - size, mx + size], [my - size, my + size],
            color=color, lw=2.5, solid_capstyle="round", zorder=12,
            clip_on=False)
    ax.plot([mx - size, mx + size], [my + size, my - size],
            color=color, lw=2.5, solid_capstyle="round", zorder=12,
            clip_on=False)


# ═══════════════════════════════════════════════════════════════════
# FIGURE 1 — DAG: magnitude confound vs direction instability
#
# Same 4 nodes in both panels. Difference: which paths are active.
#
#   Perturbation
#    /         \
# Magnitude   Biological signal
#    \         /
#     Score
#
# Panel a: both paths active (confounded)
# Panel b: magnitude path blocked by cosine normalization
# ═══════════════════════════════════════════════════════════════════
def fig1_dag_v2():
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(7.5, 4.5),
                                      gridspec_kw={"wspace": 0.40})

    PERT = (3.0, 6.0)
    MAG  = (1.0, 3.8)
    BIO  = (5.0, 3.8)
    SCORE = (3.0, 1.6)

    NW, NH = 2.1, 0.7

    for ax in (ax_a, ax_b):
        ax.set_xlim(-0.8, 6.8)
        ax.set_ylim(0.3, 7.5)
        ax.set_aspect("equal")
        ax.axis("off")

    # ── Panel a: Confounded ──────────────────────────────────────
    _arrow(ax_a, PERT[0]-0.4, PERT[1]-0.38, MAG[0]+0.4, MAG[1]+0.38,
           color=C_FRECH, lw=1.4)
    _arrow(ax_a, PERT[0]+0.4, PERT[1]-0.38, BIO[0]-0.4, BIO[1]+0.38,
           color=C_BIO, lw=1.4)
    _arrow(ax_a, MAG[0]+0.4, MAG[1]-0.38, SCORE[0]-0.4, SCORE[1]+0.38,
           color=C_FRECH, lw=1.4)
    _arrow(ax_a, BIO[0]-0.4, BIO[1]-0.38, SCORE[0]+0.4, SCORE[1]+0.38,
           color=C_BIO, lw=1.4)

    _rect_node(ax_a, *PERT, "Perturbation", w=NW, h=NH,
               fc="#F0F0F0", ec=C_DARK, fontcolor="#333", bold=True, fontsize=11)
    _rect_node(ax_a, *MAG, "Magnitude", w=NW, h=NH,
               fc="#FDDBC7", ec=C_FRECH, fontcolor=C_FRECH, fontsize=10)
    _rect_node(ax_a, *BIO, "Biological\nsignal", w=NW, h=NH+0.2,
               fc="#E8F5E9", ec=C_BIO, fontcolor=C_BIO, fontsize=10)
    _rect_node(ax_a, *SCORE, "Metric score", w=NW+0.2, h=NH,
               fc="#FDDBC7", ec=C_FRECH, fontcolor=C_FRECH, fontsize=10)

    ax_a.text(SCORE[0], SCORE[1]-0.6, "magnitude + biology entangled",
              ha="center", fontsize=8, color="#666", style="italic")

    ax_a.text(3.0, 7.35, "Confounded metrics", ha="center", fontsize=13,
              fontweight="bold", color="#333")
    ax_a.text(-0.6, 7.35, "a", fontsize=15, fontweight="bold", color="#333")

    # ── Panel b: Direction instability (blocked) ─────────────────
    _arrow(ax_b, PERT[0]-0.4, PERT[1]-0.38, MAG[0]+0.4, MAG[1]+0.38,
           color=C_MUTED, lw=0.8)
    _arrow(ax_b, PERT[0]+0.4, PERT[1]-0.38, BIO[0]-0.4, BIO[1]+0.38,
           color=C_BIO, lw=1.4)
    _blocked_arrow(ax_b, MAG[0]+0.4, MAG[1]-0.38,
                   SCORE[0]-0.4, SCORE[1]+0.38, color=C_MUTED)
    _arrow(ax_b, BIO[0]-0.4, BIO[1]-0.38, SCORE[0]+0.4, SCORE[1]+0.38,
           color=C_BIO, lw=1.4)

    _rect_node(ax_b, *PERT, "Perturbation", w=NW, h=NH,
               fc="#F0F0F0", ec=C_DARK, fontcolor="#333", bold=True, fontsize=11)
    _rect_node(ax_b, *MAG, "Magnitude", w=NW, h=NH,
               fc="#F8F8F8", ec=C_MUTED, fontcolor=C_MUTED, fontsize=10)
    _rect_node(ax_b, *BIO, "Biological\nsignal", w=NW, h=NH+0.2,
               fc="#E8F5E9", ec=C_BIO, fontcolor=C_BIO, fontsize=10)
    _rect_node(ax_b, *SCORE, "Direction\ninstability", w=NW+0.2, h=NH+0.2,
               fc="#D1E5F0", ec=C_DI, fontcolor=C_DI, fontsize=10)

    ax_b.text(SCORE[0], SCORE[1]-0.7, "biology only",
              ha="center", fontsize=8, color=C_BIO, style="italic")

    ax_b.text(3.0, 7.35, "Direction instability", ha="center", fontsize=13,
              fontweight="bold", color="#333")
    ax_b.text(-0.6, 7.35, "b", fontsize=15, fontweight="bold", color="#333")

    fig.savefig(OUT / "fig1_dag_v2.pdf")
    fig.savefig(OUT / "fig1_dag_v2.png")
    plt.close(fig)
    print(f"  fig1_dag_v2 done")


def fig1_dag_v3():
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(6.5, 3.2),
                                      gridspec_kw={"wspace": 0.35})

    NW, NH = 1.7, 0.55
    FS = 9

    PERT  = (2.5, 4.6)
    MAG   = (0.8, 2.8)
    BIO   = (4.2, 2.8)
    SCORE = (2.5, 1.0)

    for ax in (ax_a, ax_b):
        ax.set_xlim(-0.5, 5.5)
        ax.set_ylim(0.1, 5.8)
        ax.set_aspect("equal")
        ax.axis("off")

    # ── Panel a: Confounded ──────────────────────────────────────
    _arrow(ax_a, PERT[0]-0.3, PERT[1]-0.30, MAG[0]+0.3, MAG[1]+0.30,
           color=C_FRECH, lw=1.2, shrinkA=12, shrinkB=12)
    _arrow(ax_a, PERT[0]+0.3, PERT[1]-0.30, BIO[0]-0.3, BIO[1]+0.30,
           color=C_BIO, lw=1.2, shrinkA=12, shrinkB=12)
    _arrow(ax_a, MAG[0]+0.3, MAG[1]-0.30, SCORE[0]-0.3, SCORE[1]+0.30,
           color=C_FRECH, lw=1.2, shrinkA=12, shrinkB=12)
    _arrow(ax_a, BIO[0]-0.3, BIO[1]-0.30, SCORE[0]+0.3, SCORE[1]+0.30,
           color=C_BIO, lw=1.2, shrinkA=12, shrinkB=12)

    _rect_node(ax_a, *PERT, "Perturbation", w=NW, h=NH,
               fc="#F0F0F0", ec="#888", fontcolor="#444", fontsize=FS)
    _rect_node(ax_a, *MAG, "Magnitude", w=NW, h=NH,
               fc="#FDDBC7", ec=C_FRECH, fontcolor=C_FRECH, fontsize=FS)
    _rect_node(ax_a, *BIO, "Biological\nsignal", w=NW, h=NH+0.15,
               fc="#E8F5E9", ec=C_BIO, fontcolor=C_BIO, fontsize=FS)
    _rect_node(ax_a, *SCORE, "Metric score", w=NW+0.1, h=NH,
               fc="#FDDBC7", ec=C_FRECH, fontcolor=C_FRECH, fontsize=FS)

    ax_a.text(2.5, 5.55, "Confounded metrics", ha="center", fontsize=11,
              fontweight="bold", color="#333")
    ax_a.text(-0.3, 5.55, "a", fontsize=13, fontweight="bold", color="#333")

    # ── Panel b: Direction instability (dashed magnitude box) ─────
    _arrow(ax_b, PERT[0]-0.3, PERT[1]-0.30, MAG[0]+0.3, MAG[1]+0.30,
           color=C_MUTED, lw=0.8, shrinkA=12, shrinkB=12)
    _arrow(ax_b, PERT[0]+0.3, PERT[1]-0.30, BIO[0]-0.3, BIO[1]+0.30,
           color=C_BIO, lw=1.2, shrinkA=12, shrinkB=12)
    _arrow(ax_b, MAG[0]+0.3, MAG[1]-0.30, SCORE[0]-0.3, SCORE[1]+0.30,
           color=C_MUTED, lw=0.8, shrinkA=12, shrinkB=12)
    _arrow(ax_b, BIO[0]-0.3, BIO[1]-0.30, SCORE[0]+0.3, SCORE[1]+0.30,
           color=C_BIO, lw=1.2, shrinkA=12, shrinkB=12)

    _rect_node(ax_b, *PERT, "Perturbation", w=NW, h=NH,
               fc="#F0F0F0", ec="#888", fontcolor="#444", fontsize=FS)
    _rect_node(ax_b, *MAG, "Magnitude", w=NW, h=NH,
               fc="#EDEDED", ec="#999", fontcolor="#888", fontsize=FS,
               ls=(0, (4, 3)))
    _rect_node(ax_b, *BIO, "Biological\nsignal", w=NW, h=NH+0.15,
               fc="#E8F5E9", ec=C_BIO, fontcolor=C_BIO, fontsize=FS)
    _rect_node(ax_b, *SCORE, "Direction\ninstability", w=NW+0.1, h=NH+0.15,
               fc="#D1E5F0", ec=C_DI, fontcolor=C_DI, fontsize=FS)

    # Centered under Magnitude, halfway to Score
    ax_b.text(MAG[0], (MAG[1] + SCORE[1]) / 2, "cosine\nnormalization",
              ha="center", va="center", fontsize=6.5, color="#999",
              style="italic")

    ax_b.text(2.5, 5.55, "Direction instability", ha="center", fontsize=11,
              fontweight="bold", color="#333")
    ax_b.text(-0.3, 5.55, "b", fontsize=13, fontweight="bold", color="#333")

    fig.savefig(OUT / "fig1_dag_v3.pdf")
    fig.savefig(OUT / "fig1_dag_v3.png")
    plt.close(fig)
    print(f"  fig1_dag_v3 done")


def fig1_dag_v4():
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(6.5, 3.2),
                                      gridspec_kw={"wspace": 0.35})

    NW, NH = 1.7, 0.55
    FS = 9

    PERT  = (2.5, 4.6)
    MAG   = (0.8, 2.8)
    BIO   = (4.2, 2.8)
    SCORE = (2.5, 1.0)

    for ax in (ax_a, ax_b):
        ax.set_xlim(-0.5, 5.5)
        ax.set_ylim(0.1, 5.8)
        ax.set_aspect("equal")
        ax.axis("off")

    # ── Panel a: Confounded ──────────────────────────────────────
    _arrow(ax_a, PERT[0]-0.3, PERT[1]-0.30, MAG[0]+0.3, MAG[1]+0.30,
           color=C_FRECH, lw=1.2, shrinkA=12, shrinkB=12)
    _arrow(ax_a, PERT[0]+0.3, PERT[1]-0.30, BIO[0]-0.3, BIO[1]+0.30,
           color=C_BIO, lw=1.2, shrinkA=12, shrinkB=12)
    _arrow(ax_a, MAG[0]+0.3, MAG[1]-0.30, SCORE[0]-0.3, SCORE[1]+0.30,
           color=C_FRECH, lw=1.2, shrinkA=12, shrinkB=12)
    _arrow(ax_a, BIO[0]-0.3, BIO[1]-0.30, SCORE[0]+0.3, SCORE[1]+0.30,
           color=C_BIO, lw=1.2, shrinkA=12, shrinkB=12)

    _rect_node(ax_a, *PERT, "Perturbation", w=NW, h=NH,
               fc="#F0F0F0", ec="#888", fontcolor="#444", fontsize=FS)
    _rect_node(ax_a, *MAG, "Magnitude", w=NW, h=NH,
               fc="#FDDBC7", ec=C_FRECH, fontcolor=C_FRECH, fontsize=FS)
    _rect_node(ax_a, *BIO, "Biological\nsignal", w=NW, h=NH+0.15,
               fc="#E8F5E9", ec=C_BIO, fontcolor=C_BIO, fontsize=FS)
    _rect_node(ax_a, *SCORE, "Metric score", w=NW+0.1, h=NH,
               fc="#FDDBC7", ec=C_FRECH, fontcolor=C_FRECH, fontsize=FS)

    ax_a.text(2.5, 5.55, "Confounded metrics", ha="center", fontsize=11,
              fontweight="bold", color="#333")
    ax_a.text(-0.3, 5.55, "a", fontsize=13, fontweight="bold", color="#333")

    # ── Panel b: Direction instability (gapped magnitude arrows) ────
    _gapped_arrow(ax_b, PERT[0]-0.3, PERT[1]-0.30, MAG[0]+0.3, MAG[1]+0.30,
                  color=C_MUTED, lw=0.8, shrinkA=12, shrinkB=12)
    _arrow(ax_b, PERT[0]+0.3, PERT[1]-0.30, BIO[0]-0.3, BIO[1]+0.30,
           color=C_BIO, lw=1.2, shrinkA=12, shrinkB=12)
    _gapped_arrow(ax_b, MAG[0]+0.3, MAG[1]-0.30, SCORE[0]-0.3, SCORE[1]+0.30,
                  color=C_MUTED, lw=0.8, shrinkA=12, shrinkB=12)
    _arrow(ax_b, BIO[0]-0.3, BIO[1]-0.30, SCORE[0]+0.3, SCORE[1]+0.30,
           color=C_BIO, lw=1.2, shrinkA=12, shrinkB=12)

    _rect_node(ax_b, *PERT, "Perturbation", w=NW, h=NH,
               fc="#F0F0F0", ec="#888", fontcolor="#444", fontsize=FS)
    _rect_node(ax_b, *MAG, "Magnitude", w=NW, h=NH,
               fc="#F8F8F8", ec=C_MUTED, fontcolor=C_MUTED, fontsize=FS,
               ls=(0, (4, 3)))
    _rect_node(ax_b, *BIO, "Biological\nsignal", w=NW, h=NH+0.15,
               fc="#E8F5E9", ec=C_BIO, fontcolor=C_BIO, fontsize=FS)
    _rect_node(ax_b, *SCORE, "Direction\ninstability", w=NW+0.1, h=NH+0.15,
               fc="#D1E5F0", ec=C_DI, fontcolor=C_DI, fontsize=FS)

    ax_b.text(MAG[0], (MAG[1] + SCORE[1]) / 2, "cosine\nnormalization",
              ha="center", va="center", fontsize=6.5, color="#999",
              style="italic")

    ax_b.text(2.5, 5.55, "Direction instability", ha="center", fontsize=11,
              fontweight="bold", color="#333")
    ax_b.text(-0.3, 5.55, "b", fontsize=13, fontweight="bold", color="#333")

    fig.savefig(OUT / "fig1_dag_v4.pdf")
    fig.savefig(OUT / "fig1_dag_v4.png")
    plt.close(fig)
    print(f"  fig1_dag_v4 done")


# ═══════════════════════════════════════════════════════════════════
# FIGURE 6 — Cross-domain confound structures (3 panels)
#
# Same diamond layout as fig1, but the confound mechanism differs:
#   DepMap: magnitude → variance
#   Perturb-seq: magnitude → estimation quality
#   DI: magnitude blocked
# ═══════════════════════════════════════════════════════════════════
def fig6_causal_dag():
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.5),
                              gridspec_kw={"wspace": 0.30})

    PERT = (3.0, 6.0)
    MAG  = (1.0, 3.8)
    MECH = (5.0, 3.8)
    SCORE = (3.0, 1.6)

    NW, NH = 2.1, 0.65

    panels = [
        {
            "title": "DepMap",
            "label": "a",
            "mechanism": "Score\nvariance",
            "score": "Observed\ndivergence",
            "confounded": True,
        },
        {
            "title": "Perturb-seq",
            "label": "b",
            "mechanism": "Estimation\nquality",
            "score": "Geodesic\ndistance",
            "confounded": True,
        },
        {
            "title": "Direction instability",
            "label": "c",
            "mechanism": "Biological\nsignal",
            "score": "Direction\ninstability",
            "confounded": False,
        },
    ]

    for ax in axes:
        ax.set_xlim(-0.8, 6.8)
        ax.set_ylim(0.3, 7.5)
        ax.set_aspect("equal")
        ax.axis("off")

    for ax, p in zip(axes, panels):
        if p["confounded"]:
            _arrow(ax, PERT[0]-0.4, PERT[1]-0.35, MAG[0]+0.4, MAG[1]+0.35,
                   color=C_FRECH, lw=1.4)
            _arrow(ax, PERT[0]+0.4, PERT[1]-0.35, MECH[0]-0.4, MECH[1]+0.35,
                   color=C_FRECH, lw=1.4)
            _arrow(ax, MAG[0]+0.4, MAG[1]-0.35, SCORE[0]-0.4, SCORE[1]+0.35,
                   color=C_FRECH, lw=1.4)
            _arrow(ax, MECH[0]-0.4, MECH[1]-0.35, SCORE[0]+0.4, SCORE[1]+0.35,
                   color=C_FRECH, lw=1.4)

            _rect_node(ax, *PERT, "Perturbation", w=NW+0.2, h=NH,
                       fc="#F0F0F0", ec=C_DARK, fontcolor="#333",
                       bold=True, fontsize=10)
            _rect_node(ax, *MAG, "Magnitude", w=NW, h=NH,
                       fc="#FDDBC7", ec=C_FRECH, fontcolor=C_FRECH,
                       fontsize=9.5)
            _rect_node(ax, *MECH, p["mechanism"], w=NW, h=NH+0.2,
                       fc="#FDDBC7", ec=C_FRECH, fontcolor=C_FRECH,
                       fontsize=9.5)
            _rect_node(ax, *SCORE, p["score"], w=NW+0.2, h=NH+0.2,
                       fc="#FDDBC7", ec=C_FRECH, fontcolor=C_FRECH,
                       fontsize=9.5)
        else:
            _arrow(ax, PERT[0]-0.4, PERT[1]-0.35, MAG[0]+0.4, MAG[1]+0.35,
                   color=C_MUTED, lw=0.8)
            _arrow(ax, PERT[0]+0.4, PERT[1]-0.35, MECH[0]-0.4, MECH[1]+0.35,
                   color=C_BIO, lw=1.4)
            _blocked_arrow(ax, MAG[0]+0.4, MAG[1]-0.35,
                           SCORE[0]-0.4, SCORE[1]+0.35, color=C_MUTED)
            _arrow(ax, MECH[0]-0.4, MECH[1]-0.35, SCORE[0]+0.4, SCORE[1]+0.35,
                   color=C_BIO, lw=1.4)

            _rect_node(ax, *PERT, "Perturbation", w=NW+0.2, h=NH,
                       fc="#F0F0F0", ec=C_DARK, fontcolor="#333",
                       bold=True, fontsize=10)
            _rect_node(ax, *MAG, "Magnitude", w=NW, h=NH,
                       fc="#F8F8F8", ec=C_MUTED, fontcolor=C_MUTED,
                       fontsize=9.5)
            _rect_node(ax, *MECH, p["mechanism"], w=NW, h=NH+0.2,
                       fc="#E8F5E9", ec=C_BIO, fontcolor=C_BIO,
                       fontsize=9.5)
            _rect_node(ax, *SCORE, p["score"], w=NW+0.2, h=NH+0.2,
                       fc="#D1E5F0", ec=C_DI, fontcolor=C_DI,
                       fontsize=9.5)

        ax.text(3.0, 7.3, p["title"], ha="center", fontsize=11,
                fontweight="bold", color="#333")
        ax.text(-0.6, 7.3, p["label"], fontsize=13, fontweight="bold",
                color="#333")

    fig.savefig(OUT / "fig6_causal_dag.pdf")
    fig.savefig(OUT / "fig6_causal_dag.png")
    plt.close(fig)
    print(f"  fig6_causal_dag done")


# ═══════════════════════════════════════════════════════════════════
# FIGURE 5 — HDAC selectivity gradient
# ═══════════════════════════════════════════════════════════════════
def fig5_hdac():
    drugs = [
        ("Panobinostat", "Pan-HDAC", 0.14),
        ("Vorinostat", "Pan-HDAC", 0.16),
        ("Belinostat", "Pan-HDAC", 0.19),
        ("Entinostat", "Class I", 0.27),
        ("Tubacin", "HDAC6", 0.32),
        ("PCI-34051", "HDAC8", 0.45),
    ]

    sel_colors = {
        "Pan-HDAC": "#5E3C99",
        "Class I": "#B2ABD2",
        "HDAC6": "#E66101",
        "HDAC8": "#FDB863",
    }

    fig, ax = plt.subplots(figsize=(3.8, 2.4))

    y = np.arange(len(drugs))[::-1]
    for i, (name, sel, di) in enumerate(drugs):
        color = sel_colors[sel]
        ax.barh(y[i], di, height=0.6, color=color, alpha=0.90,
                edgecolor="none", zorder=3)
        ax.text(di + 0.01, y[i], f"{di:.2f}", va="center", ha="left",
                fontsize=7.5, fontweight="bold", color="0.3")

    ax.set_yticks(y)
    ax.set_yticklabels([d[0] for d in drugs], fontsize=7.5)
    ax.set_xlabel("Direction instability (D)")
    ax.set_xlim(0, 0.58)

    handles = [mpatches.Patch(color=sel_colors[s], label=s, alpha=0.90)
               for s in ["Pan-HDAC", "Class I", "HDAC6", "HDAC8"]]
    ax.legend(handles=handles, fontsize=6.5, frameon=True, facecolor="white",
              edgecolor="0.85", loc="center right",
              title="Selectivity", title_fontsize=7,
              bbox_to_anchor=(1.0, 0.55))

    fig.savefig(OUT / "fig5_hdac.pdf")
    fig.savefig(OUT / "fig5_hdac.png")
    plt.close(fig)
    print(f"  fig5_hdac done")


# ═══════════════════════════════════════════════════════════════════
# FIGURE H5 — Leave-one-out delta-rho
# ═══════════════════════════════════════════════════════════════════
def fig_h5_full_delta_rho():
    with open(RESULTS / "05c_h5_full" / "h5_full_leave_one_out_results.json") as f:
        folds = json.load(f)

    delta_rhos = [f["ts_spearman"] - f["raw_spearman"] for f in folds]
    n_drugs = [f["n_drugs"] for f in folds]
    mean_delta = np.mean(delta_rhos)
    freq_corr = spearmanr(n_drugs, delta_rhos)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5.5, 2.6),
                                    gridspec_kw={"wspace": 0.40})

    ax1.hist(delta_rhos, bins=15, color=C_DI, edgecolor="white",
             linewidth=0.5, alpha=0.85, zorder=3)
    ax1.axvline(0, color="0.6", ls="--", lw=0.5, zorder=1)
    ax1.axvline(mean_delta, color=C_FRECH, lw=1.2, zorder=4,
                label=f"Mean = {mean_delta:.2f}")
    ax1.set_xlabel(r"$\Delta\rho$ (TS $-$ raw)")
    ax1.set_ylabel("Number of folds")
    ax1.set_title(r"All 66 folds: $\Delta\rho > 0$", fontsize=9)
    ax1.legend(fontsize=7, frameon=True, facecolor="white", edgecolor="0.8")

    ax2.scatter(n_drugs, delta_rhos, s=18, color=C_DI, alpha=0.6,
                edgecolors="none", zorder=3)
    ax2.set_xlabel("Drugs per fold")
    ax2.set_xscale("log")
    ax2.set_ylabel(r"$\Delta\rho$")
    rho_str = f"{freq_corr.statistic:.2f}"
    p_str = f"{freq_corr.pvalue:.3f}"
    ax2.set_title(f"Frequency gradient\n"
                  rf"$\rho$ = {rho_str}, p = {p_str}", fontsize=9)

    fig.savefig(OUT / "fig_h5_full_delta_rho.pdf")
    fig.savefig(OUT / "fig_h5_full_delta_rho.png")
    plt.close(fig)
    print(f"  fig_h5_full_delta_rho done")


if __name__ == "__main__":
    print("Regenerating BMC MRM figures...")
    fig1_dag_v2()
    fig1_dag_v3()
    fig1_dag_v4()
    fig5_hdac()
    fig6_causal_dag()
    fig_h5_full_delta_rho()
    print(f"\nAll figures saved to {OUT}/")
