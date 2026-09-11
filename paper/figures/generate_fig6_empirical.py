"""Cross-domain comparison figure.

The JUMP-CP bars are read from the deposited Experiment 8 artifacts rather than
typed in, and asserted against them at build time: a figure that disagrees with
the table it illustrates is the failure this guards against.

The LINCS and Perturb-seq bars are still literals inherited from earlier drafts.
They are flagged here because they have no artifact behind them in this
repository; they should be bound to their source files before submission.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.size': 11,
    'font.family': 'serif',
    'mathtext.fontset': 'cm',
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
})

RESULTS = Path(__file__).resolve().parents[2] / "results" / "08_jump_cp_feature_ablation"
_j = json.loads((RESULTS / "feature_ablation_results.json").read_text())
_rows = (RESULTS / "compound_instability.csv").read_text().splitlines()[1:]
JUMP_N = int(_j["n_compounds"])
JUMP_RHO = float(_j["spearman_display"])
assert len(_rows) == JUMP_N, f"figure source disagrees with the table: {len(_rows)} vs {JUMP_N}"

# max absolute percentile-rank shift, converted to positions, from that same table
_raw = sorted(range(JUMP_N), key=lambda i: float(_rows[i].split(",")[2]))
_abl = sorted(range(JUMP_N), key=lambda i: float(_rows[i].split(",")[3]))
_rank_raw = {v: k for k, v in enumerate(_raw)}
_rank_abl = {v: k for k, v in enumerate(_abl)}
JUMP_MAX_SHIFT = max(abs(_rank_raw[i] - _rank_abl[i]) for i in range(JUMP_N))

# LINCS: recomputed from the deposited per-drug table, not typed in.
_lincs = json.loads((Path(__file__).resolve().parents[2] / "results" /
                     "01_toxicity_failure" / "toxicity_results.json").read_text())
_lr = sorted(range(len(_lincs)), key=lambda i: _lincs[i]["raw_instability"])
_lc = sorted(range(len(_lincs)), key=lambda i: _lincs[i]["corrected_instability"])
_pr = {v: k for k, v in enumerate(_lr)}
_pc = {v: k for k, v in enumerate(_lc)}
LINCS_MAX = max(abs(_pr[i] - _pc[i]) for i in range(len(_lincs)))

import numpy as _np
from scipy import stats as _st

LINCS_RHO = round(float(_st.spearmanr(
    [r["raw_instability"] for r in _lincs],
    [r["corrected_instability"] for r in _lincs]).statistic), 4)

# Perturb-seq: recomputed from the recovered per-gene distances.
_z = _np.load(Path(__file__).resolve().parents[2] / "results" /
              "07_perturbseq_correction" / "corrected_distances.npz")
_raw_d, _corr_d = _z["raw_dists"], _z["corr_dists"]
PSEQ_RHO = round(float(_st.spearmanr(_raw_d, _corr_d).statistic), 4)
PSEQ_MAX = int(_np.abs(_st.rankdata(_raw_d) - _st.rankdata(_corr_d)).max())

domains_left = ['LINCS\n(drugs)\nn = 8,949',
                'Perturb-seq\n(genes)\nn = 1,676',
                f'JUMP-CP\n(compounds)\nn = {JUMP_N:,}']

domains_right = ['LINCS\ntoxicity genes',
                 'Perturb-seq\nessential-gene\nsubspace',
                 'JUMP-CP\nmorphological\nblock']

rho_values = [LINCS_RHO, PSEQ_RHO, JUMP_RHO]
max_rank_shifts = [LINCS_MAX, PSEQ_MAX, JUMP_MAX_SHIFT]

colors = ['#2166ac', '#4393c3', '#92c5de']
edge_colors = ['#08519c', '#2171b5', '#4292c6']

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 3.6), gridspec_kw={'width_ratios': [1, 1]})

# Left panel: Spearman rho
bars1 = ax1.bar(range(3), rho_values, color=colors, edgecolor=edge_colors,
                linewidth=1.2, width=0.55, zorder=3)

for i, v in enumerate(rho_values):
    ax1.text(i, v + 0.006, f'{v:.3f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

ax1.set_xticks(range(3))
ax1.set_xticklabels(domains_left, fontsize=8, linespacing=1.15)
ax1.set_ylabel(r'Spearman $\rho$ (raw vs. corrected)', fontsize=10)
ax1.set_ylim(0.55, 1.05)
ax1.set_title('(a)  Population ranking preserved', fontsize=11, fontweight='bold', loc='left')
ax1.grid(axis='y', alpha=0.3, zorder=0)
ax1.set_axisbelow(True)

# Right panel: Max individual rank shift
bars2 = ax2.bar(range(3), max_rank_shifts, color=colors, edgecolor=edge_colors,
                linewidth=1.2, width=0.55, zorder=3)

for i, v in enumerate(max_rank_shifts):
    ax2.text(i, v + 180, f'{v:,}', ha='center', va='bottom', fontsize=10, fontweight='bold')

ax2.set_xticks(range(3))
ax2.set_xticklabels(domains_right, fontsize=8, linespacing=1.15)
ax2.set_ylabel('Max individual rank shift', fontsize=10)
ax2.set_title('(b)  Item-level confound identification', fontsize=11, fontweight='bold', loc='left')
ax2.grid(axis='y', alpha=0.3, zorder=0)
ax2.set_axisbelow(True)
ax2.set_ylim(0, max(max_rank_shifts) * 1.18)

plt.tight_layout(w_pad=2.5)

out_pdf = 'fig6_cross_domain_empirical.pdf'
out_png = 'fig6_cross_domain_empirical.png'
plt.savefig(out_pdf, bbox_inches='tight', dpi=300)
plt.savefig(out_png, bbox_inches='tight', dpi=300)
print(f'Saved {out_pdf} and {out_png}')
