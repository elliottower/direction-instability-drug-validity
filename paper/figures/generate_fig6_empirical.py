"""Generate empirical cross-domain comparison figure (replaces conceptual DAG)."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    'font.size': 11,
    'font.family': 'serif',
    'mathtext.fontset': 'cm',
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
})

domains_left = ['LINCS\n(drugs)\nn = 8,949',
                'Perturb-seq\n(genes)\nn = 1,676',
                'JUMP-CP\n(compounds)\nn = 25,254']

domains_right = ['LINCS\ntoxicity genes',
                 'Perturb-seq\nessential-gene\nsubspace',
                 'JUMP-CP\ncell-health\nfeatures']

rho_values = [0.9991, 0.91, 0.987]
max_rank_shifts = [1281, 900, 7155]

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
