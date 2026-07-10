#!/usr/bin/env python3
"""Figure: does the single-chain coordination advantage generalize from gp16 to ClpX?
Compares the design/native ratio of coordination proxies. Ratio > 1 = single-chain more
coordinated. gp16 (cp233) shows a robust 1.79x; ClpX (this study) sits at parity."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

prs = json.load(open('gp16_design/outputs/clpx_coupling/clpx_prs_anm_results.json'))
fn = json.load(open('gp16_design/outputs/clpx_coupling/clpx_force_network_results.json'))
d, n = prs['design'], prs['native']

INK = '#1b2a4a'; RED = '#c0392b'; BLUE = '#2e6db4'; GREY = '#9aa4b2'
plt.rcParams.update({'font.size': 11, 'axes.edgecolor': '#444', 'svg.fonttype': 'none'})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.2, 5.3))

# -------- Panel A: headline ratio, gp16 vs ClpX --------
labels = ['gp16 cp233\n(circular permutant,\napo vs apo, MD-validated)',
          'ClpX single-chain\n(direct fusion,\napo vs 6PP5 spiral, ENM only)']
ratios = [1.788, 1.084]
errs = [[0.0], [0.17]]  # ClpX cutoff spread 0.91-1.08 -> ~+/-0.17 around 1.0
colors = [RED, BLUE]
bars = ax1.bar(labels, ratios, color=colors, width=0.6, zorder=3, edgecolor='white')
# ClpX cutoff-sensitivity whisker (0.91 - 1.08)
ax1.plot([1, 1], [0.909, 1.084], color=INK, lw=2.5, zorder=4)
ax1.plot([0.93, 1.07], [0.909, 0.909], color=INK, lw=2.5, zorder=4)
ax1.plot([0.93, 1.07], [1.084, 1.084], color=INK, lw=2.5, zorder=4)
ax1.axhline(1.0, color=GREY, ls='--', lw=1.5, zorder=2)
ax1.text(1.46, 1.02, 'parity (no advantage)', color=GREY, fontsize=9, va='bottom', ha='right')
ax1.text(0, 1.788 + 0.03, '1.79x', ha='center', va='bottom', fontweight='bold', color=RED, fontsize=13)
ax1.text(1, 1.084 + 0.20, '~1.0x\n(0.91-1.08,\ncutoff-sensitive)', ha='center', va='bottom',
         fontweight='bold', color=BLUE, fontsize=11)
ax1.set_ylabel('PRS ATP -> neighbour-ATP coupling\n(single-chain / native)')
ax1.set_ylim(0, 2.1)
ax1.set_title('A. Primary coordination proxy: gp16 vs ClpX', fontweight='bold', loc='left')
for s in ('top', 'right'):
    ax1.spines[s].set_visible(False)

# -------- Panel B: ClpX proxies, design vs native --------
metrics = ['PRS ATP->NN\ncoupling', 'PRS inter/intra\n(whole ring)',
           'GNM ATP-NN\nco-fluct.', 'R307 arg-finger\nbetweenness']
des_vals = [d['prs_ATP_NN_coupling'], d['prs_inter_over_intra_all'],
            d['gnm_ATP_NN_crosscorr'], fn['A_design']['R307_argfinger_betweenness']]
nat_vals = [n['prs_ATP_NN_coupling'], n['prs_inter_over_intra_all'],
            n['gnm_ATP_NN_crosscorr'], fn['A_native']['R307_argfinger_betweenness']]
# normalize each metric to native = 1.0 for visual comparison
des_norm = [dv / nv for dv, nv in zip(des_vals, nat_vals)]
nat_norm = [1.0] * len(metrics)
x = np.arange(len(metrics)); w = 0.38
ax2.bar(x - w/2, nat_norm, w, label='native ClpX (6PP5)', color=GREY, zorder=3, edgecolor='white')
ax2.bar(x + w/2, des_norm, w, label='single-chain design', color=BLUE, zorder=3, edgecolor='white')
ax2.axhline(1.0, color='#555', lw=1.0, zorder=2)
for i, v in enumerate(des_norm):
    ax2.text(x[i] + w/2, v + 0.02, f'{v:.2f}x', ha='center', va='bottom', fontsize=9,
             color=(BLUE if v >= 1 else RED), fontweight='bold')
ax2.set_xticks(x); ax2.set_xticklabels(metrics, fontsize=9)
ax2.set_ylabel('design / native  (native = 1.0)')
ax2.set_ylim(0, 1.45)
ax2.legend(frameon=False, fontsize=9, loc='upper right')
ax2.set_title('B. ClpX proxies: mixed, near parity', fontweight='bold', loc='left')
for s in ('top', 'right'):
    ax2.spines[s].set_visible(False)

fig.suptitle('Does "single-chain > native coordination" generalize beyond gp16?  '
             '-> For ClpX, no: it is at parity, not 1.8x.',
             fontsize=12.5, fontweight='bold', y=1.005)
fig.tight_layout()
out = 'gp16_design/outputs/clpx_coupling/clpx_vs_gp16_coordination.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
print('wrote', out)
