#!/usr/bin/env python3
"""Analyze the integrated helical<->planar + descent-grip run: concerted vs sequential, averaged over seeds.
Reads integ_{concerted,sequential}_s*/series.json. Reports DNA net (mean±sd over seeds), coupling, and a
comparison figure (per-subunit z(t) for one seed of each mode + DNA displacement across seeds).
"""
import argparse, os, glob, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SUB_COLORS = ['#e6194B', '#3cb44b', '#4363d8', '#f58231', '#911eb4']


def load_mode(d, mode):
    runs = []
    for p in sorted(glob.glob(os.path.join(d, f'integ_{mode}_s*'))):
        j = os.path.join(p, 'series.json')
        if os.path.exists(j):
            js = json.load(open(j))
            rows = js['rows']
            runs.append(dict(t=np.array([r['t_cyc'] for r in rows]),
                             z=np.array([r['per_sub_z_A'] for r in rows]),
                             dna=np.array([r['dna_z_A'] for r in rows]) - rows[0]['dna_z_A'],
                             neng=np.array([r['n_engaged'] for r in rows]),
                             net=js['result']['dna_net_A'], nmin=js['result']['n_eng_min']))
    return runs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default='../../outputs/php_cycle/integrated')
    ap.add_argument('--out', default='../../outputs/php_cycle/integrated')
    args = ap.parse_args()
    modes = ['concerted', 'sequential']
    data = {m: load_mode(args.dir, m) for m in modes}
    summ = {}
    for m in modes:
        nets = [r['net'] for r in data[m]]
        if nets:
            summ[m] = dict(n_seeds=len(nets), dna_net_mean=round(float(np.mean(nets)), 2),
                           dna_net_sd=round(float(np.std(nets)), 2), dna_net_vals=nets,
                           n_eng_min=int(min(r['nmin'] for r in data[m])))
    json.dump(summ, open(os.path.join(args.out, 'integrated_summary.json'), 'w'), indent=2)

    fig, axes = plt.subplots(2, 2, figsize=(12, 7))
    for ci, m in enumerate(modes):
        if not data[m]:
            axes[0, ci].text(0.5, 0.5, f'{m}: no data', ha='center'); continue
        r0 = data[m][0]
        for k in range(5): axes[0, ci].plot(r0['t'], r0['z'][:, k], color=SUB_COLORS[k], lw=1.5, label=f'sub{k+1}')
        axes[0, ci].set_title(f'{m.upper()} — per-subunit axial z(t) (seed 1)\nhelical↔planar staircase cycle', fontsize=9)
        axes[0, ci].set_ylabel('subunit z (Å)'); axes[0, ci].grid(alpha=0.2)
        if ci == 0: axes[0, ci].legend(fontsize=7, ncol=5, loc='upper center')
        for r in data[m]: axes[1, ci].plot(r['t'], r['dna'], lw=1.5, alpha=0.7)
        s = summ.get(m, {})
        axes[1, ci].set_title(f'DNA axial displacement — {len(data[m])} seeds\nnet {s.get("dna_net_mean","?")}±{s.get("dna_net_sd","?")} Å', fontsize=9)
        axes[1, ci].set_ylabel('DNA z (Å)'); axes[1, ci].set_xlabel('cycles'); axes[1, ci].grid(alpha=0.2)
    fig.suptitle('INTEGRATED helical↔planar + descent-grip (real DNA-contact residues): concerted vs sequential', fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(args.out, 'integrated_compare.png'), dpi=150); plt.close(fig)
    print('wrote integrated_compare.png ; summary:', json.dumps(summ))


if __name__ == '__main__':
    main()
