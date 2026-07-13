#!/usr/bin/env python3
"""Compare the CONCERTED vs SEQUENTIAL ratchet runs — the 3D-MINFLUX discriminator.
Reads ratchet_{concerted,sequential}/series.json and makes a 2-column figure:
  row1 per-subunit axial z(t)  (concerted = 5 curves move TOGETHER; sequential = STAGGERED)
  row2 DNA axial displacement  (concerted = big steps; sequential = small substeps)
  row3 M2 coupling n_engaged
So MINFLUX would distinguish the mechanisms by (a) whether the 5 subunit traces move in phase or in a
travelling wave, and (b) whether DNA advances in one big step or several substeps per cycle.
"""
import argparse, os, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SUB_COLORS = ['#e6194B', '#3cb44b', '#4363d8', '#f58231', '#911eb4']


def load(d):
    j = json.load(open(os.path.join(d, 'series.json')))
    rows = j['rows']; res = j['result']
    t = np.array([r['t_cyc'] for r in rows])
    z = np.array([r['per_sub_z_A'] for r in rows])       # (F,5)
    dna = np.array([r['dna_z_A'] for r in rows]); dna = dna - dna[0]
    neng = np.array([r['n_engaged'] for r in rows])
    return dict(t=t, z=z, dna=dna, neng=neng, res=res)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default='../../outputs/php_cycle/ratchet')
    ap.add_argument('--out', default='../../outputs/php_cycle/ratchet')
    args = ap.parse_args()
    modes = [('concerted', 'ratchet_concerted'), ('sequential', 'ratchet_sequential')]
    data = {}
    for name, sub in modes:
        p = os.path.join(args.dir, sub)
        if os.path.exists(os.path.join(p, 'series.json')): data[name] = load(p)
    if not data:
        print("no ratchet runs found"); return

    fig, axes = plt.subplots(3, 2, figsize=(12, 8.5), sharex='col')
    for ci, (name, _sub) in enumerate(modes):
        if name not in data:
            for r in range(3): axes[r, ci].text(0.5, 0.5, f'{name}: (no data)', ha='center'); continue
        d = data[name]; net = d['res'].get('dna_net_A')
        a1, a2, a3 = axes[0, ci], axes[1, ci], axes[2, ci]
        for k in range(5): a1.plot(d['t'], d['z'][:, k], color=SUB_COLORS[k], lw=1.6, label=f'sub{k+1}')
        a1.set_title(f'{name.upper()} ratchet — per-subunit axial z(t)\n'
                     f'{"5 traces move TOGETHER" if name=="concerted" else "travelling wave / staggered"}', fontsize=9)
        a1.set_ylabel('subunit z (Å)'); a1.grid(alpha=0.2)
        if ci == 0: a1.legend(fontsize=7, ncol=5, loc='upper center')
        a2.plot(d['t'], d['dna'], color='#111', lw=2.0)
        a2.set_title(f'DNA axial displacement (net {net:+.1f} Å)', fontsize=9); a2.set_ylabel('DNA z (Å)'); a2.grid(alpha=0.2)
        a3.plot(d['t'], d['neng'], color='#276749', lw=1.8); a3.axhline(5, ls=':', c='gray'); a3.set_ylim(-0.2, 5.4)
        a3.set_title('M2 coupling n_engaged', fontsize=9); a3.set_ylabel('/5'); a3.set_xlabel('cycles'); a3.grid(alpha=0.2)
    fig.suptitle('Mechanochemical ratchet: CONCERTED vs SEQUENTIAL — the 3D-MINFLUX discriminator', fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(os.path.join(args.out, 'ratchet_compare.png'), dpi=150); plt.close(fig)
    summ = {m: {'dna_net_A': data[m]['res'].get('dna_net_A'),
                'n_eng_min': data[m]['res'].get('n_eng_min')} for m in data}
    json.dump(summ, open(os.path.join(args.out, 'ratchet_compare_summary.json'), 'w'), indent=2)
    print('wrote ratchet_compare.png ; summary:', json.dumps(summ))


if __name__ == '__main__':
    main()
