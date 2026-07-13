#!/usr/bin/env python3
"""Deep analysis of the P->H->P cycle campaign (fwd_s* ascent + rev_s* descent, ring+dsDNA implicit MD).
Uses a FIXED planar-reference ring axis (from C_plus_dna_relaxed.pdb) for ALL runs so per-subunit axial z and
DNA axial position are directly comparable across ascent and descent. Produces:
  cycle_mechanics.png  — per-subunit z, M2 coupling, DNA axial displacement vs cycle progress (fwd then rev)
  cycle_stats.png      — unequal per-subunit treads (mean±sd over seeds) + forward/reverse steered work
  cycle_summary.json   — the numbers for the report
Caveat: the staircase CV drives all subunits concurrently, so this probes DRIVABILITY + COUPLING + geometric
treads, NOT concerted-vs-sequential ordering (that is CV-imposed) and NOT absolute rate/ΔG.
"""
import argparse, os, json, glob
import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

COPY_LOS = [1, 353, 705, 1057, 1409]; COPY_LEN = 342
SUB_COLORS = ['#e6194B', '#3cb44b', '#4363d8', '#f58231', '#911eb4']


def ring_axis_from(planar_top, planar_xyz):
    sub = [planar_top.select(f"name CA and chainid 0 and resSeq {lo} to {lo + COPY_LEN - 1}") for lo in COPY_LOS]
    cents = np.array([planar_xyz[g].mean(0) for g in sub]); ctr = cents.mean(0)
    _, _, Vt = np.linalg.svd(cents - ctr); ax = Vt[2]
    if ax[2] < 0: ax = -ax
    return ax, sub


def load_run(d, axis, sub_sel_names):
    """Return per-frame lam, per-subunit axial z (5), DNA axial pos, planarity, span, n_engaged, and W_total."""
    top_pdb = os.path.join(d, 'final.pdb'); dcd = os.path.join(d, 'traj.dcd')
    ser = json.load(open(os.path.join(d, 'series.json')))
    rows = ser['rows']; res = ser['result']
    t = md.load(dcd, top=top_pdb); top = t.topology; X = t.xyz * 10.0
    sub = [top.select(f"name CA and chainid 0 and resSeq {lo} to {lo + COPY_LEN - 1}") for lo in COPY_LOS]
    ring_all = np.concatenate(sub)
    dna = np.concatenate([top.select(f"name P and chainid {ci}") for ci in range(1, top.n_chains)]) \
        if top.n_chains > 1 else np.array([], int)
    F = min(X.shape[0], len(rows))
    lam = np.array([rows[i]['lam'] for i in range(F)])
    subz = np.zeros((F, 5)); dnaz = np.zeros(F)
    for f in range(F):
        com = X[f, ring_all].mean(0)
        for k, g in enumerate(sub): subz[f, k] = (X[f, g].mean(0) - com) @ axis
        if len(dna): dnaz[f] = (X[f, dna].mean(0) - com) @ axis
    plan = np.array([rows[i]['planarity_A'] for i in range(F)])
    span = np.array([rows[i]['axial_span_A'] for i in range(F)])
    neng = np.array([rows[i]['n_engaged'] for i in range(F)])
    return dict(lam=lam, subz=subz, dnaz=dnaz, plan=plan, span=span, neng=neng,
                W=res.get('W_total_kcal', np.nan), dna_net=float(dnaz[-1] - dnaz[0]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default='../../outputs/php_cycle/cycle_campaign')
    ap.add_argument('--planar', default='inputs/C_plus_dna_relaxed.pdb')
    ap.add_argument('--out', default='../../outputs/php_cycle/cycle_campaign')
    args = ap.parse_args()
    pl = md.load(args.planar); axis, _ = ring_axis_from(pl.topology, pl.xyz[0] * 10.0)
    print("fixed ring axis:", np.round(axis, 3))

    fwd = [load_run(d, axis, None) for d in sorted(glob.glob(os.path.join(args.dir, 'fwd_s*'))) if os.path.exists(os.path.join(d, 'series.json'))]
    rev = [load_run(d, axis, None) for d in sorted(glob.glob(os.path.join(args.dir, 'rev_s*'))) if os.path.exists(os.path.join(d, 'series.json'))]
    print(f"loaded fwd={len(fwd)} rev={len(rev)} runs")
    if not fwd and not rev:
        print("no runs found"); return

    S = {}
    # per-subunit treads (fwd endpoint z), mean±sd over seeds
    if fwd:
        treads = np.array([r['subz'][-1] - r['subz'][0] for r in fwd])  # (nseed,5) final-minus-initial axial
        S['tread_mean_A'] = treads.mean(0).tolist(); S['tread_sd_A'] = treads.std(0).tolist()
        S['coupling_fwd_min'] = float(np.mean([r['neng'].min() for r in fwd]))
        S['W_fwd_mean'] = float(np.nanmean([r['W'] for r in fwd])); S['W_fwd_sd'] = float(np.nanstd([r['W'] for r in fwd]))
        S['dna_net_fwd_A'] = float(np.mean([r['dna_net'] for r in fwd])); S['dna_net_fwd_sd'] = float(np.std([r['dna_net'] for r in fwd]))
    if rev:
        S['coupling_rev_min'] = float(np.mean([r['neng'].min() for r in rev]))
        S['W_rev_mean'] = float(np.nanmean([r['W'] for r in rev])); S['W_rev_sd'] = float(np.nanstd([r['W'] for r in rev]))
        S['dna_net_rev_A'] = float(np.mean([r['dna_net'] for r in rev])); S['dna_net_rev_sd'] = float(np.std([r['dna_net'] for r in rev]))
        S['plan_return'] = float(np.mean([r['plan'][-1] for r in rev]))  # does the ring return to planar?
    json.dump(S, open(os.path.join(args.out, 'cycle_summary.json'), 'w'), indent=2)

    # --- Figure A: cycle mechanics (representative fwd + rev; cycle progress 0..1) ---
    fig, (a1, a2, a3) = plt.subplots(3, 1, figsize=(8, 8.5), sharex=True)
    def cyc_x(r, half):  # half 0 = ascent (0..0.5), 1 = descent (0.5..1)
        p = (r['lam'] - r['lam'][0]) / (r['lam'][-1] - r['lam'][0] + 1e-9)
        return 0.5 * p if half == 0 else 0.5 + 0.5 * p
    for r in fwd:
        x = cyc_x(r, 0)
        for k in range(5): a1.plot(x, r['subz'][:, k], color=SUB_COLORS[k], lw=1.0, alpha=0.35)
        a2.plot(x, r['neng'], color='#276749', lw=1.0, alpha=0.35)
        a3.plot(x, r['dnaz'] - r['dnaz'][0], color='#111', lw=1.0, alpha=0.35)
    for r in rev:
        x = cyc_x(r, 1)
        for k in range(5): a1.plot(x, r['subz'][:, k], color=SUB_COLORS[k], lw=1.0, alpha=0.35)
        a2.plot(x, r['neng'], color='#276749', lw=1.0, alpha=0.35)
        a3.plot(x, r['dnaz'] - r['dnaz'][0], color='#111', lw=1.0, alpha=0.35)
    for k in range(5): a1.plot([], [], color=SUB_COLORS[k], label=f'subunit {k+1}')
    a1.axvline(0.5, ls=':', c='gray'); a2.axvline(0.5, ls=':', c='gray'); a3.axvline(0.5, ls=':', c='gray')
    a1.set_ylabel('subunit axial z (Å)'); a1.legend(fontsize=7, ncol=5, loc='upper center')
    a1.set_title('P→H→P cycle (left half = ascent/loading, right half = descent/power stroke)\n'
                 'per-subunit axial z: graded staircase up then back down', fontsize=9)
    a2.set_ylabel('M2 engaged /5'); a2.set_ylim(-0.2, 5.4); a2.set_title('M2 coupling across the cycle', fontsize=9)
    a3.set_ylabel('DNA axial disp (Å)'); a3.set_xlabel('cycle progress'); a3.set_title('DNA axial displacement across the cycle', fontsize=9)
    for a in (a1, a2, a3): a.grid(alpha=0.2)
    fig.tight_layout(); fig.savefig(os.path.join(args.out, 'cycle_mechanics.png'), dpi=150); plt.close(fig)

    # --- Figure B: stats (treads ± sd; work fwd/rev) ---
    fig, (b1, b2) = plt.subplots(1, 2, figsize=(11, 4))
    if fwd:
        xk = np.arange(5)
        b1.bar(xk, S['tread_mean_A'], yerr=S['tread_sd_A'], color=SUB_COLORS, capsize=4, edgecolor='k')
        b1.axhline(0, color='k', lw=0.6); b1.set_xticks(xk); b1.set_xticklabels([f'sub{ k+1}' for k in range(5)])
        b1.set_ylabel('axial tread Δz over ascent (Å)')
        b1.set_title(f'Unequal per-subunit treads (mean±sd, n={len(fwd)})\ngraded, position-dependent — NOT 5 equal steps', fontsize=9)
    lbls = []; vals = []; errs = []
    if fwd: lbls.append('forward\n(ascent)'); vals.append(S['W_fwd_mean']); errs.append(S['W_fwd_sd'])
    if rev: lbls.append('reverse\n(descent)'); vals.append(S['W_rev_mean']); errs.append(S['W_rev_sd'])
    b2.bar(range(len(vals)), vals, yerr=errs, color=['#2b6cb0', '#c05621'][:len(vals)], capsize=4, edgecolor='k')
    b2.set_xticks(range(len(lbls))); b2.set_xticklabels(lbls); b2.set_ylabel('steered work (kcal/mol)')
    b2.set_title('Steered work per half-stroke (dissipative upper bound, NOT ΔG)', fontsize=9)
    for a in (b1, b2): a.grid(alpha=0.2, axis='y')
    fig.tight_layout(); fig.savefig(os.path.join(args.out, 'cycle_stats.png'), dpi=150); plt.close(fig)

    print(json.dumps(S, indent=2))


if __name__ == '__main__':
    main()
