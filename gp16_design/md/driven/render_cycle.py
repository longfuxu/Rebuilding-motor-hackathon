#!/usr/bin/env python3
"""Concatenate a forward (P->H ascent) + reverse (H->P descent) ring+DNA trajectory into ONE continuous
P->H->P cycle animation. Fixed planar-reference ring axis (axial = vertical). Panels:
  LEFT  3D side view: subunit-centroid spheres + faint Ca traces + dsDNA backbone.
  RIGHT per-subunit axial z, M2 coupling, DNA axial displacement — all vs CYCLE PROGRESS (ascent | descent).
No ping-pong (the reverse half already returns the ring to planar). Real MD both halves.

Usage: python render_cycle.py --fwd <fwd_dir> --rev <rev_dir> --planar inputs/C_plus_dna_relaxed.pdb --out <dir>
"""
import argparse, os, json
import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

COPY_LOS = [1, 353, 705, 1057, 1409]; COPY_LEN = 342
SUB_COLORS = ['#e6194B', '#3cb44b', '#4363d8', '#f58231', '#911eb4']


def axis_from_planar(planar):
    pl = md.load(planar); top = pl.topology; xyz = pl.xyz[0] * 10.0
    sub = [top.select(f"name CA and chainid 0 and resSeq {lo} to {lo + COPY_LEN - 1}") for lo in COPY_LOS]
    cents = np.array([xyz[g].mean(0) for g in sub]); _, _, Vt = np.linalg.svd(cents - cents.mean(0))
    ax = Vt[2]
    if ax[2] < 0: ax = -ax
    return ax


def load(d, axis, nsub):
    t = md.load(os.path.join(d, 'traj.dcd'), top=os.path.join(d, 'final.pdb'))
    top = t.topology; X = t.xyz * 10.0
    ser = json.load(open(os.path.join(d, 'series.json')))['rows']
    sub = [top.select(f"name CA and chainid 0 and resSeq {lo} to {lo + COPY_LEN - 1}") for lo in COPY_LOS]
    ring = np.concatenate(sub)
    dna = np.concatenate([top.select(f"name P and chainid {ci}") for ci in range(1, top.n_chains)]) if top.n_chains > 1 else np.array([], int)
    F = min(X.shape[0], len(ser))
    Y = np.zeros((F, X.shape[1], 3)); subz = np.zeros((F, 5)); dnaz = np.zeros(F)
    e3 = axis; tmp = np.array([1.0, 0, 0]); tmp = tmp - tmp.dot(e3) * e3; e1 = tmp / np.linalg.norm(tmp); e2 = np.cross(e3, e1)
    R = np.column_stack([e1, e2, e3])
    for f in range(F):
        com = X[f, ring].mean(0); Y[f] = (X[f] - com) @ R
        for k, g in enumerate(sub): subz[f, k] = Y[f, g].mean(0)[2]
        if len(dna): dnaz[f] = Y[f, dna].mean(0)[2]
    neng = np.array([ser[i]['n_engaged'] for i in range(F)])
    return dict(Y=Y, subz=subz, dnaz=dnaz, neng=neng, sub=sub, dna_traces=[top.select(f"name P and chainid {ci}") for ci in range(1, top.n_chains)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fwd', required=True); ap.add_argument('--rev', required=True)
    ap.add_argument('--planar', default='inputs/C_plus_dna_relaxed.pdb')
    ap.add_argument('--out', default='../../outputs/php_cycle/cycle_campaign'); ap.add_argument('--nframes', type=int, default=100)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    axis = axis_from_planar(args.planar)
    F = load(args.fwd, axis, 5); Rv = load(args.rev, axis, 5)
    # concat cycle: forward frames then reverse frames
    Y = np.concatenate([F['Y'], Rv['Y']]); subz = np.concatenate([F['subz'], Rv['subz']])
    neng = np.concatenate([F['neng'], Rv['neng']])
    # DNA displacement PER HALF (each relative to its own start) — the reverse run seeds from a fixed helical
    # endpoint, not fwd's exact last frame, so a global cumulative would carry a spurious junction jump.
    dnaz = np.concatenate([F['dnaz'] - F['dnaz'][0], Rv['dnaz'] - Rv['dnaz'][0]])
    nF = len(F['Y']); N = len(Y); prog = np.arange(N) / (N - 1)
    sub = F['sub']; dna_traces = [s for s in F['dna_traces'] if len(s)]
    # subsample
    step = max(1, N // args.nframes); idx = np.arange(0, N, step)

    fig = plt.figure(figsize=(13, 6.2))
    gs = fig.add_gridspec(3, 2, width_ratios=[1.15, 1.0], hspace=0.5, wspace=0.22)
    ax3 = fig.add_subplot(gs[:, 0], projection='3d')
    aZ = fig.add_subplot(gs[0, 1]); aC = fig.add_subplot(gs[1, 1]); aD = fig.add_subplot(gs[2, 1])
    allY = Y[:, np.concatenate(sub)].reshape(-1, 3); mid = 0.5 * (allY.min(0) + allY.max(0)); rng = float((allY.max(0) - allY.min(0)).max()) * 0.6
    zlo, zhi = subz.min() - 3, subz.max() + 3

    def draw(fi):
        f = idx[fi]; ax3.clear(); aZ.clear(); aC.clear(); aD.clear()
        for k, g in enumerate(sub):
            p = Y[f, g]; ax3.plot(p[:, 0], p[:, 1], p[:, 2], color=SUB_COLORS[k], lw=0.6, alpha=0.3)
            c = p.mean(0); ax3.scatter([c[0]], [c[1]], [c[2]], color=SUB_COLORS[k], s=150, edgecolors='k', linewidths=0.7, depthshade=False)
        for tr in dna_traces:
            p = Y[f, tr]; ax3.plot(p[:, 0], p[:, 1], p[:, 2], color='#111', lw=3, alpha=0.9)
        ax3.set_xlim(mid[0]-rng, mid[0]+rng); ax3.set_ylim(mid[1]-rng, mid[1]+rng); ax3.set_zlim(zlo, zhi); ax3.set_axis_off()
        ax3.view_init(elev=6, azim=(-75 + 16*np.sin(fi/9.0)))
        phase = 'ASCENT (P→H, loading)' if f < nF else 'DESCENT (H→P, return to planar)'
        ax3.set_title(f"gp16 single-chain ring + dsDNA — full P→H→P cycle\n{phase}  ({prog[f]*100:3.0f}%)", fontsize=10)
        for panel in (aZ, aC, aD): panel.axvspan(0, 50, color='#2b6cb0', alpha=0.05); panel.axvspan(50, 100, color='#c05621', alpha=0.05); panel.axvline(prog[f]*100, color='gray', ls=':', lw=1)
        for k in range(5):
            aZ.plot(prog*100, subz[:, k], color=SUB_COLORS[k], lw=1.3, alpha=0.55); aZ.scatter([prog[f]*100], [subz[f, k]], color=SUB_COLORS[k], s=28, edgecolors='k', linewidths=0.4, zorder=5)
        aZ.set_ylim(zlo, zhi); aZ.set_title('per-subunit axial z — staircase up then down', fontsize=8.5); aZ.grid(alpha=0.2); aZ.set_xticklabels([])
        aC.plot(prog*100, neng, color='#276749', lw=1.8); aC.scatter([prog[f]*100], [neng[f]], color='#276749', s=34, zorder=5, edgecolors='k', linewidths=0.4)
        aC.axhline(5, ls=':', c='gray'); aC.set_ylim(-0.2, 5.4); aC.set_title('M2 coupling n_engaged', fontsize=8.5); aC.grid(alpha=0.2); aC.set_xticklabels([])
        aD.plot(prog*100, dnaz, color='#111', lw=1.8); aD.scatter([prog[f]*100], [dnaz[f]], color='#111', s=34, zorder=5)
        aD.set_ylabel('Å'); aD.set_title('DNA axial displacement', fontsize=8.5); aD.set_xlabel('cycle progress (%)  |  ascent · descent'); aD.grid(alpha=0.2)

    anim = FuncAnimation(fig, draw, frames=len(idx), interval=90)
    anim.save(os.path.join(args.out, 'cycle.gif'), writer=PillowWriter(fps=12)); print("wrote cycle.gif")
    try: anim.save(os.path.join(args.out, 'cycle.mp4'), writer='ffmpeg', fps=12, dpi=120); print("wrote cycle.mp4")
    except Exception as e: print("mp4 skipped:", repr(e)[:100])
    plt.close(fig)
    print(f"DNA per-half disp: ascent end {F['dnaz'][-1]-F['dnaz'][0]:+.2f} A ; descent end {Rv['dnaz'][-1]-Rv['dnaz'][0]:+.2f} A "
          f"(each rel. to its own start; junction NOT continuous across seeds)")


if __name__ == '__main__':
    main()
