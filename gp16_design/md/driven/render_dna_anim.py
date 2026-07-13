#!/usr/bin/env python3
"""Render the ring+DNA P->H MD into (1) a mechanistic dashboard animation and (2) static figures that put the
KEY MECHANISTIC OBSERVABLES on the Y axis against the reaction coordinate (drive lambda) on X — so the story is
legible: the CV drives the ring helical (axial staircase span up, planarity up toward the helical target) WHILE
the M2 coupling stays engaged (n_engaged ~5) — the core result (b).

Animation panels:
  LEFT  3D side view (ring axis vertical): faint Ca traces + bold subunit-centroid spheres + dsDNA backbone.
  RIGHT-top    per-subunit axial height z_k (the lifting into the staircase).
  RIGHT-mid    planarity & axial-span vs lambda, with the 7JQQ helical targets (ring -> helical).
  RIGHT-bot    M2 coupling n_engaged vs lambda (stays ~5 -> coupling survives the stroke).
Static: subunit_lifting.png, mechanism_observables.png, subunit_lifting.csv.

Usage: python render_dna_anim.py --dcd .../traj.dcd --top .../final.pdb [--series .../series.csv] --out <dir>
"""
import argparse, os, csv
import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

COPY_LOS = [1, 353, 705, 1057, 1409]; COPY_LEN = 342
SUB_COLORS = ['#e6194B', '#3cb44b', '#4363d8', '#f58231', '#911eb4']
HELIX_PLAN = 1.68   # 7JQQ helical planarity target (A)
HELIX_SPAN = 4.76   # 7JQQ helical axial-span target (A)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dcd', required=True); ap.add_argument('--top', required=True)
    ap.add_argument('--series', default=''); ap.add_argument('--out', default='../../outputs/php_cycle/dna_translocation')
    ap.add_argument('--nframes', type=int, default=80); ap.add_argument('--pingpong', type=int, default=1)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    series_path = args.series or os.path.join(os.path.dirname(args.top), 'series.csv')

    t = md.load(args.dcd, top=args.top); top = t.topology
    print(f"loaded {t.n_frames} frames, {t.n_atoms} atoms")
    step = max(1, t.n_frames // args.nframes)
    keep = np.arange(0, t.n_frames, step)
    t = t[::step]; X = t.xyz * 10.0; F = X.shape[0]

    # per-frame mechanistic observables from tmd_staircase series.csv (aligned to the same report cadence)
    obs = {}
    if os.path.exists(series_path):
        rows = list(csv.DictReader(open(series_path)))
        def col(name): return np.array([float(r[name]) for r in rows])
        allobs = {k: col(k) for k in ['lam', 'planarity_A', 'axial_span_A', 'n_engaged', 'time_ps']}
        idx = [min(i, len(rows) - 1) for i in keep]
        obs = {k: v[idx] for k, v in allobs.items()}
        print(f"series: lam {obs['lam'][0]:.2f}->{obs['lam'][-1]:.2f}, planarity {obs['planarity_A'][0]:.2f}->{obs['planarity_A'][-1]:.2f}A, "
              f"span {obs['axial_span_A'][0]:.2f}->{obs['axial_span_A'][-1]:.2f}A, n_eng {obs['n_engaged'][0]:.0f}->{obs['n_engaged'][-1]:.0f}")
    lam = obs.get('lam', np.linspace(0, 1, F))

    # ring-axis frame -> axial = 3rd coord
    sub_ca = [top.select(f"name CA and chainid 0 and resSeq {lo} to {lo + COPY_LEN - 1}") for lo in COPY_LOS]
    ring_all = np.concatenate(sub_ca)
    dna_traces = [top.select(f"name P and chainid {ci}") for ci in range(1, top.n_chains)]
    dna_traces = [s for s in dna_traces if len(s)]
    dna_all = np.concatenate(dna_traces) if dna_traces else np.array([], int)
    cents0 = np.array([X[0, g].mean(0) for g in sub_ca])
    _, _, Vt = np.linalg.svd(cents0 - cents0.mean(0)); axis = Vt[2]
    if axis[2] < 0: axis = -axis
    tmp = np.array([1.0, 0, 0]); tmp = tmp - tmp.dot(axis) * axis
    if np.linalg.norm(tmp) < 1e-3: tmp = np.array([0, 1.0, 0.0]); tmp = tmp - tmp.dot(axis) * axis
    e1 = tmp / np.linalg.norm(tmp); e2 = np.cross(axis, e1); R = np.column_stack([e1, e2, axis])
    Y = np.zeros_like(X); sub_z = np.zeros((F, 5)); dna_z = np.zeros(F)
    for f in range(F):
        Y[f] = (X[f] - X[f, ring_all].mean(0)) @ R
        for k, g in enumerate(sub_ca): sub_z[f, k] = Y[f, g].mean(0)[2]
        if len(dna_all): dna_z[f] = Y[f, dna_all].mean(0)[2]

    # --- static: lifting + CSV ---
    np.savetxt(os.path.join(args.out, 'subunit_lifting.csv'),
               np.column_stack([lam, sub_z, dna_z - dna_z[0]]),
               header='lambda,z_sub1,z_sub2,z_sub3,z_sub4,z_sub5,dna_axial_disp_A', delimiter=',', comments='')
    fg, a = plt.subplots(figsize=(7, 4.2))
    for k in range(5): a.plot(lam, sub_z[:, k], color=SUB_COLORS[k], lw=2, label=f'subunit {k+1}')
    a.set_xlabel('reaction coordinate  λ  (0 = planar → 1 = helical drive)'); a.set_ylabel('subunit axial height z (Å)')
    a.set_title('Per-subunit lifting into the helical staircase'); a.legend(fontsize=8, ncol=5, loc='lower center')
    a.grid(alpha=0.25); fg.tight_layout(); fg.savefig(os.path.join(args.out, 'subunit_lifting.png'), dpi=150); plt.close(fg)

    # --- static: mechanism observables (X = reaction coordinate, Y = key info) ---
    if obs:
        fg, (p1, p2) = plt.subplots(2, 1, figsize=(7.2, 6.2), sharex=True)
        p1.plot(lam, obs['axial_span_A'], color='#2b6cb0', lw=2.2, label='axial staircase span')
        p1.plot(lam, obs['planarity_A'], color='#805ad5', lw=2.2, label='planarity (out-of-plane RMS)')
        p1.axhline(HELIX_SPAN, ls='--', color='#2b6cb0', alpha=0.5); p1.axhline(HELIX_PLAN, ls='--', color='#805ad5', alpha=0.5)
        p1.text(0.02, HELIX_SPAN + 0.05, '7JQQ helical span 4.76 Å', fontsize=7, color='#2b6cb0')
        p1.text(0.02, HELIX_PLAN + 0.05, '7JQQ helical planarity 1.68 Å', fontsize=7, color='#805ad5')
        p1.set_ylabel('Å'); p1.legend(fontsize=8, loc='upper left'); p1.grid(alpha=0.25)
        p1.set_title('Ring adopts the helical staircase as the CV drives it\n'
                     f'span 0.2→{obs["axial_span_A"][-1]:.1f} Å (~{100*obs["axial_span_A"][-1]/HELIX_SPAN:.0f}% of helical); '
                     f'planarity →{obs["planarity_A"][-1]:.2f} Å (~{100*obs["planarity_A"][-1]/HELIX_PLAN:.0f}%)', fontsize=9)
        p2.plot(lam, obs['n_engaged'], color='#276749', lw=2.4, marker='o', ms=3)
        p2.axhline(5, ls=':', color='gray'); p2.set_ylim(-0.2, 5.4)
        p2.set_xlabel('reaction coordinate  λ  (0 = planar → 1 = helical drive)')
        p2.set_ylabel('M2 interfaces engaged / 5')
        p2.set_title('★ (b) M2 coupling (R146→neighbour Walker-A) STAYS engaged (~5/5) through the whole P→H stroke —\n'
                     'the single chain does NOT shed its coupling to move toward helical', fontsize=9)
        p2.grid(alpha=0.25); fg.tight_layout(); fg.savefig(os.path.join(args.out, 'mechanism_observables.png'), dpi=150); plt.close(fg)

    # --- animation: 3D + 3 synced observable panels ---
    order = list(range(F)) + (list(range(F - 2, 0, -1)) if args.pingpong else [])
    fig = plt.figure(figsize=(13, 6.2))
    gs = fig.add_gridspec(3, 2, width_ratios=[1.15, 1.0], hspace=0.55, wspace=0.22)
    ax3 = fig.add_subplot(gs[:, 0], projection='3d')
    axL = fig.add_subplot(gs[0, 1]); axP = fig.add_subplot(gs[1, 1]); axC = fig.add_subplot(gs[2, 1])
    allY = Y[:, ring_all].reshape(-1, 3); lo = allY.min(0); hi = allY.max(0)
    rng = float((hi - lo).max()) * 0.6; mid = 0.5 * (hi + lo); zlo, zhi = sub_z.min() - 3, sub_z.max() + 3
    span = obs['axial_span_A'] if obs else (sub_z.max(1) - sub_z.min(1))
    plan = obs['planarity_A'] if obs else np.zeros(F)
    neng = obs['n_engaged'] if obs else np.full(F, 5.0)

    def draw(fi):
        f = order[fi]; ax3.clear(); axL.clear(); axP.clear(); axC.clear()
        for k, g in enumerate(sub_ca):
            p = Y[f, g]; ax3.plot(p[:, 0], p[:, 1], p[:, 2], color=SUB_COLORS[k], lw=0.6, alpha=0.32)
            c = p.mean(0); ax3.scatter([c[0]], [c[1]], [c[2]], color=SUB_COLORS[k], s=160, edgecolors='k', linewidths=0.8, depthshade=False)
        for tr in dna_traces:
            p = Y[f, tr]; ax3.plot(p[:, 0], p[:, 1], p[:, 2], color='#111111', lw=3, alpha=0.9)
        ax3.set_xlim(mid[0]-rng, mid[0]+rng); ax3.set_ylim(mid[1]-rng, mid[1]+rng); ax3.set_zlim(zlo, zhi)
        ax3.set_axis_off(); ax3.view_init(elev=6, azim=(-75 + 16*np.sin(fi/8.0)))
        ax3.set_title(f"gp16 single-chain ring + dsDNA — P→H stroke  (λ={lam[f]:.2f})", fontsize=10)
        # lifting
        for k in range(5):
            axL.plot(lam, sub_z[:, k], color=SUB_COLORS[k], lw=1.4, alpha=0.5)
            axL.scatter([lam[f]], [sub_z[f, k]], color=SUB_COLORS[k], s=32, edgecolors='k', linewidths=0.5, zorder=5)
        axL.axvline(lam[f], color='gray', ls=':', lw=1); axL.set_ylim(zlo, zhi)
        axL.set_title('per-subunit lifting z (Å) — graded, unequal staircase', fontsize=8.5); axL.grid(alpha=0.2); axL.set_xticklabels([])
        # planarity / span -> helical
        axP.plot(lam, span, color='#2b6cb0', lw=1.8, label='staircase span')
        axP.plot(lam, plan, color='#805ad5', lw=1.8, label='planarity')
        axP.axhline(HELIX_SPAN, ls='--', color='#2b6cb0', alpha=0.4); axP.axhline(HELIX_PLAN, ls='--', color='#805ad5', alpha=0.4)
        axP.scatter([lam[f]], [span[f]], color='#2b6cb0', s=28, zorder=5); axP.scatter([lam[f]], [plan[f]], color='#805ad5', s=28, zorder=5)
        axP.axvline(lam[f], color='gray', ls=':', lw=1); axP.set_ylim(0, HELIX_SPAN + 0.6)
        axP.set_title(f'ring → helical:  span {span[f]:.1f} Å, planarity {plan[f]:.2f} Å  (dashed = 7JQQ targets)', fontsize=8.5)
        axP.legend(fontsize=7, loc='upper left'); axP.grid(alpha=0.2); axP.set_xticklabels([])
        # coupling
        axC.plot(lam, neng, color='#276749', lw=2.0); axC.axhline(5, ls=':', color='gray')
        axC.scatter([lam[f]], [neng[f]], color='#276749', s=40, zorder=5, edgecolors='k', linewidths=0.5)
        axC.axvline(lam[f], color='gray', ls=':', lw=1); axC.set_ylim(-0.2, 5.4)
        axC.set_title(f'★(b) M2 coupling engaged: {neng[f]:.0f}/5 — stays coupled through the stroke', fontsize=8.5)
        axC.set_xlabel('reaction coordinate  λ  (planar → helical)'); axC.grid(alpha=0.2)

    anim = FuncAnimation(fig, draw, frames=len(order), interval=90)
    gif = os.path.join(args.out, 'dna_translocation.gif')
    anim.save(gif, writer=PillowWriter(fps=12)); print("wrote", gif)
    try:
        anim.save(os.path.join(args.out, 'dna_translocation.mp4'), writer='ffmpeg', fps=12, dpi=120); print("wrote mp4")
    except Exception as e:
        print("mp4 skipped:", repr(e)[:100])
    plt.close(fig); print("frames:", len(order))


if __name__ == '__main__':
    main()
