#!/usr/bin/env python3
"""Analyze the A100 umbrella-sampling output (window_data.json) into the Goal-4 answers:
  (a) is H a reachable minimum without unfolding      -> shape of G(xi), G(H)-G(P), barrier
  (b) does M2 (R146->neighbour Walker-A) coupling hold -> n_engaged vs xi
  (c) pathway concerted vs sequential, equal vs unequal-> per-subunit axial z vs xi
  (d) reversible?                                       -> forward-vs-reverse G (hysteresis)

Runs local (Mac), free. Partial-safe: works on however many windows are present. WHAM is done
per direction here (the driver's built-in pmf.json is forward-only). Writes figures + summary md.

Usage: python analyze_umbrella.py [--data outputs/php_cycle/C3_umbrella/umb/window_data.json]
                                  [--out outputs/php_cycle/C3_umbrella]
"""
import argparse, json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

kT = 0.0019872041 * 300.0  # Boltzmann kcal/mol/K * 300 K = 0.5962 kcal/mol


def wham(windows, k_eff, nbins=61, niter=2000):
    """1-D WHAM over a set of windows (each has 'xi' sample list and umbrella center 'lam')."""
    allxi = np.concatenate([np.asarray(w['xi'], float) for w in windows])
    edges = np.linspace(allxi.min() - .05, allxi.max() + .05, nbins + 1)
    ctrs = 0.5 * (edges[:-1] + edges[1:])
    hist = np.array([np.histogram(w['xi'], edges)[0] for w in windows], float)
    N = hist.sum(1)
    centers = np.array([w['lam'] for w in windows])          # umbrella minima are at xi=lam
    bias = 0.5 * k_eff * (ctrs[None, :] - centers[:, None]) ** 2 / kT   # dimensionless
    F = np.zeros(len(windows)); num = hist.sum(0)
    for _ in range(niter):
        denom = (N[:, None] * np.exp(F[:, None] - bias)).sum(0)
        P = np.where(denom > 0, num / denom, 0.0)
        s = P.sum();  P = P / s if s > 0 else P
        Fn = -np.log((P[None, :] * np.exp(-bias)).sum(1) + 1e-300)
        if np.max(np.abs(Fn - F)) < 1e-7:
            F = Fn; break
        F = Fn
    G = -kT * np.log(P + 1e-300)
    G -= np.nanmin(G[np.isfinite(G)])
    return ctrs, G


def basin_stats(ctrs, G):
    gP = float(G[np.argmin(np.abs(ctrs - 0.0))])
    gH = float(G[np.argmin(np.abs(ctrs - 1.0))])
    mid = (ctrs > 0.1) & (ctrs < 0.9) & np.isfinite(G)
    barrier = float(np.nanmax(G[mid]) - min(gP, gH)) if mid.any() else float('nan')
    xi_barrier = float(ctrs[mid][np.argmax(G[mid])]) if mid.any() else float('nan')
    return gP, gH, barrier, xi_barrier


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', default='../../outputs/php_cycle/C3_umbrella/umb/window_data.json')
    ap.add_argument('--out', default='../../outputs/php_cycle/C3_umbrella')
    args = ap.parse_args()
    d = json.load(open(args.data))
    k_eff = d['k_eff']; W = d['windows']
    os.makedirs(args.out, exist_ok=True)
    fwd = [w for w in W if w.get('pass_') == 'fwd']
    rev = [w for w in W if w.get('pass_') == 'rev']
    print(f"loaded {len(W)} windows: fwd={len(fwd)} rev={len(rev)} k_eff={k_eff:.0f}")

    # ---- PMF (a) + reversibility (d) ----
    fig, ax = plt.subplots(figsize=(6, 4.2))
    summ = {'k_eff': k_eff, 'n_fwd': len(fwd), 'n_rev': len(rev)}
    for name, ws, col in [('forward', fwd, '#2b6cb0'), ('reverse', rev, '#c05621')]:
        if len(ws) < 3:
            continue
        ctrs, G = wham(ws, k_eff)
        gP, gH, bar, xib = basin_stats(ctrs, G)
        summ[name] = dict(G_P=gP, G_H=gH, dG_H_minus_P=gH - gP, barrier=bar, xi_barrier=xib)
        ax.plot(ctrs, G, '-', color=col, lw=2, label=f"{name}: ΔG(H−P)={gH-gP:+.1f}, barrier≈{bar:.1f} kcal/mol")
        np.savetxt(os.path.join(args.out, f'pmf_{name}.csv'),
                   np.column_stack([ctrs, G]), header='xi,G_kcal', delimiter=',', comments='')
    ax.axvline(0, ls=':', c='gray'); ax.axvline(1, ls=':', c='gray')
    ax.set_xlabel('ξ  (0 = planar / P,  1 = 7JQQ helical staircase / H)')
    ax.set_ylabel('G(ξ)  (kcal/mol)')
    ax.set_title('gp16 single-chain ring — P→H free-energy landscape (umbrella/WHAM)')
    ax.legend(fontsize=8, loc='upper center'); fig.tight_layout()
    fig.savefig(os.path.join(args.out, 'pmf.png'), dpi=150); plt.close(fig)

    # ---- coupling (b): n_engaged vs xi (forward) ----
    if fwd:
        xs = [w['xi_mean'] for w in fwd]; ne = [w['n_engaged_mean'] for w in fwd]
        fig, ax = plt.subplots(figsize=(6, 3.6))
        ax.plot(xs, ne, 'o-', color='#276749'); ax.set_ylim(-0.2, 5.2)
        ax.set_xlabel('ξ'); ax.set_ylabel('M2 interfaces engaged (R146→WalkerA <8Å)')
        ax.set_title('M2 coupling vs P→H progress'); fig.tight_layout()
        fig.savefig(os.path.join(args.out, 'coupling_vs_xi.png'), dpi=150); plt.close(fig)
        summ['coupling'] = dict(n_engaged_at_P=ne[0], n_engaged_at_H=ne[-1],
                                n_engaged_min=float(np.min(ne)))

    # ---- pathway (c): per-subunit axial z vs xi (forward) ----
    if fwd:
        xs = np.array([w['xi_mean'] for w in fwd])
        Z = np.array([w['per_sub_axial_nm'] for w in fwd]) * 10.0  # -> Angstrom
        order = np.argsort(xs); xs = xs[order]; Z = Z[order]
        fig, ax = plt.subplots(figsize=(6, 4.0))
        for k in range(Z.shape[1]):
            ax.plot(xs, Z[:, k], 'o-', ms=3, label=f'subunit {k+1}')
        ax.set_xlabel('ξ'); ax.set_ylabel('axial z of subunit centroid (Å, rel. ring COM)')
        ax.set_title('Per-subunit axial displacement vs ξ  (concerted vs sequential)')
        ax.legend(fontsize=8, ncol=5, loc='lower center'); fig.tight_layout()
        fig.savefig(os.path.join(args.out, 'per_subunit_z.png'), dpi=150); plt.close(fig)
        # crude concerted-vs-sequential metric: spread of the xi at which each subunit reaches half its travel
        halfxi = []
        for k in range(Z.shape[1]):
            zk = Z[:, k]; lo, hi = zk[0], zk[-1]
            if abs(hi - lo) < 1e-6: halfxi.append(float('nan')); continue
            frac = (zk - lo) / (hi - lo)
            idx = int(np.argmin(np.abs(frac - 0.5)))
            halfxi.append(float(xs[idx]))
        summ['pathway'] = dict(per_subunit_half_xi=halfxi,
                               half_xi_spread=float(np.nanmax(halfxi) - np.nanmin(halfxi)),
                               travel_A=[round(float(Z[-1, k] - Z[0, k]), 2) for k in range(Z.shape[1])])

    json.dump(summ, open(os.path.join(args.out, 'analysis_summary.json'), 'w'), indent=2)
    # ---- markdown ----
    L = ['# Umbrella-sampling analysis — gp16 single-chain P→H landscape', '',
         f"windows: forward={len(fwd)}, reverse={len(rev)}; k_eff={k_eff:.0f} kJ/mol/nm²; kT={kT:.3f} kcal/mol", '']
    if 'forward' in summ:
        f = summ['forward']
        L += ['## (a) Is H a reachable minimum? — G(ξ) shape',
              f"- forward: **ΔG(H−P) = {f['dG_H_minus_P']:+.1f} kcal/mol**, barrier ≈ **{f['barrier']:.1f} kcal/mol** at ξ≈{f['xi_barrier']:.2f}",
              f"  - G(P)={f['G_P']:.1f}, G(H)={f['G_H']:.1f} kcal/mol", '']
    if 'reverse' in summ:
        r = summ['reverse']
        L += ['## (d) Reversible? — forward vs reverse (hysteresis)',
              f"- reverse: ΔG(H−P) = {r['dG_H_minus_P']:+.1f} kcal/mol, barrier ≈ {r['barrier']:.1f} kcal/mol",
              f"- hysteresis in ΔG(H−P): {abs(summ.get('forward',{}).get('dG_H_minus_P',float('nan')) - r['dG_H_minus_P']):.1f} kcal/mol", '']
    if 'coupling' in summ:
        c = summ['coupling']
        L += ['## (b) Does M2 coupling survive?',
              f"- n_engaged: P={c['n_engaged_at_P']:.1f} → H={c['n_engaged_at_H']:.1f} (min along path {c['n_engaged_min']:.1f} of 5)", '']
    if 'pathway' in summ:
        p = summ['pathway']
        L += ['## (c) Concerted vs sequential / equal vs unequal',
              f"- per-subunit half-travel ξ: {[round(x,2) for x in p['per_subunit_half_xi']]}",
              f"- spread = **{p['half_xi_spread']:.2f}** (≈0 ⇒ concerted; large ⇒ sequential/staggered)",
              f"- per-subunit axial travel (Å): {p['travel_A']}  (unequal ⇒ graded treads → unequal substeps)", '']
    L += ['## Figures', '- `pmf.png` — G(ξ) forward (+reverse)',
          '- `coupling_vs_xi.png` — M2 engagement vs ξ', '- `per_subunit_z.png` — per-subunit axial z vs ξ', '',
          '_Caveat: idealized monotonic staircase CV; implicit solvent; 7JQQ is a 3-ATP partial-occupancy '
          'intermediate (~4.8 Å), not the full ~35 Å stroke. This is the relative-thermodynamics/mechanism '
          'landscape, NOT an absolute rate._']
    open(os.path.join(args.out, 'C3_UMBRELLA_RESULTS.md'), 'w').write('\n'.join(L))
    print('wrote:', args.out, '-> pmf.png, coupling_vs_xi.png, per_subunit_z.png, analysis_summary.json, C3_UMBRELLA_RESULTS.md')
    print(json.dumps(summ, indent=1))


if __name__ == '__main__':
    main()
