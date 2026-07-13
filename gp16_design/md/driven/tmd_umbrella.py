#!/usr/bin/env python3
"""Track B (reliable-timescale version) — UMBRELLA SAMPLING of the design gp16 ring along the
per-subunit AXIAL STAIRCASE reaction coordinate, to get the P<->H free-energy landscape G(xi)
WITHOUT simulating the real 10-80 ms (impossible for MD). Answers, as MECHANISM + relative
thermodynamics (NOT the absolute rate):
  (a) is H a reachable minimum from P without unfolding   -> shape of G(xi)
  (b) does M2 (R146->neighbour Walker-A) coupling survive  -> n_engaged per window
  (c) pathway: concerted vs sequential, equal vs unequal   -> per-subunit axial z vs xi
  (d) reversible?                                           -> forward- vs reverse-seeded G (hysteresis)

Method: reuse the proven CustomCentroidBondForce staircase CV, but HOLD lam fixed per window
(umbrella) instead of ramping. Windows seeded sequentially (i from i-1's end) for overlap.
Preemption-safe: window_data.json + ckpt.chk written after EVERY window. Implicit GBSA-OBC2 ~28k atoms.
xi = (sum_k proj_k*z_k)/(sum_k z_k^2)  (0=planar, 1=full 7JQQ staircase); k_eff = kcv*sum_k z_k^2.
"""
import argparse, os, time, json
import numpy as np
from tmd_staircase import log, COPY_LOS, COPY_LEN, HELIX_AXIAL_SPAN, ca_groups, m2_indices, ring_geometry

def m2_min_dists(P, guan, walk):
    d = []
    for k in range(5):
        nb = (k + 1) % 5
        if not guan[k] or not walk[nb]: d.append(99.0); continue
        A = P[guan[k]][:, None, :]; B = P[walk[nb]][None, :, :]
        d.append(float(np.sqrt(((A - B) ** 2).sum(-1)).min()))
    return d

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', default='C_start.pdb')
    ap.add_argument('--out', required=True)
    ap.add_argument('--nwin', type=int, default=25)
    ap.add_argument('--lam_min', type=float, default=-0.10)
    ap.add_argument('--lam_max', type=float, default=1.10)
    ap.add_argument('--kcv', type=float, default=30000.0)
    ap.add_argument('--equil_ps', type=float, default=400.0)
    ap.add_argument('--sample_ps', type=float, default=1600.0)
    ap.add_argument('--report_ps', type=float, default=4.0)
    ap.add_argument('--timestep_fs', type=float, default=3.0)
    ap.add_argument('--both', type=int, default=1)
    ap.add_argument('--passes', default='auto', choices=['auto', 'fwd', 'rev', 'both'],
                    help="'auto' honors --both; 'fwd'/'rev' run one direction; 'both' = fwd then rev")
    ap.add_argument('--resume_state', default='',
                    help='portable OpenMM State XML to resume from (e.g. state_fwd_end.xml); skips minimize/equil. '
                         'Use with --passes rev to run the reverse pass seeded from the forward endpoint.')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--gb_cutoff_nm', type=float, default=2.0)
    ap.add_argument('--platform', default='auto')
    args = ap.parse_args()
    from openmm import app, unit, Platform, LangevinMiddleIntegrator, CustomCentroidBondForce
    os.makedirs(args.out, exist_ok=True)
    json.dump(vars(args), open(os.path.join(args.out, 'params.json'), 'w'), indent=1)

    pdb = app.PDBFile(args.input); top, pos = pdb.topology, pdb.positions
    groups, m = ca_groups(top); guan, walk = m2_indices(m)
    all_ca = [i for g in groups for i in g]
    P0 = np.array(pos.value_in_unit(unit.angstrom))
    cents0 = np.array([P0[g].mean(0) for g in groups]); ctr0 = cents0.mean(0)
    _, _, Vt = np.linalg.svd(cents0 - ctr0); n_axis = Vt[2]
    if n_axis[2] < 0: n_axis = -n_axis
    amp = HELIX_AXIAL_SPAN / 2.0
    z_target_nm = (np.linspace(-amp, amp, 5)) / 10.0
    sumz2 = float((z_target_nm ** 2).sum()); k_eff = args.kcv * sumz2
    log(f"axis={np.round(n_axis,3)} z_target_nm={np.round(z_target_nm,3)} sumz2={sumz2:.4f} k_eff={k_eff:.0f}")

    ff = app.ForceField('amber14-all.xml', 'implicit/obc2.xml')
    kw = dict(constraints=app.HBonds, soluteDielectric=1.0, solventDielectric=78.5,
              hydrogenMass=1.5 * unit.amu,
              nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=args.gb_cutoff_nm * unit.nanometer)
    system = ff.createSystem(top, **kw)
    cf = CustomCentroidBondForce(2, "0.5*kcv*(nx*(x1-x2)+ny*(y1-y2)+nz*(z1-z2) - lam*zt)^2")
    cf.addGlobalParameter('kcv', args.kcv); cf.addGlobalParameter('lam', 0.0)
    for nm, v in [('nx', n_axis[0]), ('ny', n_axis[1]), ('nz', n_axis[2])]: cf.addGlobalParameter(nm, float(v))
    cf.addPerBondParameter('zt')
    gidx = [cf.addGroup([int(i) for i in g]) for g in groups]; gcom = cf.addGroup([int(i) for i in all_ca])
    for k in range(5): cf.addBond([gidx[k], gcom], [float(z_target_nm[k])])
    system.addForce(cf)

    integ = LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond, args.timestep_fs * unit.femtoseconds)
    integ.setRandomNumberSeed(args.seed)
    avail = {Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())}
    pn = args.platform if args.platform != 'auto' else next(p for p in ['CUDA','OpenCL','CPU','Reference'] if p in avail)
    try: sim = app.Simulation(top, system, integ, Platform.getPlatformByName(pn), {'Precision': 'mixed'})
    except Exception: sim = app.Simulation(top, system, integ, Platform.getPlatformByName(pn))
    log(f"platform={pn} atoms={top.getNumAtoms()}")
    sim.context.setPositions(pos)
    if args.resume_state:
        from openmm import XmlSerializer
        with open(args.resume_state) as fh:
            sim.context.setState(XmlSerializer.deserialize(fh.read()))
        log(f"resumed from portable state {args.resume_state} (skipped minimize/equil)")
    else:
        sim.minimizeEnergy(maxIterations=2000)
        sim.context.setVelocitiesToTemperature(300 * unit.kelvin, args.seed)

    spp = int(round(1000.0 / args.timestep_fs)); rep = max(1, int(args.report_ps * spp))
    def xi_and_obs():
        st = sim.context.getState(getPositions=True); P = np.array(st.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
        cents = np.array([P[g].mean(0) for g in groups]); com = P[all_ca].mean(0)
        proj = (cents - com) @ n_axis
        xi = float((proj * z_target_nm).sum() / sumz2)
        dd = m2_min_dists(P * 10.0, guan, walk); neng = int(sum(x < 8.0 for x in dd))
        rad, radcv, plan = ring_geometry(cents * 10.0)
        return xi, proj.tolist(), dd, neng, plan, rad

    lams = np.linspace(args.lam_min, args.lam_max, args.nwin)
    if args.passes == 'auto':
        passes = [('fwd', lams)] + ([('rev', lams[::-1])] if args.both else [])
    elif args.passes == 'fwd':
        passes = [('fwd', lams)]
    elif args.passes == 'rev':
        passes = [('rev', lams[::-1])]
    else:
        passes = [('fwd', lams), ('rev', lams[::-1])]
    log(f"passes={[p for p, _ in passes]}")
    windows = []
    for pname, seq in passes:
        for wi, lam in enumerate(seq):
            sim.context.setParameter('lam', float(lam))
            sim.step(int(args.equil_ps * spp))
            xis, projs, m2s, obs = [], [], [], []
            for _ in range(int(args.sample_ps / args.report_ps)):
                sim.step(rep); xi, proj, dd, neng, plan, rad = xi_and_obs()
                xis.append(xi); projs.append(proj); m2s.append(neng); obs.append([plan, rad])
            windows.append(dict(pass_=pname, win=wi, lam=float(lam), k_eff=k_eff,
                xi_mean=float(np.mean(xis)), xi=xis, per_sub_axial_nm=np.mean(projs, 0).tolist(),
                n_engaged_mean=float(np.mean(m2s)), planarity=float(np.mean([o[0] for o in obs])),
                radius=float(np.mean([o[1] for o in obs]))))
            json.dump({'params': vars(args), 'k_eff': k_eff, 'windows': windows},
                      open(os.path.join(args.out, 'window_data.json'), 'w'))
            sim.saveCheckpoint(os.path.join(args.out, 'ckpt.chk'))
            log(f"[{pname}] win {wi:2d} lam={lam:+.2f} xi={windows[-1]['xi_mean']:+.2f} "
                f"n_eng={windows[-1]['n_engaged_mean']:.1f} plan={windows[-1]['planarity']:.2f} rad={windows[-1]['radius']:.1f}")
        sim.saveState(os.path.join(args.out, f'state_{pname}_end.xml'))
        log(f"saved portable end-state state_{pname}_end.xml (resume seed for the opposite pass)")

    try:
        fwd = [w for w in windows if w['pass_'] == 'fwd']
        allxi = np.concatenate([np.array(w['xi']) for w in fwd])
        edges = np.linspace(allxi.min() - .05, allxi.max() + .05, 61); ctrs = 0.5 * (edges[:-1] + edges[1:]); kT = 2.494
        hist = np.array([np.histogram(w['xi'], edges)[0] for w in fwd], float)
        N = hist.sum(1); centers = np.array([w['lam'] for w in fwd])
        bias = 0.5 * k_eff * (ctrs[None, :] - centers[:, None]) ** 2 / kT
        F = np.zeros(len(fwd)); num = hist.sum(0)
        for _ in range(500):
            denom = (N[:, None] * np.exp(F[:, None] - bias)).sum(0)
            P = np.where(denom > 0, num / denom, 0); P /= P.sum()
            Fn = -np.log((P[None, :] * np.exp(-bias)).sum(1) + 1e-300)
            if np.max(np.abs(Fn - F)) < 1e-6: F = Fn; break
            F = Fn
        G = -kT * np.log(P + 1e-300); G -= np.nanmin(G[np.isfinite(G)])
        json.dump({'xi': ctrs.tolist(), 'G_kcal': G.tolist()}, open(os.path.join(args.out, 'pmf.json'), 'w'))
        gP = G[np.argmin(np.abs(ctrs))]; gH = G[np.argmin(np.abs(ctrs - 1))]
        bar = float(np.nanmax(G[(ctrs > .1) & (ctrs < .9)]) - min(gP, gH))
        log(f"WHAM: G(P)={gP:.1f} G(H)={gH:.1f} kcal/mol barrier~{bar:.1f}")
    except Exception as e:
        log("WHAM failed (raw window_data.json saved):", repr(e))
    log("UMBRELLA DONE")

if __name__ == '__main__': main()
