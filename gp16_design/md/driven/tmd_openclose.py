#!/usr/bin/env python3
"""Track B / step T1 of md/STEERED_TARGETED_MD_PLAN.md — targeted MD: drive the design
gp16 ring planar -> helical and read the work + whether M2 coupling stays engaged.

The reaction coordinate is CA-RMSD to a HELICAL MORPH TARGET of the design's own atoms:
each of the 5 subunits is rigidly translated to impose the 7JQQ-amplitude staircase
(axial span ~4.76 A, radius +2.2 A) measured from B_7jqq_helical.pdb. A monotonic
staircase in copy order (copy1->..->copy5) puts the seam at the single chain's one
NON-covalent interface (copy5->copy1) -- the physically right place to open.

Method: OpenMM RMSDForce (optimal-aligned) wrapped in a CustomCVForce harmonic
restraint U = 1/2 k (rmsd - r0)^2. r0 is ramped RMSD0 -> ~0 at constant speed (steered/
targeted MD). Forward Jarzynski work W = sum_i [U(x_i;r0_{i+1}) - U(x_i;r0_i)], accumulated
analytically. Implicit GBSA-OBC2 (design-only ~28k atoms) -> runs on Mac OpenCL/CPU.

Readouts per frame: r0, actual RMSD, restraint energy, accumulated work, ring planarity,
ring radius, and per-interface R146->neighbour-WalkerA min distance + n_engaged (M2).
"""
import argparse, os, sys, time, json
import numpy as np

def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

# design copy layout (absolute resSeq), from analyze.py / score_m2.py
COPY_LOS = [1, 353, 705, 1057, 1409]
COPY_LEN = 342                      # copy span lo..lo+341
R146_OFF = 254                      # R146 at lo+254
WALK_OFF = range(132, 140)          # Walker-A at lo+132..lo+139
GUAN = ('NE', 'CZ', 'NH1', 'NH2')

# 7JQQ helical amplitude (measured from B_7jqq_helical.pdb, see endpoint geometry probe):
#   planar radius ~25.8 A (design), helical radius ~29.1 A  -> +2.2 A radial (use design+dR)
#   axial span 4.76 A -> monotonic staircase amplitude per subunit below
HELIX_AXIAL_SPAN = 4.76             # A (7JQQ)
HELIX_DRADIUS    = 2.2              # A radial expansion planar->helical (26.9->29.1)


def load_topology_positions(pdb):
    from openmm.app import PDBFile
    p = PDBFile(pdb)
    return p.topology, p.positions


def _resmap(top):
    """resSeq(int) -> OpenMM Residue, over the single design chain."""
    ch = list(top.chains())[0]
    resmap = {}
    for r in ch.residues():
        try:
            resmap[int(r.id)] = r
        except (TypeError, ValueError):
            pass
    return resmap


def subunit_ca_indices(top):
    """Return (list of 5 arrays of CA atom indices per copy, resmap)."""
    resmap = _resmap(top)
    ca = []
    for lo in COPY_LOS:
        hi = lo + COPY_LEN - 1
        idx = [a.index for rn in range(lo, hi + 1) if rn in resmap
               for a in resmap[rn].atoms() if a.name == 'CA']
        ca.append(np.array(idx, dtype=int))
    return ca, resmap


def subunit_all_indices(top):
    """Return list of 5 arrays of ALL atom indices per copy (for rigid translation)."""
    resmap = _resmap(top)
    allidx = []
    for lo in COPY_LOS:
        hi = lo + COPY_LEN - 1
        idx = [a.index for rn in range(lo, hi + 1) if rn in resmap
               for a in resmap[rn].atoms()]
        allidx.append(np.array(idx, dtype=int))
    return allidx


def m2_atom_indices(top, resmap):
    """Per copy: guanidinium atom indices of R146 and heavy Walker-A atoms."""
    guan, walk = [], []
    for lo in COPY_LOS:
        r146 = resmap.get(lo + R146_OFF)
        g = [a.index for a in r146.atoms() if a.name in GUAN] if r146 else []
        w = [a.index for off in WALK_OFF if (lo + off) in resmap
             for a in resmap[lo + off].atoms()
             if a.element is not None and a.element.symbol != 'H']
        guan.append(np.array(g, dtype=int)); walk.append(np.array(w, dtype=int))
    return guan, walk


def build_helical_target(pos_nm, ca_idx, all_idx):
    """Return target positions (nm array Nx3): each subunit rigidly translated so its
    centroid moves onto a 7JQQ-amplitude helical staircase (axial) + radial expansion."""
    P = np.asarray(pos_nm, dtype=float) * 10.0   # nm -> A
    cents = np.array([P[ca_idx[k]].mean(0) for k in range(5)])
    ctr = cents.mean(0)
    X = cents - ctr
    _, _, Vt = np.linalg.svd(X)
    axis = Vt[2]                                   # ring normal
    # in-plane radial unit vector per subunit
    oop = X @ axis
    inplane = X - np.outer(oop, axis)
    radii = np.linalg.norm(inplane, axis=1)
    radial_hat = inplane / radii[:, None]
    # target axial staircase (monotonic in copy order; seam at copy5->copy1 non-covalent iface)
    amp = HELIX_AXIAL_SPAN / 2.0
    zt = np.linspace(-amp, amp, 5)                 # A, [-2.38..+2.38]
    target = P.copy()
    disp = []
    for k in range(5):
        new_cent = ctr + radial_hat[k] * (radii[k] + HELIX_DRADIUS) + axis * zt[k]
        d = new_cent - cents[k]
        disp.append(np.linalg.norm(d))
        target[all_idx[k]] += d                    # rigid translation of whole subunit
    return target / 10.0, np.array(disp)           # back to nm, per-subunit disp (A)


def ring_geometry(cent):
    C = np.asarray(cent); ctr = C.mean(0); Xc = C - ctr
    _, _, Vt = np.linalg.svd(Xc); normal = Vt[2]
    oop = Xc @ normal; planarity = float(np.sqrt((oop ** 2).mean()))
    inplane = Xc - np.outer(oop, normal); radii = np.linalg.norm(inplane, axis=1)
    return float(radii.mean()), float(radii.std() / radii.mean()), planarity


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', default='../openmm_validation/trajectories/C/C_start.pdb')
    ap.add_argument('--out', required=True)
    ap.add_argument('--k', type=float, default=5000.0)         # kJ/mol/nm^2
    ap.add_argument('--pull_ps', type=float, default=300.0)
    ap.add_argument('--equil_ps', type=float, default=20.0)
    ap.add_argument('--timestep_fs', type=float, default=2.0)  # HBonds+HMR->4 ok; 2 safe for TMD
    ap.add_argument('--hmr', type=int, default=1)
    ap.add_argument('--report_ps', type=float, default=2.0)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--r0_floor_nm', type=float, default=0.0)
    ap.add_argument('--platform', default='auto')
    ap.add_argument('--minimize', type=int, default=1)
    ap.add_argument('--gb_cutoff_nm', type=float, default=0.0)   # >0 => CutoffNonPeriodic GB (fast)
    args = ap.parse_args()

    from openmm import app, unit, Platform, LangevinMiddleIntegrator
    from openmm import RMSDForce, CustomCVForce

    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    log(f"=== TMD planar->helical  out={args.out}  seed={args.seed}  k={args.k}  pull={args.pull_ps}ps ===")

    top, pos = load_topology_positions(args.input)
    ca_idx, resmap = subunit_ca_indices(top)
    all_idx = subunit_all_indices(top)
    guan, walk = m2_atom_indices(top, resmap)
    ca_all = np.concatenate(ca_idx)
    log(f"atoms={top.getNumAtoms()}  CA/subunit={[len(c) for c in ca_idx]}")

    pos_nm = np.array(pos.value_in_unit(unit.nanometer))
    target_nm, disp = build_helical_target(pos_nm, ca_idx, all_idx)
    log(f"per-subunit morph displacement (A): {['%.2f'%d for d in disp]}")

    # ---- build implicit system ----
    ff = app.ForceField('amber14-all.xml', 'implicit/obc2.xml')
    kw = dict(constraints=app.HBonds, soluteDielectric=1.0, solventDielectric=78.5)
    if args.gb_cutoff_nm and args.gb_cutoff_nm > 0:   # NoCutoff GB is O(N^2), very slow
        kw['nonbondedMethod'] = app.CutoffNonPeriodic
        kw['nonbondedCutoff'] = args.gb_cutoff_nm * unit.nanometer
    else:
        kw['nonbondedMethod'] = app.NoCutoff
    if args.hmr:
        kw['hydrogenMass'] = 1.5 * unit.amu
    system = ff.createSystem(top, **kw)

    # ---- RMSD restraint (optimal-aligned) over all CA ----
    rmsd = RMSDForce(target_nm, [int(i) for i in ca_all])
    cv = CustomCVForce('0.5*k*(rmsd - r0)^2')
    cv.addCollectiveVariable('rmsd', rmsd)
    cv.addGlobalParameter('k', args.k)
    cv.addGlobalParameter('r0', 1.0)   # set below
    system.addForce(cv)

    integ = LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond,
                                     args.timestep_fs * unit.femtoseconds)
    integ.setRandomNumberSeed(args.seed)

    avail = {Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())}
    plat_name = args.platform
    if plat_name == 'auto':
        plat_name = next(p for p in ['CUDA', 'OpenCL', 'CPU', 'Reference'] if p in avail)
    plat = Platform.getPlatformByName(plat_name)
    try:
        sim = app.Simulation(top, system, integ, plat, {'Precision': 'mixed'})
    except Exception:
        sim = app.Simulation(top, system, integ, plat)
    log(f"platform={plat_name}")
    sim.context.setPositions(pos)

    # measure initial RMSD to target (r0 must start there so restraint is satisfied)
    rmsd0 = cv.getCollectiveVariableValues(sim.context)[0]   # nm
    log(f"initial CA-RMSD to helical target = {rmsd0*10:.2f} A")
    sim.context.setParameter('r0', rmsd0)

    if args.minimize:
        log("minimizing (restraint on, r0=rmsd0)...")
        sim.minimizeEnergy(maxIterations=2000)

    sim.context.setVelocitiesToTemperature(300 * unit.kelvin, args.seed)
    steps_per_ps = int(round(1000.0 / args.timestep_fs))
    equil_steps = int(args.equil_ps * steps_per_ps)
    if equil_steps > 0:
        log(f"equilibrating {args.equil_ps} ps at r0=rmsd0...")
        sim.step(equil_steps)

    # ---- steered ramp r0: rmsd0 -> floor ----
    pull_steps = int(args.pull_ps * steps_per_ps)
    report_steps = max(1, int(args.report_ps * steps_per_ps))
    n_reports = pull_steps // report_steps
    r0_start = float(cv.getCollectiveVariableValues(sim.context)[0])
    sim.context.setParameter('r0', r0_start)   # sync applied restraint center to the pull's start
    r0_end = args.r0_floor_nm
    dcd = app.DCDReporter(os.path.join(args.out, 'traj.dcd'), report_steps)
    sim.reporters.append(dcd)

    rows = []
    W = 0.0                      # kJ/mol, forward Jarzynski work
    r0 = r0_start
    log(f"steered pull: {args.pull_ps} ps, r0 {r0_start*10:.2f} -> {r0_end*10:.2f} A, "
        f"{n_reports} reports")
    for i in range(n_reports):
        r0_next = r0_start + (r0_end - r0_start) * ((i + 1) / n_reports)
        # advance dynamics at current r0
        sim.step(report_steps)
        st = sim.context.getState(getEnergy=True, getPositions=True)
        cur_rmsd = float(cv.getCollectiveVariableValues(sim.context)[0])   # nm
        # discrete Jarzynski work increment: U(x; r0_next) - U(x; r0)
        dW = 0.5 * args.k * ((cur_rmsd - r0_next) ** 2 - (cur_rmsd - r0) ** 2)
        W += dW
        # geometry + M2 readouts
        P = st.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        cents = [P[ca_idx[k]].mean(0) for k in range(5)]
        radius, radcv, plan = ring_geometry(cents)
        dmins = []
        for k in range(5):
            acc = (k + 1) % 5
            g, w = guan[k], walk[acc]
            if len(g) and len(w):
                dd = np.linalg.norm(P[g][:, None, :] - P[w][None, :, :], axis=2)
                dmins.append(float(dd.min()))
            else:
                dmins.append(float('nan'))
        n_eng = int(np.sum([d < 8.0 for d in dmins if d == d]))
        Ucv = 0.5 * args.k * (cur_rmsd - r0_next) ** 2
        rows.append(dict(step=(i + 1) * report_steps,
                         time_ps=(i + 1) * report_steps * args.timestep_fs / 1000.0,
                         r0_A=r0_next * 10, rmsd_A=cur_rmsd * 10, U_kJ=Ucv,
                         W_kJ=W, W_kcal=W / 4.184, planarity_A=plan, radius_A=radius,
                         radiusCV=radcv, n_engaged=n_eng,
                         dmin_A=[round(d, 2) for d in dmins]))
        sim.context.setParameter('r0', r0_next)
        r0 = r0_next
        if (i + 1) % max(1, n_reports // 10) == 0:
            log(f"  {rows[-1]['time_ps']:.0f}ps  r0={r0_next*10:.2f}A rmsd={cur_rmsd*10:.2f}A "
                f"plan={plan:.2f}A n_eng={n_eng} W={W/4.184:.1f}kcal/mol")

    # final structure
    finpos = sim.context.getState(getPositions=True).getPositions()
    with open(os.path.join(args.out, 'final.pdb'), 'w') as fh:
        app.PDBFile.writeFile(top, finpos, fh, keepIds=True)

    result = dict(out=args.out, seed=args.seed, k=args.k, pull_ps=args.pull_ps,
                  rmsd0_A=rmsd0 * 10, rmsd_final_A=rows[-1]['rmsd_A'],
                  planarity_final_A=rows[-1]['planarity_A'],
                  radius_final_A=rows[-1]['radius_A'],
                  n_engaged_start=rows[0]['n_engaged'], n_engaged_final=rows[-1]['n_engaged'],
                  n_engaged_min=int(min(r['n_engaged'] for r in rows)),
                  W_total_kcal=W / 4.184, per_subunit_morph_disp_A=[round(float(d), 2) for d in disp],
                  wall_min=(time.time() - t0) / 60.0)
    with open(os.path.join(args.out, 'series.json'), 'w') as fh:
        json.dump(dict(result=result, rows=rows), fh, indent=2)
    with open(os.path.join(args.out, 'series.csv'), 'w') as fh:
        fh.write('time_ps,r0_A,rmsd_A,U_kJ,W_kcal,planarity_A,radius_A,radiusCV,n_engaged\n')
        for r in rows:
            fh.write(f"{r['time_ps']:.2f},{r['r0_A']:.3f},{r['rmsd_A']:.3f},{r['U_kJ']:.2f},"
                     f"{r['W_kcal']:.3f},{r['planarity_A']:.3f},{r['radius_A']:.3f},"
                     f"{r['radiusCV']:.4f},{r['n_engaged']}\n")
    log(f"DONE  W={W/4.184:.1f} kcal/mol  rmsd {rmsd0*10:.2f}->{rows[-1]['rmsd_A']:.2f} A  "
        f"planarity {rows[0]['planarity_A']:.2f}->{rows[-1]['planarity_A']:.2f} A  "
        f"n_eng {rows[0]['n_engaged']}->{rows[-1]['n_engaged']}  ({(time.time()-t0)/60:.1f} min)")
    print(f"RESULT_TMD seed={args.seed} W_kcal={W/4.184:.2f} "
          f"planarity_final={rows[-1]['planarity_A']:.2f} n_eng_final={rows[-1]['n_engaged']}", flush=True)


if __name__ == '__main__':
    main()
