#!/usr/bin/env python3
"""Track B / step T1 (fast, GPU-native) — targeted MD driving the design gp16 ring
planar -> helical along an AXIAL STAIRCASE reaction coordinate.

Why this CV (not RMSD): the endpoint geometry probe showed P->H is an axial staircase
(planarity 0.05->1.80 A, subunit-centroid axial span 0.1->4.8 A) with only ~2 A radial
change. So the essential coordinate is the per-subunit axial displacement. RMSDForce
works but on OpenCL it syncs GPU<->CPU every step (~40x slowdown); CustomCentroidBondForce
is GPU-native and translation-invariant here, so multi-seed Jarzynski is tractable on the Mac.

Drive: for each subunit k, restrain n.(centroid_k - ring_COM) toward lam*z_k, with lam
ramped 0->1 (steered).  n = ring axis; z_k = monotonic staircase [-a..+a], a = 7JQQ
axial-span/2, so the seam falls at the copy5->copy1 NON-covalent interface. Forward
Jarzynski work W = sum_i sum_k [U_k(x_i; lam_{i+1}) - U_k(x_i; lam_i)] accumulated
analytically. Implicit GBSA-OBC2, ~28k atoms.

Readouts per frame: lam, mean|axial proj|, ring planarity, radius, radiusCV, accumulated
work, and per-interface R146->neighbour-WalkerA min distance + n_engaged (M2 coupling).
"""
import argparse, os, time, json
import numpy as np

def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

COPY_LOS = [1, 353, 705, 1057, 1409]
COPY_LEN = 342
R146_OFF = 254
WALK_OFF = range(132, 140)
GUAN = ('NE', 'CZ', 'NH1', 'NH2')
HELIX_AXIAL_SPAN = 4.76   # A, from B_7jqq_helical endpoint geometry


def _resmap(top):
    ch = list(top.chains())[0]
    m = {}
    for r in ch.residues():
        try: m[int(r.id)] = r
        except (TypeError, ValueError): pass
    return m


def ca_groups(top):
    m = _resmap(top); groups = []
    for lo in COPY_LOS:
        hi = lo + COPY_LEN - 1
        groups.append([a.index for rn in range(lo, hi + 1) if rn in m
                       for a in m[rn].atoms() if a.name == 'CA'])
    return groups, m


def m2_indices(m):
    guan, walk = [], []
    for lo in COPY_LOS:
        r = m.get(lo + R146_OFF)
        guan.append([a.index for a in r.atoms() if a.name in GUAN] if r else [])
        walk.append([a.index for off in WALK_OFF if (lo + off) in m
                     for a in m[lo + off].atoms()
                     if a.element is not None and a.element.symbol != 'H'])
    return guan, walk


def ring_geometry(cents):
    C = np.asarray(cents); ctr = C.mean(0); X = C - ctr
    _, _, Vt = np.linalg.svd(X); nrm = Vt[2]
    oop = X @ nrm; planarity = float(np.sqrt((oop ** 2).mean()))
    inpl = X - np.outer(oop, nrm); radii = np.linalg.norm(inpl, axis=1)
    return float(radii.mean()), float(radii.std() / radii.mean()), planarity


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', default='../openmm_validation/trajectories/C/C_start.pdb')
    ap.add_argument('--out', required=True)
    ap.add_argument('--kcv', type=float, default=30000.0)   # kJ/mol/nm^2 per subunit
    ap.add_argument('--pull_ps', type=float, default=200.0)
    ap.add_argument('--equil_ps', type=float, default=20.0)
    ap.add_argument('--timestep_fs', type=float, default=2.0)
    ap.add_argument('--hmr', type=int, default=1)
    ap.add_argument('--report_ps', type=float, default=2.0)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--minimize', type=int, default=1)
    ap.add_argument('--gb_cutoff_nm', type=float, default=2.0)   # GB cutoff (NoCutoff GB is O(N^2), very slow)
    ap.add_argument('--platform', default='auto')
    ap.add_argument('--lam0', type=float, default=0.0, help='start of the lambda ramp (0=planar, 1=helical). Reverse H->P: --lam0 1 --lam1 0')
    ap.add_argument('--lam1', type=float, default=1.0, help='end of the lambda ramp')
    ap.add_argument('--pos_from', default='', help='PDB (same atoms as --input) to take INITIAL POSITIONS from, '
                    'while the CV axis/z_target are still computed from --input. For a reverse H->P pull: '
                    '--input <planar>.pdb --pos_from <helical_endpoint>.pdb --lam0 1 --lam1 0')
    ap.add_argument('--triangle', type=int, default=0, help='if 1, ramp lam0->lam1->lam0 in ONE continuous run '
                    '(a complete P->H->P cycle; DNA position stays continuous so a net-per-cycle is meaningful)')
    args = ap.parse_args()

    from openmm import app, unit, Platform, LangevinMiddleIntegrator, CustomCentroidBondForce
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    log(f"=== TMD staircase planar->helical  out={args.out} seed={args.seed} "
        f"kcv={args.kcv} pull={args.pull_ps}ps ===")

    pdb = app.PDBFile(args.input)
    top, pos = pdb.topology, pdb.positions
    sim_pos = app.PDBFile(args.pos_from).positions if args.pos_from else pos  # CV geometry from --input; start coords from here
    groups, m = ca_groups(top)
    guan, walk = m2_indices(m)
    all_ca = [i for g in groups for i in g]
    P0 = np.array(pos.value_in_unit(unit.angstrom))
    cents0 = np.array([P0[g].mean(0) for g in groups])
    ctr0 = cents0.mean(0); X0 = cents0 - ctr0
    _, _, Vt = np.linalg.svd(X0); n_axis = Vt[2]
    if n_axis[2] < 0: n_axis = -n_axis
    # target staircase (A) -> nm; monotonic so seam = copy5->copy1 non-covalent interface
    amp = HELIX_AXIAL_SPAN / 2.0
    z_target_A = np.linspace(-amp, amp, 5)
    z_target_nm = z_target_A / 10.0
    proj0 = (cents0 - ctr0) @ n_axis
    log(f"axis={np.round(n_axis,3)}  start axial proj (A)={np.round(proj0,2)}  "
        f"target staircase (A)={np.round(z_target_A,2)}")

    ff = app.ForceField('amber14-all.xml', 'implicit/obc2.xml')
    kw = dict(constraints=app.HBonds, soluteDielectric=1.0, solventDielectric=78.5)
    if args.gb_cutoff_nm and args.gb_cutoff_nm > 0:
        kw['nonbondedMethod'] = app.CutoffNonPeriodic
        kw['nonbondedCutoff'] = args.gb_cutoff_nm * unit.nanometer
    else:
        kw['nonbondedMethod'] = app.NoCutoff
    if args.hmr: kw['hydrogenMass'] = 1.5 * unit.amu
    system = ff.createSystem(top, **kw)

    # CustomCentroidBondForce: bond k connects (subunit_k, ring_COM); drive axial proj.
    cf = CustomCentroidBondForce(
        2, "0.5*kcv*(nx*(x1-x2)+ny*(y1-y2)+nz*(z1-z2) - lam*zt)^2")
    cf.addGlobalParameter('kcv', args.kcv)
    cf.addGlobalParameter('lam', 0.0)
    for nm, v in [('nx', n_axis[0]), ('ny', n_axis[1]), ('nz', n_axis[2])]:
        cf.addGlobalParameter(nm, float(v))
    cf.addPerBondParameter('zt')
    gidx = [cf.addGroup([int(i) for i in g]) for g in groups]
    gcom = cf.addGroup([int(i) for i in all_ca])
    for k in range(5):
        cf.addBond([gidx[k], gcom], [float(z_target_nm[k])])
    system.addForce(cf)

    integ = LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond,
                                     args.timestep_fs * unit.femtoseconds)
    integ.setRandomNumberSeed(args.seed)
    avail = {Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())}
    pn = args.platform if args.platform != 'auto' else \
        next(p for p in ['CUDA', 'OpenCL', 'CPU', 'Reference'] if p in avail)
    try:
        sim = app.Simulation(top, system, integ, Platform.getPlatformByName(pn),
                             {'Precision': 'mixed'})
    except Exception:
        sim = app.Simulation(top, system, integ, Platform.getPlatformByName(pn))
    log(f"platform={pn}  atoms={top.getNumAtoms()}")
    sim.context.setPositions(sim_pos)
    sim.context.setParameter('lam', args.lam0)   # set CV target to the ramp start BEFORE minimize (lam0=1 for reverse H->P)
    if args.minimize:
        log(f"minimizing (lam={args.lam0})...")
        sim.minimizeEnergy(maxIterations=2000)
    sim.context.setVelocitiesToTemperature(300 * unit.kelvin, args.seed)
    spp = int(round(1000.0 / args.timestep_fs))
    if args.equil_ps > 0:
        log(f"equilibrating {args.equil_ps} ps at lam={args.lam0}...")
        sim.step(int(args.equil_ps * spp))

    pull_steps = int(args.pull_ps * spp)
    rep = max(1, int(args.report_ps * spp))
    nrep = pull_steps // rep
    dcd = app.DCDReporter(os.path.join(args.out, 'traj.dcd'), rep)
    sim.reporters.append(dcd)
    zt = z_target_nm
    rows, W, lam = [], 0.0, args.lam0
    log(f"steered pull {args.pull_ps} ps, lam {args.lam0}->{args.lam1}{'->'+str(args.lam0) if args.triangle else ''}, {nrep} reports")
    for i in range(nrep):
        p = (i + 1) / nrep
        if args.triangle:   # lam0 -> lam1 (first half) -> lam0 (second half): one continuous P->H->P
            lam_next = args.lam0 + (args.lam1 - args.lam0) * (2 * p if p <= 0.5 else 2 * (1 - p))
        else:
            lam_next = args.lam0 + (args.lam1 - args.lam0) * p
        sim.step(rep)
        st = sim.context.getState(getPositions=True)
        P = st.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        cents = np.array([P[g].mean(0) for g in groups])
        com = P[all_ca].mean(0)
        proj_nm = ((cents - com) @ n_axis) / 10.0    # nm
        dW = float(np.sum(0.5 * args.kcv * ((proj_nm - lam_next * zt) ** 2
                                            - (proj_nm - lam * zt) ** 2)))
        W += dW
        radius, radcv, plan = ring_geometry(cents)
        dmins = []
        for k in range(5):
            g, w = guan[k], walk[(k + 1) % 5]
            if g and w:
                dd = np.linalg.norm(P[np.array(g)][:, None] - P[np.array(w)][None], axis=2)
                dmins.append(float(dd.min()))
            else:
                dmins.append(float('nan'))
        n_eng = int(np.sum([d < 8.0 for d in dmins if d == d]))
        proj_A = (cents - com) @ n_axis
        axspan = float(proj_A.max() - proj_A.min())
        rows.append(dict(time_ps=(i + 1) * rep * args.timestep_fs / 1000.0, lam=lam_next,
                         W_kcal=W / 4.184, planarity_A=plan, radius_A=radius, radiusCV=radcv,
                         axial_span_A=axspan, n_engaged=n_eng,
                         dmin_A=[round(d, 2) for d in dmins]))
        sim.context.setParameter('lam', lam_next); lam = lam_next
        if (i + 1) % max(1, nrep // 10) == 0:
            log(f"  {rows[-1]['time_ps']:.0f}ps lam={lam_next:.2f} plan={plan:.2f}A "
                f"axspan={axspan:.2f}A n_eng={n_eng} W={W/4.184:.1f}kcal/mol")

    with open(os.path.join(args.out, 'final.pdb'), 'w') as fh:
        app.PDBFile.writeFile(top, sim.context.getState(getPositions=True).getPositions(),
                              fh, keepIds=True)
    res = dict(out=args.out, seed=args.seed, kcv=args.kcv, pull_ps=args.pull_ps,
               planarity_start=rows[0]['planarity_A'], planarity_final=rows[-1]['planarity_A'],
               axial_span_start=rows[0]['axial_span_A'], axial_span_final=rows[-1]['axial_span_A'],
               radius_start=rows[0]['radius_A'], radius_final=rows[-1]['radius_A'],
               n_eng_start=rows[0]['n_engaged'], n_eng_final=rows[-1]['n_engaged'],
               n_eng_min=int(min(r['n_engaged'] for r in rows)), W_total_kcal=W / 4.184,
               target_staircase_A=[round(float(z), 2) for z in z_target_A],
               wall_min=(time.time() - t0) / 60.0)
    json.dump(dict(result=res, rows=rows), open(os.path.join(args.out, 'series.json'), 'w'), indent=2)
    with open(os.path.join(args.out, 'series.csv'), 'w') as fh:
        fh.write('time_ps,lam,W_kcal,planarity_A,radius_A,radiusCV,axial_span_A,n_engaged\n')
        for r in rows:
            fh.write(f"{r['time_ps']:.2f},{r['lam']:.3f},{r['W_kcal']:.3f},{r['planarity_A']:.3f},"
                     f"{r['radius_A']:.3f},{r['radiusCV']:.4f},{r['axial_span_A']:.3f},{r['n_engaged']}\n")
    log(f"DONE W={W/4.184:.1f} kcal/mol  planarity {rows[0]['planarity_A']:.2f}->{rows[-1]['planarity_A']:.2f}A "
        f"axspan {rows[0]['axial_span_A']:.2f}->{rows[-1]['axial_span_A']:.2f}A  "
        f"n_eng {rows[0]['n_engaged']}->{rows[-1]['n_engaged']}  ({(time.time()-t0)/60:.1f} min)")
    print(f"RESULT_TMD seed={args.seed} W_kcal={W/4.184:.2f} plan_final={rows[-1]['planarity_A']:.2f} "
          f"axspan_final={rows[-1]['axial_span_A']:.2f} n_eng {rows[0]['n_engaged']}->{rows[-1]['n_engaged']}", flush=True)


if __name__ == '__main__':
    main()
