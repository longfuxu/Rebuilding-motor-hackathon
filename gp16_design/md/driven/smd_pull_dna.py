#!/usr/bin/env python3
"""Track A / steps S1-S2 of md/STEERED_TARGETED_MD_PLAN.md — steered MD = in-silico
optical tweezers: pull the dsDNA along the ring channel axis and read the force + work.

Reaction coordinate: axial separation of the dsDNA COM from the ring COM,
  s = n . (COM_dna - COM_ring),   n = channel axis (from the 5 subunit CA centroids).
A moving harmonic restraint U = 1/2 kpull (s - s0(t))^2 pulls at constant velocity
(s0 ramps s_init -> s_init + pull). The ring is held by a weak CA positional restraint
(kanchor) -- this both keeps the lab axis valid AND mirrors the optical-tweezers geometry
(the motor/capsid is held on a bead while DNA is pulled).

Readouts: pulling force = kpull*(s0 - s) [reported in pN], displacement, forward Jarzynski
work W = sum_i [U(x_i;s0_{i+1}) - U(x_i;s0_i)]; and M2 (R146->neighbour WalkerA) tracking.
Run N seeds -> <F(x)>, Jarzynski/Crooks dG; compare design vs native, stall vs ~57 pN.

Modes:
  implicit : GBSA-OBC2 (local proof-of-concept; implicit handles DNA electrostatics poorly
             -> a PROXY force number, flagged as such).
  explicit : TIP3P + PME + 0.15 M + Mg, staged minimize->NVT->NPT->pull (the real S1/S2,
             ~250k atoms -> GCP A100). DNA=amber14/DNA.OL15, RNA=amber14/RNA.OL3.

KJ_PER_MOL_NM_TO_PN = 1.66054  (1 kJ/mol/nm = 1.66054 pN)
"""
import argparse, os, time, json
import numpy as np

def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

COPY_LOS = [1, 353, 705, 1057, 1409]
COPY_LEN = 342
R146_OFF = 254
WALK_OFF = range(132, 140)
GUAN = ('NE', 'CZ', 'NH1', 'NH2')
PN_PER_KJMOLNM = 1.66054


def protein_ca_and_m2(top):
    """Return (list of 5 CA-index lists, all_ca, guan[5], walk[5]) for the design chain A."""
    ch = None
    for c in top.chains():
        if c.id == 'A' or ch is None:
            ch = c; break
    m = {}
    for r in ch.residues():
        try: m[int(r.id)] = r
        except (TypeError, ValueError): pass
    groups, guan, walk = [], [], []
    for lo in COPY_LOS:
        hi = lo + COPY_LEN - 1
        groups.append([a.index for rn in range(lo, hi + 1) if rn in m
                       for a in m[rn].atoms() if a.name == 'CA'])
        r = m.get(lo + R146_OFF)
        guan.append([a.index for a in r.atoms() if a.name in GUAN] if r else [])
        walk.append([a.index for off in WALK_OFF if (lo + off) in m
                     for a in m[lo + off].atoms()
                     if a.element is not None and a.element.symbol != 'H'])
    all_ca = [i for g in groups for i in g]
    return groups, all_ca, guan, walk


def dna_atom_indices(top):
    idx = []
    for ch in top.chains():
        if ch.id in ('F', 'G'):
            idx += [a.index for r in ch.residues() for a in r.atoms()]
    return idx


def build_system(top, mode, hmr, positions=None):
    from openmm import app, unit, MonteCarloBarostat
    if mode == 'implicit':
        ff = app.ForceField('amber14-all.xml', 'implicit/obc2.xml')
        kw = dict(nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=2.0 * unit.nanometer,
                  constraints=app.HBonds, soluteDielectric=1.0, solventDielectric=78.5)
        if hmr: kw['hydrogenMass'] = 1.5 * unit.amu
        return ff, ff.createSystem(top, **kw), top, positions
    else:  # explicit
        ff = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
        mod = app.Modeller(top, positions)
        mod.addSolvent(ff, model='tip3p', padding=1.5 * unit.nanometer,
                       ionicStrength=0.15 * unit.molar, neutralize=True)
        kw = dict(nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometer,
                  constraints=app.HBonds)
        if hmr: kw['hydrogenMass'] = 1.5 * unit.amu
        system = ff.createSystem(mod.topology, **kw)
        # NVT (no barostat): a MonteCarloBarostat rescales positions and FIGHTS the absolute-position
        # solute-settle restraints -> NaN. The box is already at solvation density; NVT is standard for
        # SMD pulls. (Both prior explicit NaNs traced to the barostat.)
        return ff, system, mod.topology, mod.positions


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', default='inputs/C_plus_dna.pdb')
    ap.add_argument('--out', required=True)
    ap.add_argument('--mode', choices=['implicit', 'explicit'], default='implicit')
    ap.add_argument('--pull_A', type=float, default=17.0)     # ~5 bp translocation
    ap.add_argument('--pull_ps', type=float, default=500.0)
    ap.add_argument('--kpull', type=float, default=4000.0)    # kJ/mol/nm^2 on the CV
    ap.add_argument('--kanchor', type=float, default=800.0)   # kJ/mol/nm^2 on protein CA
    ap.add_argument('--equil_ps', type=float, default=50.0)
    ap.add_argument('--report_ps', type=float, default=2.0)
    ap.add_argument('--timestep_fs', type=float, default=2.0)
    ap.add_argument('--hmr', type=int, default=1)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--minimize', type=int, default=1)
    ap.add_argument('--platform', default='auto')
    args = ap.parse_args()

    from openmm import (app, unit, Platform, LangevinMiddleIntegrator,
                        CustomCentroidBondForce, CustomExternalForce)
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    log(f"=== SMD pull DNA  out={args.out} mode={args.mode} seed={args.seed} "
        f"pull={args.pull_A}A/{args.pull_ps}ps kpull={args.kpull} ===")

    pdb = app.PDBFile(args.input)
    # protonate DNA (crystal DNA lacks H); protein already has H
    ff0 = app.ForceField('amber14-all.xml', 'implicit/obc2.xml')
    mod0 = app.Modeller(pdb.topology, pdb.positions)
    mod0.addHydrogens(ff0)
    top, positions = mod0.topology, mod0.positions

    groups, all_ca, guan, walk = protein_ca_and_m2(top)
    dna = dna_atom_indices(top)
    log(f"protein CA={len(all_ca)}  DNA atoms={len(dna)}  total atoms={top.getNumAtoms()}")

    ff, system, top, positions = build_system(top, args.mode, args.hmr, positions)
    natom = system.getNumParticles()
    log(f"system particles={natom} ({args.mode})")

    P0 = np.array(positions.value_in_unit(unit.angstrom))
    cents0 = np.array([P0[g].mean(0) for g in groups])
    ctr0 = cents0.mean(0); X0 = cents0 - ctr0
    _, _, Vt = np.linalg.svd(X0); n_axis = Vt[2]
    if n_axis[2] < 0: n_axis = -n_axis
    dna_com0 = P0[dna].mean(0)
    s_init = float((dna_com0 - ctr0) @ n_axis) / 10.0     # nm
    log(f"axis={np.round(n_axis,3)}  s_init={s_init*10:.2f} A  pull to {s_init*10+args.pull_A:.2f} A")

    # ---- ring hold: weak CA positional restraint (holds lab frame; OT motor-on-bead) ----
    anchor = CustomExternalForce('0.5*kanchor*((x-ax)^2+(y-ay)^2+(z-az)^2)')
    anchor.addGlobalParameter('kanchor', args.kanchor)
    for p in ('ax', 'ay', 'az'): anchor.addPerParticleParameter(p)
    for i in all_ca:
        x, y, z = positions[i].value_in_unit(unit.nanometer)
        anchor.addParticle(int(i), [x, y, z])
    system.addForce(anchor)

    # ---- pull CV: axial separation of DNA COM from ring COM ----
    pull = CustomCentroidBondForce(
        2, "0.5*kpull*(nx*(x1-x2)+ny*(y1-y2)+nz*(z1-z2) - s0)^2")
    pull.addGlobalParameter('kpull', args.kpull)
    pull.addGlobalParameter('s0', s_init)
    for nm, v in [('nx', n_axis[0]), ('ny', n_axis[1]), ('nz', n_axis[2])]:
        pull.addGlobalParameter(nm, float(v))
    # equal weights -> geometric centroid, so the restraint matches the Python COM readout
    # (P[dna].mean(0)); default mass-weighting would restrain a DIFFERENT point than we measure.
    g_dna = pull.addGroup([int(i) for i in dna], [1.0] * len(dna))
    g_ring = pull.addGroup([int(i) for i in all_ca], [1.0] * len(all_ca))
    pull.addBond([g_dna, g_ring], [])
    system.addForce(pull)

    # solvent-settling restraint (explicit): hold ALL solute heavy atoms to their initial positions
    # while fresh water/ions relax around them, then release. Prevents the NaN from stepping a freshly
    # solvated 1.16M-atom box straight to 300 K dynamics.
    WATER_IONS = {'HOH', 'WAT', 'NA', 'CL', 'K', 'MG', 'NA+', 'CL-'}
    settle = CustomExternalForce('0.5*ksettle*((x-sx)^2+(y-sy)^2+(z-sz)^2)')
    settle.addGlobalParameter('ksettle', 0.0)
    for p in ('sx', 'sy', 'sz'):
        settle.addPerParticleParameter(p)
    if args.mode == 'explicit':
        solute_heavy = [a.index for a in top.atoms()
                        if a.residue.name not in WATER_IONS
                        and a.element is not None and a.element.symbol != 'H']
        for i in solute_heavy:
            x, y, z = positions[i].value_in_unit(unit.nanometer)
            settle.addParticle(int(i), [x, y, z])
        system.addForce(settle)
        log(f"settle restraint on {len(solute_heavy)} solute heavy atoms (off until after minimize)")

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
    log(f"platform={pn}")
    sim.context.setPositions(positions)

    spp = int(round(1000.0 / args.timestep_fs))
    if args.minimize:
        log("minimizing...")
        sim.minimizeEnergy(maxIterations=5000 if args.mode == 'explicit' else 3000)
    # staged solvent-settling for explicit: solute held, gradual heat 50->300 K, then release
    if args.mode == 'explicit':
        log("solvent-settling: solute restrained, heating 50->300 K, then 50 ps NPT, then release...")
        sim.context.setParameter('ksettle', 1000.0)
        sim.context.setVelocitiesToTemperature(50 * unit.kelvin, args.seed)
        for T in (50, 100, 150, 200, 250, 300):
            integ.setTemperature(T * unit.kelvin)
            sim.step(int(5 * spp))          # 5 ps per temperature rung
        sim.step(int(50 * spp))             # 50 ps at 300 K, solute still held (water/box settle)
        sim.context.setParameter('ksettle', 0.0)   # release the solute
        log("solute released; proceeding to free equilibration")
    sim.context.setVelocitiesToTemperature(300 * unit.kelvin, args.seed)
    if args.equil_ps > 0:
        log(f"equilibrating {args.equil_ps} ps (s0 fixed)...")
        sim.step(int(args.equil_ps * spp))

    # reset s0 to the equilibrated value so the pull starts from rest
    st = sim.context.getState(getPositions=True)
    P = st.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    s_start = float((P[dna].mean(0) - P[all_ca].mean(0)) @ n_axis) / 10.0
    sim.context.setParameter('s0', s_start)
    s_end = s_start + args.pull_A / 10.0

    pull_steps = int(args.pull_ps * spp)
    rep = max(1, int(args.report_ps * spp))
    nrep = pull_steps // rep
    sim.reporters.append(app.DCDReporter(os.path.join(args.out, 'traj.dcd'), rep))
    rows, W, s0 = [], 0.0, s_start
    log(f"pull {args.pull_ps} ps, s0 {s_start*10:.2f}->{s_end*10:.2f} A, {nrep} reports")
    for i in range(nrep):
        s0_next = s_start + (s_end - s_start) * ((i + 1) / nrep)
        sim.step(rep)
        st = sim.context.getState(getPositions=True)
        P = st.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        s = float((P[dna].mean(0) - P[all_ca].mean(0)) @ n_axis) / 10.0   # nm
        # force under the restraint center that ACTUALLY acted during this step (s0), not s0_next.
        # NOTE: single-snapshot force has ~sqrt(kBT*kpull) thermal noise (~150+ pN at kpull=4000);
        # a meaningful mean needs block-averaging over many samples/seeds -- see RESULTS caveats.
        force_pn = args.kpull * (s0 - s) * PN_PER_KJMOLNM                  # pN along axis
        dW = 0.5 * args.kpull * ((s - s0_next) ** 2 - (s - s0) ** 2)
        W += dW
        dmins = []
        for k in range(5):
            g, w = guan[k], walk[(k + 1) % 5]
            if g and w:
                dd = np.linalg.norm(P[np.array(g)][:, None] - P[np.array(w)][None], axis=2)
                dmins.append(float(dd.min()))
            else:
                dmins.append(float('nan'))
        n_eng = int(np.sum([d < 8.0 for d in dmins if d == d]))
        rows.append(dict(time_ps=(i + 1) * rep * args.timestep_fs / 1000.0,
                         disp_A=(s - s_start) * 10, s0_A=s0_next * 10, s_A=s * 10,
                         force_pN=force_pn, W_kcal=W / 4.184, n_engaged=n_eng,
                         dmin_A=[round(d, 2) for d in dmins]))
        sim.context.setParameter('s0', s0_next); s0 = s0_next
        if (i + 1) % max(1, nrep // 10) == 0:
            log(f"  {rows[-1]['time_ps']:.0f}ps disp={rows[-1]['disp_A']:.2f}A "
                f"F={force_pn:.0f}pN n_eng={n_eng} W={W/4.184:.1f}kcal/mol")

    with open(os.path.join(args.out, 'final.pdb'), 'w') as fh:
        app.PDBFile.writeFile(top, sim.context.getState(getPositions=True).getPositions(),
                              fh, keepIds=True)
    forces = np.array([r['force_pN'] for r in rows])
    res = dict(out=args.out, mode=args.mode, seed=args.seed, pull_A=args.pull_A,
               kpull=args.kpull, kanchor=args.kanchor,
               mean_force_pN=float(forces.mean()), max_force_pN=float(forces.max()),
               plateau_force_pN=float(forces[len(forces)//2:].mean()),
               W_total_kcal=W / 4.184, n_eng_start=rows[0]['n_engaged'],
               n_eng_final=rows[-1]['n_engaged'], n_eng_min=int(min(r['n_engaged'] for r in rows)),
               wall_min=(time.time() - t0) / 60.0,
               caveat=('implicit GBSA underestimates DNA electrostatics -> PROXY force; '
                       'use explicit+PME on GCP for the real number' if args.mode == 'implicit'
                       else 'explicit TIP3P+PME'))
    json.dump(dict(result=res, rows=rows), open(os.path.join(args.out, 'series.json'), 'w'), indent=2)
    with open(os.path.join(args.out, 'series.csv'), 'w') as fh:
        fh.write('time_ps,disp_A,s0_A,s_A,force_pN,W_kcal,n_engaged\n')
        for r in rows:
            fh.write(f"{r['time_ps']:.2f},{r['disp_A']:.3f},{r['s0_A']:.3f},{r['s_A']:.3f},"
                     f"{r['force_pN']:.2f},{r['W_kcal']:.3f},{r['n_engaged']}\n")
    log(f"DONE  <F>={forces.mean():.0f}pN plateau={res['plateau_force_pN']:.0f}pN "
        f"maxF={forces.max():.0f}pN W={W/4.184:.1f}kcal/mol  n_eng {rows[0]['n_engaged']}->"
        f"{rows[-1]['n_engaged']}  ({(time.time()-t0)/60:.1f} min)")
    print(f"RESULT_SMD mode={args.mode} seed={args.seed} meanF_pN={forces.mean():.1f} "
          f"plateauF_pN={res['plateau_force_pN']:.1f} W_kcal={W/4.184:.2f}", flush=True)


if __name__ == '__main__':
    main()
