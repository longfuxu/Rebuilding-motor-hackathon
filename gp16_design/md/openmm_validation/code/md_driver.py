#!/usr/bin/env python3
"""OpenMM MD driver for the gp16 ring physics-validation (A=apo, B=7JQQ helical, C=design).

Runs on a Colab A100. Modes:
  bench_implicit / bench_explicit : time a short run -> print ns/day, exit (Phase 0).
  implicit                        : GBSA-OBC2 minimize->heat->NVT production (Phase 1).
  explicit                        : TIP3P + PME, minimize->NVT->NPT->production (Phase 2, optional).

Incremental: DCD appended, StateData CSV appended, checkpoint every --ckpt_ps.
Lockfile-guarded (O_EXCL) so a double-launch self-terminates.
"""
import argparse, os, sys, time, math

def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

def load_and_prep(inp, pH=7.0):
    """PDBFixer: fill missing sidechain atoms (NOT missing loops), add H at pH."""
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    fixer = PDBFixer(filename=inp)
    fixer.findMissingResidues()
    fixer.missingResidues = {}          # do NOT model disordered/terminal missing loops
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()            # missing heavy sidechain atoms within present residues
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH)
    return fixer

def build_implicit(topology, timestep_fs, hmr, gb_cutoff_nm=0.0):
    from openmm import app, unit, LangevinMiddleIntegrator
    ff = app.ForceField('amber14-all.xml', 'implicit/obc2.xml')
    kw = dict(constraints=app.HBonds, soluteDielectric=1.0, solventDielectric=78.5)
    if gb_cutoff_nm and gb_cutoff_nm > 0:
        kw['nonbondedMethod'] = app.CutoffNonPeriodic
        kw['nonbondedCutoff'] = gb_cutoff_nm*unit.nanometer
    else:
        kw['nonbondedMethod'] = app.NoCutoff
    if hmr:
        kw['hydrogenMass'] = 1.5*unit.amu
    system = ff.createSystem(topology, **kw)
    integ = LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond,
                                     timestep_fs*unit.femtoseconds)
    return ff, system, integ

def build_explicit(fixer, timestep_fs, hmr):
    from openmm import app, unit, LangevinMiddleIntegrator, MonteCarloBarostat
    ff = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
    modeller = app.Modeller(fixer.topology, fixer.positions)
    log("solvating (TIP3P, 1.0 nm pad, 0.15 M, neutralize)...")
    modeller.addSolvent(ff, model='tip3p', padding=1.0*unit.nanometer,
                        ionicStrength=0.15*unit.molar, neutralize=True)
    log(f"solvated atoms: {modeller.topology.getNumAtoms()}")
    kw = dict(nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometer,
              constraints=app.HBonds)
    if hmr:
        kw['hydrogenMass'] = 1.5*unit.amu
    system = ff.createSystem(modeller.topology, **kw)
    system.addForce(MonteCarloBarostat(1*unit.bar, 300*unit.kelvin, 25))
    integ = LangevinMiddleIntegrator(300*unit.kelvin, 1.0/unit.picosecond,
                                     timestep_fs*unit.femtoseconds)
    return ff, system, integ, modeller

def make_sim(topology, system, integ, platform='auto'):
    from openmm import app, Platform
    avail = {Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())}
    if platform == 'auto':
        order = ['CUDA', 'OpenCL', 'CPU', 'Reference']
        platform = next(p for p in order if p in avail)
    plat = Platform.getPlatformByName(platform)
    for props in ({'Precision': 'mixed'}, {}):
        try:
            sim = app.Simulation(topology, system, integ, plat, props)
            log(f"using platform: {platform} props={props}")
            return sim
        except Exception:
            continue
    return app.Simulation(topology, system, integ, plat)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--system', required=True, choices=['A','B','C'])
    ap.add_argument('--input', required=True)
    ap.add_argument('--mode', required=True,
                    choices=['bench_implicit','bench_explicit','implicit','explicit'])
    ap.add_argument('--ns', type=float, default=10.0)
    ap.add_argument('--equil_ps', type=float, default=200.0)
    ap.add_argument('--report_ps', type=float, default=20.0)
    ap.add_argument('--ckpt_ps', type=float, default=500.0)
    ap.add_argument('--timestep_fs', type=float, default=4.0)
    ap.add_argument('--hmr', type=int, default=1)
    ap.add_argument('--out', required=True)
    ap.add_argument('--bench_steps', type=int, default=2000)
    ap.add_argument('--min_iters', type=int, default=5000)
    ap.add_argument('--gb_cutoff_nm', type=float, default=0.0)
    ap.add_argument('--platform', default='auto')
    args = ap.parse_args()

    from openmm import app, unit, Platform
    os.makedirs(args.out, exist_ok=True)
    lock = os.path.join(args.out, '.lock')
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        os.write(fd, str(os.getpid()).encode()); os.close(fd)
    except FileExistsError:
        log(f"lock exists ({lock}); another driver owns {args.out}. Exiting."); sys.exit(0)

    t0 = time.time()
    log(f"=== system {args.system} mode {args.mode} input {args.input} ===")
    log("platforms available:",
        [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())],
        "requested:", args.platform)

    fixer = load_and_prep(args.input)
    nat_prot = fixer.topology.getNumAtoms()
    log(f"prepped protein atoms (with H): {nat_prot}")

    implicit = args.mode in ('bench_implicit','implicit')
    if implicit:
        ff, system, integ = build_implicit(fixer.topology, args.timestep_fs, args.hmr,
                                           args.gb_cutoff_nm)
        top, pos = fixer.topology, fixer.positions
    else:
        ff, system, integ, modeller = build_explicit(fixer, args.timestep_fs, args.hmr)
        top, pos = modeller.topology, modeller.positions
    n_atoms = top.getNumAtoms()
    log(f"system atoms: {n_atoms}  dof-timestep {args.timestep_fs} fs  hmr={args.hmr}")

    sim = make_sim(top, system, integ, args.platform)
    sim.context.setPositions(pos)

    # write the exact protonated starting topology (defines residue indexing for analysis)
    start_pdb = os.path.join(args.out, f'{args.system}_start.pdb')
    with open(start_pdb, 'w') as fh:
        app.PDBFile.writeFile(top, pos, fh, keepIds=True)
    log(f"wrote {start_pdb}")

    log(f"minimizing (maxIterations={args.min_iters})...")
    sim.minimizeEnergy(maxIterations=args.min_iters)
    minpos = sim.context.getState(getPositions=True).getPositions()
    with open(os.path.join(args.out, f'{args.system}_min.pdb'), 'w') as fh:
        app.PDBFile.writeFile(top, minpos, fh, keepIds=True)

    ts = args.timestep_fs * unit.femtoseconds
    steps_per_ns = int(round(1.0/ (args.timestep_fs*1e-6)))   # fs->ns
    report_steps = max(1, int(round(args.report_ps/1000.0 * steps_per_ns)))
    ckpt_steps   = max(report_steps, int(round(args.ckpt_ps/1000.0 * steps_per_ns)))

    # ---------------- benchmark mode ----------------
    if args.mode.startswith('bench'):
        sim.context.setVelocitiesToTemperature(300*unit.kelvin)
        log(f"benchmark: warming 200 steps then timing {args.bench_steps} steps...")
        sim.step(200)
        st = time.time()
        sim.step(args.bench_steps)
        el = time.time() - st
        ns_done = args.bench_steps * args.timestep_fs * 1e-6
        nsday = ns_done / (el/86400.0)
        log(f"BENCH system={args.system} mode={args.mode} atoms={n_atoms} "
            f"steps={args.bench_steps} dt={args.timestep_fs}fs wall={el:.1f}s "
            f"-> {nsday:.1f} ns/day")
        print(f"RESULT_NSDAY {args.system} {args.mode} {n_atoms} {nsday:.2f}", flush=True)
        os.remove(lock); return

    # ---------------- equilibration ----------------
    sim.context.setVelocitiesToTemperature(300*unit.kelvin)
    equil_steps = int(round(args.equil_ps/1000.0 * steps_per_ns))
    log(f"equilibration: {args.equil_ps} ps ({equil_steps} steps) NVT@300K...")
    if equil_steps > 0:
        sim.reporters.append(app.StateDataReporter(
            os.path.join(args.out, f'{args.system}_equil.csv'), max(1,report_steps),
            step=True, time=True, potentialEnergy=True, temperature=True, speed=True))
        sim.step(equil_steps)
        sim.reporters.clear()

    # ---------------- production ----------------
    prod_steps = int(round(args.ns * steps_per_ns))
    dcd = os.path.join(args.out, f'{args.system}_prod.dcd')
    csv = os.path.join(args.out, f'{args.system}_prod.csv')
    chk = os.path.join(args.out, f'{args.system}_prod.chk')
    append = os.path.exists(dcd)
    sim.reporters.append(app.DCDReporter(dcd, report_steps, append=append))
    sim.reporters.append(app.StateDataReporter(
        csv, report_steps, step=True, time=True, potentialEnergy=True,
        kineticEnergy=True, temperature=True, speed=True, elapsedTime=True,
        append=append))
    sim.reporters.append(app.StateDataReporter(
        sys.stdout, ckpt_steps, step=True, time=True, temperature=True,
        speed=True, progress=True, remainingTime=True, totalSteps=prod_steps))
    sim.reporters.append(app.CheckpointReporter(chk, ckpt_steps))
    log(f"production: {args.ns} ns ({prod_steps} steps), frame every {args.report_ps} ps")
    sim.step(prod_steps)

    finpos = sim.context.getState(getPositions=True).getPositions()
    with open(os.path.join(args.out, f'{args.system}_final.pdb'), 'w') as fh:
        app.PDBFile.writeFile(top, finpos, fh, keepIds=True)
    log(f"DONE system {args.system} in {(time.time()-t0)/60:.1f} min")
    os.remove(lock)

if __name__ == '__main__':
    main()
