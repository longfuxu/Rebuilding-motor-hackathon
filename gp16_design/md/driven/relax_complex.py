#!/usr/bin/env python3
"""Relax the rigidly-threaded design+dsDNA complex into a USABLE SMD start structure.

Rigid threading (build_dna_complex.py) tops out near ~1 A clearance = steric overlap
(the design pore differs from 7JQQ's). This script does the standard "insert then relax":
build the implicit system, hold the protein ring with a CA positional restraint (so the
DNA is not ejected and the ring frame is preserved), and energy-minimize -> the local
protein-DNA overlaps are pushed apart into normal contacts. Reports the full-atom
clearance before/after and writes C_plus_dna_relaxed.pdb.

This addresses the CRITICAL audit finding (0.427 A overlap in the raw complex). The
relaxed structure is the S0 the explicit-solvent SMD should start from.
"""
import argparse, os, time
import numpy as np

def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

def full_atom_clearance(top, pos_A):
    prot = [a.index for a in top.atoms() if a.residue.chain.id == 'A'
            and a.element is not None and a.element.symbol != 'H']
    dna = [a.index for a in top.atoms() if a.residue.chain.id in ('F', 'G')
           and a.element is not None and a.element.symbol != 'H']
    P, D = pos_A[prot], pos_A[dna]
    mn = 1e9
    for i in range(0, len(P), 2000):
        mn = min(mn, float(np.linalg.norm(P[i:i+2000][:, None, :] - D[None, :, :], axis=2).min()))
    return mn

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', default='inputs/C_plus_dna.pdb')
    ap.add_argument('--out', default='inputs/C_plus_dna_relaxed.pdb')
    ap.add_argument('--kanchor', type=float, default=800.0)   # kJ/mol/nm^2 on protein CA
    ap.add_argument('--min_iters', type=int, default=8000)
    ap.add_argument('--settle_ps', type=float, default=20.0)  # short restrained NVT to settle contacts
    ap.add_argument('--platform', default='auto')
    args = ap.parse_args()

    from openmm import (app, unit, Platform, LangevinMiddleIntegrator, CustomExternalForce)
    t0 = time.time()
    pdb = app.PDBFile(args.input)
    ff = app.ForceField('amber14-all.xml', 'implicit/obc2.xml')
    mod = app.Modeller(pdb.topology, pdb.positions)
    mod.addHydrogens(ff)
    top, pos = mod.topology, mod.positions
    P0 = np.array(pos.value_in_unit(unit.angstrom))
    cl0 = full_atom_clearance(top, P0)
    log(f"loaded {top.getNumAtoms()} atoms; full-atom clearance BEFORE = {cl0:.3f} A")

    system = ff.createSystem(top, nonbondedMethod=app.CutoffNonPeriodic,
                             nonbondedCutoff=2.0 * unit.nanometer, constraints=app.HBonds,
                             hydrogenMass=1.5 * unit.amu)
    # hold protein CA (keeps ring frame + prevents DNA ejection during settle)
    ca = [a.index for a in top.atoms() if a.name == 'CA' and a.residue.chain.id == 'A']
    anchor = CustomExternalForce('0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
    anchor.addGlobalParameter('k', args.kanchor)
    for p in ('x0', 'y0', 'z0'): anchor.addPerParticleParameter(p)
    for i in ca:
        x, y, z = pos[i].value_in_unit(unit.nanometer)
        anchor.addParticle(int(i), [x, y, z])
    system.addForce(anchor)

    integ = LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond, 2 * unit.femtoseconds)
    avail = {Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())}
    pn = args.platform if args.platform != 'auto' else \
        next(p for p in ['CUDA', 'OpenCL', 'CPU', 'Reference'] if p in avail)
    try:
        sim = app.Simulation(top, system, integ, Platform.getPlatformByName(pn), {'Precision': 'mixed'})
    except Exception:
        sim = app.Simulation(top, system, integ, Platform.getPlatformByName(pn))
    sim.context.setPositions(pos)
    log(f"platform={pn}; minimizing (maxIters={args.min_iters}) to relieve the overlap...")
    sim.minimizeEnergy(maxIterations=args.min_iters)
    Pmin = sim.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    log(f"clearance after MINIMIZE = {full_atom_clearance(top, Pmin):.3f} A")

    if args.settle_ps > 0:
        sim.context.setVelocitiesToTemperature(300 * unit.kelvin, 1)
        sim.step(int(args.settle_ps * 500))
        log(f"settled {args.settle_ps} ps (restrained NVT)")

    finpos = sim.context.getState(getPositions=True).getPositions()
    Pf = np.array(finpos.value_in_unit(unit.angstrom))
    clf = full_atom_clearance(top, Pf)
    with open(args.out, 'w') as fh:
        app.PDBFile.writeFile(top, finpos, fh, keepIds=True)
    log(f"DONE clearance {cl0:.3f} -> {clf:.3f} A  wrote {args.out}  ({(time.time()-t0)/60:.1f} min)")
    print(f"RESULT_RELAX clearance_before={cl0:.3f} clearance_after={clf:.3f} usable={clf>=2.2}", flush=True)

if __name__ == '__main__':
    main()
