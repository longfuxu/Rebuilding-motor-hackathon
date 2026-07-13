#!/usr/bin/env python3
"""C — minimal MECHANOCHEMICAL RATCHET model: put the ATP-driven power stroke + a timed DNA GRIP into the
single-chain ring, and SEE the DNA translocate. Unlike the plain staircase CV (which drives all subunits
together and lets DNA slip), here:
  * each subunit k has its OWN axial power stroke z_k(t): down (power) then up (reset) each cycle;
  * a per-subunit GRIP spring couples the dsDNA axial position to subunit k ONLY during k's down-stroke
    (grip ON = carry DNA), and releases during the up-stroke (grip OFF = let go) — a directional ratchet;
  * a CLOCK sets the firing order:  --mode concerted  = all 5 fire together (DNA steps once/cycle)
                                    --mode sequential = fire 1→2→3→4→5 in turn (DNA steps in substeps).
This is a MODEL (we impose grip+stroke — that is the point): it tests whether the ring's conformational
cycle CAN be coupled to processive DNA motion, and what the per-subunit z(t) + DNA-step choreography looks
like for each mechanism (the MINFLUX discriminator). It does NOT prove the real motor's timing.

Output: traj.dcd + series.json (per frame: cycle phase, per-subunit z, DNA axial pos, n_engaged).
Implicit GBSA, ring+dsDNA (~31k atoms). Run on A100 (OpenCL ok — CustomCentroidBondForce is native).
"""
import argparse, os, time, json
import numpy as np
from tmd_staircase import log, COPY_LOS, COPY_LEN, HELIX_AXIAL_SPAN, ca_groups, m2_indices, ring_geometry


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', default='inputs/C_plus_dna_relaxed.pdb')
    ap.add_argument('--out', required=True)
    ap.add_argument('--mode', choices=['concerted', 'sequential'], default='sequential')
    ap.add_argument('--ncycles', type=int, default=4)
    ap.add_argument('--cycle_ps', type=float, default=400.0)     # duration of one full P->H->P per subunit
    ap.add_argument('--step_A', type=float, default=3.4)         # DNA axial advance per cycle (1 bp rise)
    ap.add_argument('--kcv', type=float, default=30000.0)        # subunit axial restraint
    ap.add_argument('--kgrip', type=float, default=8000.0)       # DNA<->subunit grip spring (ON during down)
    ap.add_argument('--equil_ps', type=float, default=50.0)
    ap.add_argument('--report_ps', type=float, default=4.0)
    ap.add_argument('--timestep_fs', type=float, default=3.0)
    ap.add_argument('--gb_cutoff_nm', type=float, default=2.0)
    ap.add_argument('--platform', default='auto')
    ap.add_argument('--seed', type=int, default=1)
    args = ap.parse_args()
    from openmm import app, unit, Platform, LangevinMiddleIntegrator, CustomCentroidBondForce
    os.makedirs(args.out, exist_ok=True); t0 = time.time()
    json.dump(vars(args), open(os.path.join(args.out, 'params.json'), 'w'), indent=1)

    pdb = app.PDBFile(args.input); top, pos = pdb.topology, pdb.positions
    groups, m = ca_groups(top); guan, walk = m2_indices(m)
    all_ca = [i for g in groups for i in g]
    dna_idx = [a.index for a in top.atoms() if a.residue.chain.id in ('F', 'G')
               and a.element is not None and a.element.symbol != 'H']
    P0 = np.array(pos.value_in_unit(unit.angstrom))
    cents0 = np.array([P0[g].mean(0) for g in groups]); ctr0 = cents0.mean(0)
    _, _, Vt = np.linalg.svd(cents0 - ctr0); n_axis = Vt[2]
    if n_axis[2] < 0: n_axis = -n_axis
    z_base = np.zeros(5)        # subunits rest at planar (axial 0); the power stroke dips each one DOWN by --step_A
    log(f"mode={args.mode} axis={np.round(n_axis,3)} rest=planar step={args.step_A}A x {args.ncycles} cycles")

    ff = app.ForceField('amber14-all.xml', 'implicit/obc2.xml')
    system = ff.createSystem(top, constraints=app.HBonds, hydrogenMass=1.5 * unit.amu,
                             nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=args.gb_cutoff_nm * unit.nanometer)

    # per-subunit axial power-stroke restraint: drive subunit_k axial(rel ring COM) toward zt_k (global param)
    cf = CustomCentroidBondForce(2, "0.5*kcv*(nx*(x1-x2)+ny*(y1-y2)+nz*(z1-z2) - zt)^2")
    cf.addGlobalParameter('kcv', args.kcv)
    for nm, v in [('nx', n_axis[0]), ('ny', n_axis[1]), ('nz', n_axis[2])]: cf.addGlobalParameter(nm, float(v))
    cf.addPerBondParameter('zt')
    gidx = [cf.addGroup([int(i) for i in g]) for g in groups]; gcom = cf.addGroup([int(i) for i in all_ca])
    zt_params = []
    for k in range(5):
        b = cf.addBond([gidx[k], gcom], [float(z_base[k])]); zt_params.append(b)
    system.addForce(cf)

    # GRIP: spring DNA-COM axial toward (subunit_k axial + doff_k). doff_k is CAUGHT at the current offset the
    # moment subunit k engages (catch-and-carry), so each grip holds the DNA where it is and carries it down;
    # released on the up-stroke -> a directional hand-over-hand ratchet. gf needs its OWN groups.
    gf = CustomCentroidBondForce(2, "0.5*kg*(nx*(x1-x2)+ny*(y1-y2)+nz*(z1-z2) - doff)^2")
    gf.addGlobalParameter('nx', float(n_axis[0])); gf.addGlobalParameter('ny', float(n_axis[1])); gf.addGlobalParameter('nz', float(n_axis[2]))
    gf.addPerBondParameter('kg'); gf.addPerBondParameter('doff')
    gdna = gf.addGroup([int(i) for i in dna_idx])
    gf_sub = [gf.addGroup([int(i) for i in groups[k]]) for k in range(5)]
    for k in range(5):
        gf.addBond([gdna, gf_sub[k]], [0.0, 0.0])
    system.addForce(gf)

    integ = LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond, args.timestep_fs * unit.femtoseconds)
    integ.setRandomNumberSeed(args.seed)
    avail = {Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())}
    pn = args.platform if args.platform != 'auto' else next(p for p in ['CUDA', 'OpenCL', 'CPU', 'Reference'] if p in avail)
    try: sim = app.Simulation(top, system, integ, Platform.getPlatformByName(pn), {'Precision': 'mixed'})
    except Exception: sim = app.Simulation(top, system, integ, Platform.getPlatformByName(pn))
    log(f"platform={pn} atoms={top.getNumAtoms()}")
    sim.context.setPositions(pos); sim.minimizeEnergy(maxIterations=2000)
    sim.context.setVelocitiesToTemperature(300 * unit.kelvin, args.seed)
    spp = int(round(1000.0 / args.timestep_fs))
    if args.equil_ps > 0: sim.step(int(args.equil_ps * spp))

    def stroke(phase):   # 0->1 down (power), 1->0 up (reset); phase in [0,1)
        return 1.0 - abs(1.0 - 2.0 * (phase % 1.0))
    def gripped(phase):  # grip ON during the down (power) half
        return (phase % 1.0) < 0.5

    rep = max(1, int(args.report_ps * spp)); nrep = int(args.ncycles * args.cycle_ps / args.report_ps)
    dcd = app.DCDReporter(os.path.join(args.out, 'traj.dcd'), rep); sim.reporters.append(dcd)
    rows = []; was_on = [False] * 5; doff = [0.0] * 5
    log(f"ratchet {args.ncycles} cycles, {nrep} reports")
    for i in range(nrep):
        st = sim.context.getState(getPositions=True); P = np.array(st.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
        cents = np.array([P[g].mean(0) for g in groups]); com = P[all_ca].mean(0)
        proj = (cents - com) @ n_axis                                 # per-subunit axial (nm)
        dna_z_nm = float((P[dna_idx].mean(0) - com) @ n_axis)
        tcyc = i * args.report_ps / args.cycle_ps
        for k in range(5):
            ph = tcyc - (k / 5.0 if args.mode == 'sequential' else 0.0)     # sequential = phase-lag per subunit
            depth = stroke(ph)                                        # 0..1 into the power stroke
            cf.setBondParameters(zt_params[k], [gidx[k], gcom], [float(z_base[k] - depth * (args.step_A / 10.0))])
            on = gripped(ph)
            if on and not was_on[k]:
                doff[k] = dna_z_nm - float(proj[k])                   # CATCH: grip the DNA where it currently is
            gf.setBondParameters(k, [gdna, gf_sub[k]], [args.kgrip if on else 0.0, doff[k]])
            was_on[k] = on
        cf.updateParametersInContext(sim.context); gf.updateParametersInContext(sim.context)
        dmins = []
        for k in range(5):
            g, w = guan[k], walk[(k + 1) % 5]
            dmins.append(float(np.sqrt(((P[g][:, None] - P[w][None]) ** 2 * 100).sum(-1)).min()) if g and w else 99.0)
        neng = int(sum(d < 8.0 for d in dmins))
        rows.append(dict(t_cyc=round(tcyc, 3), dna_z_A=round(dna_z_nm * 10, 2),
                         per_sub_z_A=[round(float(x) * 10, 2) for x in proj], n_engaged=neng,
                         planarity_A=round(ring_geometry(cents * 10.0)[2], 2)))
        sim.step(rep)
        if (i + 1) % max(1, nrep // 10) == 0:
            log(f"  cyc {tcyc:.2f} dna_z {dna_z_nm*10:+.1f}A n_eng {neng} plan {rows[-1]['planarity_A']}")
    dna0 = rows[0]['dna_z_A']; dnaf = rows[-1]['dna_z_A']
    res = dict(mode=args.mode, ncycles=args.ncycles, step_A=args.step_A, dna_z_start=dna0, dna_z_final=dnaf,
               dna_net_A=round(dnaf - dna0, 2), n_eng_start=rows[0]['n_engaged'], n_eng_final=rows[-1]['n_engaged'],
               n_eng_min=min(r['n_engaged'] for r in rows), wall_min=round((time.time() - t0) / 60, 1))
    json.dump(dict(result=res, rows=rows), open(os.path.join(args.out, 'series.json'), 'w'), indent=1)
    log(f"RATCHET_DONE mode={args.mode} DNA net {res['dna_net_A']:+.1f}A (expect ~{args.step_A*args.ncycles:.1f} if it ratchets) "
        f"n_eng {res['n_eng_start']}->{res['n_eng_final']} (min {res['n_eng_min']})")


if __name__ == '__main__':
    main()
