#!/usr/bin/env python3
"""INTEGRATED helical<->planar + grip model — the closest thing here to the field's φ29 mechanism, and the
honest next step after the separate umbrella / ratchet runs.

  * the WHOLE ring cycles P->H->P via the collective staircase CV (subunit k axial target = lam_k(t)*z_staircase[k];
    lam triangles 0->1->0 per cycle = form the helical lock-washer, then collapse to planar) — this IS the
    helical<->planar model (NOT the bobbing hand-over-hand of tmd_ratchet);
  * a GRIP engages the dsDNA only during the H->P DESCENT (the power stroke) and releases on the P->H ascent,
    catch-and-carry, and it is anchored on the REAL DNA-contacting ring residues (<6 A from DNA), not on whole
    subunits — addressing reviewer concern (iv, real grip residues);
  * --mode concerted (all subunits descend together, symmetric) vs sequential (staggered descent);
  * multiple cycles + multiple seeds (reviewer concern ii, replicates).

Question it answers: does the helical<->planar cycle, WITH a physical descent-phase grip, translocate DNA — and
does it need asymmetry (sequential) or does the symmetric concerted collapse suffice? Output: traj.dcd +
series.json (per frame: cycle phase, per-subunit z, DNA axial pos, n_engaged, planarity).

STILL a driven/implicit model (concerns i,iii,iv-solvent remain — see PUBLICATION_READINESS.md); this run
removes the whole-subunit-grip and no-replicate objections and puts the mechanism in the right (helical<->planar)
frame. Run on A100.
"""
import argparse, os, time, json
import numpy as np
from tmd_staircase import log, COPY_LOS, COPY_LEN, HELIX_AXIAL_SPAN, ca_groups, m2_indices, ring_geometry


def dna_contact_residues(top, P0, dna_idx, groups_res, cutoff=6.0):
    """Per subunit, the ring heavy-atom indices whose residue has any atom < cutoff of DNA (real grip anchors)."""
    dnaP = P0[np.array(dna_idx)]
    anchors = []
    for blk in groups_res:      # blk = list of (resid, [atom indices]) for that subunit
        idxs = []
        for rid, aidx in blk:
            d = np.sqrt(((P0[np.array(aidx)][:, None, :] - dnaP[None, :, :]) ** 2).sum(-1)).min()
            if d < cutoff: idxs += aidx
        anchors.append(idxs)
    return anchors


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', default='inputs/C_plus_dna_relaxed.pdb')
    ap.add_argument('--out', required=True)
    ap.add_argument('--mode', choices=['concerted', 'sequential'], default='concerted')
    ap.add_argument('--ncycles', type=int, default=4)
    ap.add_argument('--cycle_ps', type=float, default=400.0)
    ap.add_argument('--kcv', type=float, default=30000.0)
    ap.add_argument('--kgrip', type=float, default=12000.0)
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
    # per-subunit residue->atoms (heavy) for grip-anchor detection
    ch = list(top.chains())[0]; resmap = {}
    for r in ch.residues():
        try: rid = int(r.id)
        except (TypeError, ValueError): continue
        resmap[rid] = [a.index for a in r.atoms() if a.element is not None and a.element.symbol != 'H']
    groups_res = []
    for lo in COPY_LOS:
        groups_res.append([(rn, resmap[rn]) for rn in range(lo, lo + COPY_LEN) if rn in resmap])
    anchors = dna_contact_residues(top, P0, dna_idx, groups_res)
    log(f"grip anchors / subunit (real DNA-contact residues): {[len(a) for a in anchors]} atoms")

    cents0 = np.array([P0[g].mean(0) for g in groups]); ctr0 = cents0.mean(0)
    _, _, Vt = np.linalg.svd(cents0 - ctr0); n_axis = Vt[2]
    if n_axis[2] < 0: n_axis = -n_axis
    amp = HELIX_AXIAL_SPAN / 2.0
    z_stair = (np.linspace(-amp, amp, 5)) / 10.0     # helical staircase target (nm), per subunit
    log(f"mode={args.mode} axis={np.round(n_axis,3)} z_stair_nm={np.round(z_stair,3)} x{args.ncycles} cyc")

    ff = app.ForceField('amber14-all.xml', 'implicit/obc2.xml')
    system = ff.createSystem(top, constraints=app.HBonds, hydrogenMass=1.5 * unit.amu,
                             nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=args.gb_cutoff_nm * unit.nanometer)

    # collective staircase CV: subunit_k axial -> zt_k (global-ish; set per report)
    cf = CustomCentroidBondForce(2, "0.5*kcv*(nx*(x1-x2)+ny*(y1-y2)+nz*(z1-z2) - zt)^2")
    cf.addGlobalParameter('kcv', args.kcv)
    for nm, v in [('nx', n_axis[0]), ('ny', n_axis[1]), ('nz', n_axis[2])]: cf.addGlobalParameter(nm, float(v))
    cf.addPerBondParameter('zt')
    gidx = [cf.addGroup([int(i) for i in g]) for g in groups]; gcom = cf.addGroup([int(i) for i in all_ca])
    zt_b = [cf.addBond([gidx[k], gcom], [0.0]) for k in range(5)]
    system.addForce(cf)

    # grip: DNA-COM axial <-> subunit_k DNA-contact-anchor axial, ON during descent, catch-and-carry
    gf = CustomCentroidBondForce(2, "0.5*kg*(nx*(x1-x2)+ny*(y1-y2)+nz*(z1-z2) - doff)^2")
    for nm, v in [('nx', n_axis[0]), ('ny', n_axis[1]), ('nz', n_axis[2])]: gf.addGlobalParameter(nm, float(v))
    gf.addPerBondParameter('kg'); gf.addPerBondParameter('doff')
    gdna = gf.addGroup([int(i) for i in dna_idx])
    ganch = [gf.addGroup([int(i) for i in (anchors[k] if anchors[k] else groups[k])]) for k in range(5)]
    gb = [gf.addBond([gdna, ganch[k]], [0.0, 0.0]) for k in range(5)]
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

    def lam_of(phase):    # triangle 0->1->0: 0 at planar, 1 at helical peak, back to 0
        return 1.0 - abs(1.0 - 2.0 * (phase % 1.0))
    def descending(phase):    # H->P power stroke = second half of the cycle
        return (phase % 1.0) >= 0.5

    rep = max(1, int(args.report_ps * spp)); nrep = int(args.ncycles * args.cycle_ps / args.report_ps)
    dcd = app.DCDReporter(os.path.join(args.out, 'traj.dcd'), rep); sim.reporters.append(dcd)
    rows = []; was_on = [False] * 5; doff = [0.0] * 5
    log(f"integrated {args.ncycles} cyc, {nrep} reports, mode={args.mode}")
    for i in range(nrep):
        st = sim.context.getState(getPositions=True); P = np.array(st.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
        cents = np.array([P[g].mean(0) for g in groups]); com = P[all_ca].mean(0)
        proj = (cents - com) @ n_axis
        anchz = [float((P[np.array(anchors[k] if anchors[k] else groups[k])].mean(0) - com) @ n_axis) for k in range(5)]
        dna_z_nm = float((P[dna_idx].mean(0) - com) @ n_axis)
        tcyc = i * args.report_ps / args.cycle_ps
        for k in range(5):
            ph = tcyc - (k / 5.0 if args.mode == 'sequential' else 0.0)
            cf.setBondParameters(zt_b[k], [gidx[k], gcom], [float(lam_of(ph) * z_stair[k])])
            on = descending(ph)
            if on and not was_on[k]:
                doff[k] = dna_z_nm - anchz[k]            # catch DNA at the current anchor offset
            gf.setBondParameters(gb[k], [gdna, ganch[k]], [args.kgrip if on else 0.0, doff[k]])
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
    dnet = rows[-1]['dna_z_A'] - rows[0]['dna_z_A']
    res = dict(mode=args.mode, model='helical<->planar + descent-grip', ncycles=args.ncycles,
               dna_net_A=round(dnet, 2), n_eng_start=rows[0]['n_engaged'], n_eng_final=rows[-1]['n_engaged'],
               n_eng_min=min(r['n_engaged'] for r in rows), seed=args.seed, wall_min=round((time.time() - t0) / 60, 1))
    json.dump(dict(result=res, rows=rows), open(os.path.join(args.out, 'series.json'), 'w'), indent=1)
    log(f"INTEGRATED_DONE mode={args.mode} DNA net {dnet:+.1f}A n_eng {res['n_eng_start']}->{res['n_eng_final']} (min {res['n_eng_min']})")


if __name__ == '__main__':
    main()
