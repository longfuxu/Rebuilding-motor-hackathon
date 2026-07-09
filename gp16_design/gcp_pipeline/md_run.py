import sys, json, time
import numpy as np
from openmm.app import *
from openmm import *
from openmm.unit import *
from pdbfixer import PDBFixer

INP = sys.argv[1]; TAG = sys.argv[2]; NSTEPS = int(sys.argv[3]) if len(sys.argv) > 3 else 150000
t0 = time.time()


def kabsch_rmsd(P, Q):
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    H = Pc.T @ Qc; U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T)); D = np.diag([1, 1, d]); R = Vt.T @ D @ U.T
    return float(np.sqrt(((Pc @ R.T - Qc) ** 2).sum(1).mean()))


fixer = PDBFixer(filename=INP)
fixer.findMissingResidues(); fixer.findMissingAtoms(); fixer.addMissingAtoms(); fixer.addMissingHydrogens(7.0)
ff = ForceField('amber14-all.xml', 'implicit/gbn2.xml')
modeller = Modeller(fixer.topology, fixer.positions)
system = ff.createSystem(modeller.topology, nonbondedMethod=CutoffNonPeriodic,
                         nonbondedCutoff=2.0 * nanometer, constraints=HBonds)
integ = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 0.002 * picoseconds)
try:
    plat = Platform.getPlatformByName('CUDA')
except Exception:
    plat = Platform.getPlatformByName('CPU')
sim = Simulation(modeller.topology, system, integ, plat)
sim.context.setPositions(modeller.positions)
sim.minimizeEnergy()
ca_idx = [a.index for a in modeller.topology.atoms() if a.name == 'CA']
start = np.array(sim.context.getState(getPositions=True).getPositions().value_in_unit(angstrom))
ca0 = start[ca_idx]


def snap():
    pos = np.array(sim.context.getState(getPositions=True).getPositions().value_in_unit(angstrom))
    ca = pos[ca_idx]
    rmsd = kabsch_rmsd(ca, ca0)
    ctr = ca.mean(0); rad = float(np.linalg.norm(ca - ctr, axis=1).mean())  # breathing proxy
    return round(rmsd, 2), round(rad, 1)


sim.context.setVelocitiesToTemperature(300 * kelvin)
r0, rad0 = snap()
traj = [{'ps': 0.0, 'ca_rmsd': r0, 'ring_radius': rad0}]
chunk = NSTEPS // 40
for i in range(40):
    sim.step(chunk)
    r, rad = snap()
    traj.append({'ps': round((i + 1) * chunk * 0.002, 1), 'ca_rmsd': r, 'ring_radius': rad})
    json.dump({'tag': TAG, 'natoms': system.getNumParticles(), 'nca': len(ca_idx),
               'traj': traj, 'elapsed': round(time.time() - t0, 1), 'platform': plat.getName()},
              open(f'/content/{TAG}.md.json', 'w'))
open(f'/content/{TAG}.md.DONE', 'w').write('done')
print('MD_DONE', TAG, traj[-1], f'{round(time.time()-t0)}s')
