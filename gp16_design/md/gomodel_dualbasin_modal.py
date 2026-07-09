"""Cα dual-basin structure-based (Go) model of the gp16 ring, on Modal GPU.

Two energy basins: apo (closed planar) and 7JQQ (3-ATP helical intermediate). Simulating from apo lets the ring
sample the planar<->helical opening that static predictors can't produce and all-atom MD can't reach.

Prep (already done, local): md/gomodel_inputs.npz  (Xa, Xj in nm; backbone bonds; dual-basin contacts).
Run:   modal run gp16_design/md/gomodel_dualbasin_modal.py          # after `modal token new`
Output: md/gomodel_out.npz (opening + RMSD-to-basins timeseries) — inspect locally.

STATUS: v1, syntax-checked locally (no OpenMM here). Tunables that likely need a sweep on first run:
eps (contact depth), T (temperature), and the run length. Start weak-ish so the ring can cross between basins.
"""
import numpy as np
try:
    import modal
    app = modal.App("gp16-gomodel-dualbasin")
    image = modal.Image.debian_slim(python_version="3.11").pip_install("openmm", "numpy")
    GPU = "H100"
except Exception:                       # allow local import / syntax check without modal installed
    modal = None

def _simulate(Xa, Xj, bonds, contacts, n, eps=2.5, temp=300.0, n_steps=5_000_000, report=10000):
    import openmm as mm, openmm.unit as u, numpy as np
    N = len(Xa)
    sysm = mm.System()
    for _ in range(N): sysm.addParticle(110.0 * u.amu)
    # backbone pseudobonds (r0 from apo)
    bb = mm.HarmonicBondForce()
    for i, j in bonds:
        r0 = float(np.linalg.norm(Xa[j] - Xa[i]))
        bb.addBond(int(i), int(j), r0 * u.nanometer, 20000.0 * u.kilojoule_per_mole / u.nanometer**2)
    sysm.addForce(bb)
    # dual-basin native contacts: wells at BOTH the apo (r1) and 7JQQ (r2) distances
    cb = mm.CustomBondForce("-eps*(exp(-(r-r1)^2/(2*w2)) + exp(-(r-r2)^2/(2*w2)))")
    cb.addGlobalParameter("eps", eps); cb.addGlobalParameter("w2", 0.05**2)   # 0.05 nm well width
    cb.addPerBondParameter("r1"); cb.addPerBondParameter("r2")
    for i, j, ra, rj in contacts:
        cb.addBond(int(i), int(j), [float(ra), float(rj)])
    sysm.addForce(cb)
    # excluded volume on non-native pairs
    ev = mm.CustomNonbondedForce("erep*(sig/r)^12")
    ev.addGlobalParameter("erep", 1.0); ev.addGlobalParameter("sig", 0.4)
    ev.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic); ev.setCutoffDistance(0.8 * u.nanometer)
    for _ in range(N): ev.addParticle([])
    for i, j in bonds: ev.addExclusion(int(i), int(j))
    for i, j, ra, rj in contacts: ev.addExclusion(int(i), int(j))
    sysm.addForce(ev)
    integ = mm.LangevinMiddleIntegrator(temp * u.kelvin, 1.0 / u.picosecond, 0.010 * u.picosecond)
    try:    plat = mm.Platform.getPlatformByName("CUDA")
    except Exception: plat = mm.Platform.getPlatformByName("Reference")
    ctx = mm.Context(sysm, integ, plat)
    ctx.setPositions(Xa * u.nanometer)
    mm.LocalEnergyMinimizer.minimize(ctx)
    # ring frame for the opening observable
    cent0 = np.array([Xa[k*n:(k+1)*n].mean(0) for k in range(5)]); ctr = cent0.mean(0)
    _, _, Vt = np.linalg.svd(cent0 - ctr); normal = Vt[2]
    out = []
    for s in range(0, n_steps, report):
        integ.step(report)
        X = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(u.nanometer)
        cent = np.array([X[k*n:(k+1)*n].mean(0) for k in range(5)])
        oop = np.sqrt((np.array([abs((cent[k]-cent.mean(0))@normal) for k in range(5)])**2).mean())
        rmsd_apo = np.sqrt(((X-Xa)**2).sum(1).mean()); rmsd_jqq = np.sqrt(((X-Xj)**2).sum(1).mean())
        out.append((s, oop, rmsd_apo, rmsd_jqq))
    return np.array(out)

if modal is not None:
    @app.function(gpu=GPU, image=image, timeout=3600)
    def run(Xa, Xj, bonds, contacts, n, eps, temp, n_steps):
        return _simulate(Xa, Xj, bonds, contacts, n, eps, temp, n_steps)

    @app.local_entrypoint()
    def main(eps: float = 2.5, temp: float = 300.0, n_steps: int = 5_000_000):
        d = np.load("gp16_design/md/gomodel_inputs.npz")
        res = run.remote(d["Xa"], d["Xj"], d["bonds"], d["contacts"], int(d["n"]), eps, temp, n_steps)
        np.savez("gp16_design/md/gomodel_out.npz", traj=res)
        print("opening (out-of-plane, nm) vs step — apo->helical if it rises and RMSD-to-7JQQ falls:")
        for s, oop, ra, rj in res[::max(1, len(res)//20)]:
            print(f"  step {int(s):>8}: opening {oop:.2f}  rmsd_apo {ra:.2f}  rmsd_7jqq {rj:.2f}")
