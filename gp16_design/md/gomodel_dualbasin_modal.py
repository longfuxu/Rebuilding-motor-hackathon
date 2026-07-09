"""Cα dual-basin structure-based (Go) model of the gp16 ring, on Modal GPU.

Two energy basins: apo (closed planar) and 7JQQ (3-ATP helical intermediate). Simulating from apo lets the ring
sample the planar<->helical opening that static predictors can't produce and all-atom MD can't reach.

Prep (already done, local): md/gomodel_inputs.npz  (Xa, Xj in nm; backbone bonds; dual-basin contacts).
Run (FOREGROUND only): modal run gp16_design/md/gomodel_dualbasin_modal.py   # token saved in ~/.modal.toml
  IMPORTANT: run in the FOREGROUND (do NOT launch as a detached/background shell job). This uses a synchronous
  .remote()/.map() call, which Modal CANCELS once the client stops polling (a background job drops polling) — that
  is exactly what happened on 2026-07-08 (call cancelled). For true background use, refactor run/sweep to .spawn()
  and poll the FunctionCall id later. Smoke test (`--n-steps 20000`) validated the pipeline end-to-end on GPU.
Output: md/gomodel_out.npz / gomodel_sweep.npz — inspect locally. STATUS: pipeline works (smoke test 20k steps OK),
  but the 6×1M-step scan HIT THE 3600s FUNCTION TIMEOUT — this Cα model (15,231 dual-basin contacts + a 15k-exclusion
  excluded-volume force) is SLOW, not cheap. Before a real run, OPTIMIZE: single (eps,temp), <=200-500k steps, and
  drop/cheapen the excluded-volume CustomNonbondedForce (or coarsen the contact set / raise the contact cutoff filter).

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
    def sweep(Xa, Xj, bonds, contacts, n, combos, n_steps):
        return {c: _simulate(Xa, Xj, bonds, contacts, n, c[0], c[1], n_steps) for c in combos}

    @app.function(gpu=GPU, image=image, timeout=3600)
    def run(Xa, Xj, bonds, contacts, n, eps, temp, n_steps):
        return _simulate(Xa, Xj, bonds, contacts, n, eps, temp, n_steps)

    @app.local_entrypoint()
    def main(n_steps: int = 1_000_000, eps: float = 0.0, temp: float = 0.0):
        d = np.load("gp16_design/md/gomodel_inputs.npz")
        args = (d["Xa"], d["Xj"], d["bonds"], d["contacts"], int(d["n"]))
        if eps > 0:                                   # single production run at a chosen (eps,temp)
            tr = run.remote(*args, eps, temp, n_steps)
            np.savez("gp16_design/md/gomodel_out.npz", traj=tr, eps=eps, temp=temp)
            print(f"# production eps={eps} temp={temp} n_steps={n_steps}")
            print("  step      opening  rmsd_apo  rmsd_7jqq   (want opening UP, rmsd_7jqq DOWN = apo->helical)")
            for s, o, ra, rj in tr[::max(1, len(tr)//25)]:
                print(f"  {int(s):>8}  {o:8.2f}  {ra:8.2f}  {rj:9.2f}")
        else:                                          # eps/temp scan
            combos = [(2.0, 300.0), (2.0, 450.0), (5.0, 300.0), (5.0, 450.0), (10.0, 450.0), (10.0, 700.0)]
            res = sweep.remote(*args, combos, n_steps)
            np.savez("gp16_design/md/gomodel_sweep.npz",
                     **{f"eps{c[0]}_T{c[1]}": tr for c, tr in res.items()})
            print(f"# eps/temp scan, {n_steps} steps each. start->end. want opening UP and rmsd_7jqq DOWN (apo->helical):")
            print(f"  {'eps':>5} {'temp':>5} | {'opening':>14} | {'rmsd_apo':>13} | {'rmsd_7jqq':>15} | toward7jqq")
            for c, tr in res.items():
                o0, oe = tr[0][1], tr[-1][1]; ra0, ra = tr[0][2], tr[-1][2]; rj0, rj = tr[0][3], tr[-1][3]
                print(f"  {c[0]:>5} {c[1]:>5} | {o0:5.1f} -> {oe:5.1f} | {ra0:4.1f} -> {ra:4.1f} | {rj0:5.1f} -> {rj:5.1f} | {rj0-rj:+5.1f}")
