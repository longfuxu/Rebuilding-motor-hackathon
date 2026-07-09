#!/usr/bin/env python3
"""Compute stability discriminators per system and draft the one-sentence verdict.
Reads results/<sys>/<sys>_timeseries.csv (+ _contacts, _summary)."""
import os, csv, json, sys
import numpy as np

RES = sys.argv[1] if len(sys.argv) > 1 else "results"
LAB = {"A": "apo closed ring (Boltz)", "B": "7JQQ ATP-helical (ligands stripped)",
       "C": "design cp233_int15_inter10"}

def load(s):
    p = os.path.join(RES, s, f"{s}_timeseries.csv")
    if not os.path.exists(p): return None
    rows = list(csv.DictReader(open(p)))
    d = {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
    cp = os.path.join(RES, s, f"{s}_contacts.csv")
    if os.path.exists(cp):
        cr = list(csv.DictReader(open(cp)))
        d["_ct"] = np.array([float(r["nc_total"]) for r in cr])
    return d

def last(a, f=0.8): return a[int(f*len(a)):]

print(f"{'sys':>3} | {'RMSD f/plateau':>16} | {'radCV t0->f':>13} | {'nEng t0->f':>11} | "
      f"{'iface mean t0->f':>17} | {'contacts t0->f':>15}")
print("-"*95)
out = {}
for s in "ABC":
    d = load(s)
    if d is None:
        print(f"{s:>3} | (no data)"); continue
    n = len(d["time_ns"])
    rmsd = d["ca_rmsd_A"]; radcv = d["radius_CV"]; neng = d["n_engaged"]
    D = np.vstack([d[f"d_if{k+1}_A"] for k in range(5)])
    ifmean = D.mean(0)
    ct = d.get("_ct")
    rmsd_plateau = last(rmsd).mean()
    # linear slope of radius_CV over last 60% (drift up = opening)
    t = d["time_ns"]; sel = slice(int(0.4*n), n)
    slope_radcv = np.polyfit(t[sel], radcv[sel], 1)[0]  # per ns
    row = dict(
        rmsd_final=float(rmsd[-1]), rmsd_plateau=float(rmsd_plateau),
        radcv_t0=float(radcv[0]), radcv_final=float(last(radcv).mean()),
        radcv_slope_per_ns=float(slope_radcv),
        neng_t0=int(neng[0]), neng_final=float(last(neng).mean()),
        ifmean_t0=float(ifmean[0]), ifmean_final=float(last(ifmean).mean()),
        ct_t0=float(ct[0]) if ct is not None else None,
        ct_final=float(last(ct).mean()) if ct is not None else None,
        t_ns=float(t[-1]))
    out[s] = row
    ctxt = f"{row['ct_t0']:.0f}->{row['ct_final']:.0f}" if ct is not None else "n/a"
    print(f"{s:>3} | {row['rmsd_final']:5.2f}/{row['rmsd_plateau']:5.2f} Å    | "
          f"{row['radcv_t0']:.3f}->{row['radcv_final']:.3f} | "
          f"{row['neng_t0']}->{row['neng_final']:4.1f}    | "
          f"{row['ifmean_t0']:5.2f}->{row['ifmean_final']:5.2f} Å      | {ctxt:>15}")

# --- draft verdict ---
if all(k in out for k in "AC"):
    A, C = out["A"], out["C"]
    def closed(r):  # heuristic stability test
        return (r["neng_final"] >= 4.0 and r["radcv_final"] < 0.12
                and r["radcv_slope_per_ns"] < 0.02 and r["rmsd_plateau"] < 8.0)
    cA, cC = closed(A), closed(C)
    print("\n--- draft verdict inputs ---")
    print(f"A stable/closed: {cA} | C stable/closed: {cC}")
    ct_keep = (C["ct_final"]/C["ct_t0"]) if C.get("ct_t0") else float('nan')
    print(f"C interface-contact retention: {ct_keep*100:.0f}%  "
          f"radiusCV drift {C['radcv_slope_per_ns']*1000:.1f}e-3 /ns  "
          f"engaged {C['neng_t0']}->{C['neng_final']:.1f}/5")
    if cC and cA:
        v = ("PASS — under identical implicit-solvent MD the design ring stays closed and "
             "symmetric on par with the genuine apo ring (interfaces engaged, radius_CV low, "
             "RMSD plateaus); it does not open or collapse.")
    elif cA and not cC:
        v = ("FAIL/CONCERN — the apo ring stays closed but the design opens/destabilizes "
             "(interface dissociation / radius_CV drift / RMSD climb), so the design is not "
             "physically robust under MD.")
    else:
        v = ("INCONCLUSIVE — apo reference itself is not stable under this protocol; "
             "cannot cleanly judge the design (revisit protocol/length).")
    print("\nVERDICT:", v)
json.dump(out, open(os.path.join(RES, "verdict.json"), "w"), indent=2)
print("\nwrote", os.path.join(RES, "verdict.json"))
