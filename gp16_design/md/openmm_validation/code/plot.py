#!/usr/bin/env python3
"""Plot gp16 ring MD readouts for A (apo), B (7JQQ helical), C (design).
Reads results/<sys>/<sys>_{timeseries,contacts,rmsf}.csv + _summary.json.
Colorblind-safe: A=blue (apo), B=red (7JQQ ATP-helical), C=green (design)."""
import os, csv, json, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = sys.argv[1] if len(sys.argv) > 1 else "results"
OUT = sys.argv[2] if len(sys.argv) > 2 else RES
COL = {"A": "#2166AC", "B": "#B2182B", "C": "#1B7837"}
LAB = {"A": "A  apo closed ring (Boltz)",
       "B": "B  7JQQ ATP-bound helical",
       "C": "C  design cp233_int15_inter10"}

def read_ts(s):
    p = os.path.join(RES, s, f"{s}_timeseries.csv")
    if not os.path.exists(p): return None
    rows = list(csv.DictReader(open(p)))
    d = {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
    return d

def read_ct(s):
    p = os.path.join(RES, s, f"{s}_contacts.csv")
    if not os.path.exists(p): return None
    rows = list(csv.DictReader(open(p)))
    t = np.array([float(r["time_ns"]) for r in rows])
    tot = np.array([float(r["nc_total"]) for r in rows])
    return t, tot

def read_rmsf(s):
    p = os.path.join(RES, s, f"{s}_rmsf.csv")
    if not os.path.exists(p): return None
    rows = list(csv.DictReader(open(p)))
    return (np.array([int(r["idx"]) for r in rows]),
            np.array([float(r["rmsf_A"]) for r in rows]))

TS = {s: read_ts(s) for s in "ABC"}
have = [s for s in "ABC" if TS[s] is not None]
print("systems with data:", have)

# ---------- Figure 1: overview 2x3 ----------
fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))
fig.suptitle("gp16 ring physics validation — OpenMM implicit-solvent MD (GBSA-OBC2, A100)",
             fontsize=13, fontweight="bold")

def plot_line(a, key, ylabel, title, band=False):
    for s in have:
        d = TS[s]
        a.plot(d["time_ns"], d[key], color=COL[s], lw=1.6, label=LAB[s])
    a.set_xlabel("time (ns)"); a.set_ylabel(ylabel); a.set_title(title)
    a.grid(alpha=0.25)

plot_line(ax[0,0], "ca_rmsd_A", "Cα-RMSD (Å)", "(a) backbone drift vs t=0")
plot_line(ax[0,1], "radius_CV", "radius CV", "(b) radial symmetry (0 = perfect ring)")
ax[0,1].axhline(0.35, color="gray", ls="--", lw=1, label="compact-ring threshold")

# (c) mean sequential-interface R146->WalkerA distance ± spread
axc = ax[0,2]
for s in have:
    d = TS[s]
    D = np.vstack([d[f"d_if{k+1}_A"] for k in range(5)])
    m = D.mean(0); lo = D.min(0); hi = D.max(0)
    axc.plot(d["time_ns"], m, color=COL[s], lw=1.8, label=LAB[s])
    axc.fill_between(d["time_ns"], lo, hi, color=COL[s], alpha=0.15)
axc.axhline(8.0, color="gray", ls="--", lw=1)
axc.text(0.02, 8.3, "engaged < 8 Å", fontsize=8, color="gray", transform=axc.get_yaxis_transform())
axc.set_xlabel("time (ns)"); axc.set_ylabel("R146→WalkerA min-dist (Å)")
axc.set_title("(c) ring closure — mean±range of 5 interfaces")

plot_line(ax[1,0], "n_engaged", "# engaged interfaces (<8 Å)", "(d) # closed interfaces")
ax[1,0].set_ylim(-0.3, 5.3)

# (e) RMSF
axe = ax[1,1]
for s in have:
    r = read_rmsf(s)
    if r is None: continue
    axe.plot(r[0], r[1], color=COL[s], lw=1.0, label=LAB[s])
axe.set_xlabel("Cα index (concatenated subunits)"); axe.set_ylabel("RMSF (Å)")
axe.set_title("(e) per-residue flexibility")

# (f) interface contacts
axf = ax[1,2]
for s in have:
    c = read_ct(s)
    if c is None: continue
    axf.plot(c[0], c[1], color=COL[s], lw=1.6, label=LAB[s])
axf.set_xlabel("time (ns)"); axf.set_ylabel("inter-subunit residue contacts (<4.5 Å)")
axf.set_title("(f) total interface contacts")

ax[0,0].legend(fontsize=8, loc="upper left")
fig.tight_layout(rect=[0,0,1,0.96])
f1 = os.path.join(OUT, "fig1_overview.png")
fig.savefig(f1, dpi=150); print("wrote", f1)

# ---------- Figure 2: per-interface closure small multiples ----------
fig2, ax2 = plt.subplots(1, 3, figsize=(15, 4.2), sharey=True)
for j, s in enumerate(have):
    d = TS[s]
    for k in range(5):
        ax2[j].plot(d["time_ns"], d[f"d_if{k+1}_A"], lw=1.2, label=f"iface {k+1}")
    ax2[j].axhline(8.0, color="gray", ls="--", lw=1)
    ax2[j].set_title(LAB[s]); ax2[j].set_xlabel("time (ns)"); ax2[j].grid(alpha=0.25)
    ax2[j].legend(fontsize=7, ncol=2)
ax2[0].set_ylabel("R146→WalkerA min-dist (Å)")
fig2.suptitle("Per-interface ring-closure traces (each line = one sequential seam)", fontweight="bold")
fig2.tight_layout(rect=[0,0,1,0.94])
f2 = os.path.join(OUT, "fig2_per_interface.png")
fig2.savefig(f2, dpi=150); print("wrote", f2)

# ---------- summary table ----------
print("\n=== SUMMARY (last-20% means) ===")
hdr = f"{'sys':>3} {'t_ns':>6} {'RMSD_f':>7} {'nEng_t0':>7} {'nEng_f':>7} {'radCV_t0':>8} {'radCV_f':>8}"
print(hdr)
rowsout = []
for s in have:
    p = os.path.join(RES, s, f"{s}_summary.json")
    j = json.load(open(p)) if os.path.exists(p) else {}
    line = (f"{s:>3} {j.get('t_total_ns',0):>6.2f} {j.get('ca_rmsd_final',0):>7.2f} "
            f"{j.get('n_engaged_t0',0):>7} {j.get('n_engaged_mean_last20pct',0):>7.2f} "
            f"{j.get('radius_CV_t0',0):>8.3f} {j.get('radius_CV_mean_last20pct',0):>8.3f}")
    print(line); rowsout.append((s,j))
json.dump({s: j for s,j in rowsout}, open(os.path.join(OUT,"summary_all.json"),"w"), indent=2)
