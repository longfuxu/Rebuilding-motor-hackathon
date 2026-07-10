#!/usr/bin/env python3
"""Atlas decision-boundary figure: head-to-tail gap (x) vs substrate-channel fouling (y),
each ring motor a point coloured by predicted single-chain method. Anchors (validated)
are ringed; the fold-probed points (Rho/gp16 direct) are annotated."""
import os, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
rows = {r["name"]: r for r in json.load(open(os.path.join(BASE, "atlas_predictions.csv.json")))} \
    if os.path.exists(os.path.join(BASE, "atlas_predictions.csv.json")) else None

# read atlas_predictions.csv
import csv
rows = list(csv.DictReader(open(os.path.join(BASE, "atlas_predictions.csv"))))
GAP_DIRECT = 38.0
COL = {"CP": "#d1495b", "direct": "#2e8b57", "diffusion": "#3a6ea5"}
ANCHORS = {"gp16", "ClpX", "gp17"}

fig, ax = plt.subplots(figsize=(10.5, 5.6))
# manual y-offsets to declutter the off-channel cluster
YJIT = {"ClpX":0.10,"gp17":-0.10,"HslU":0.0,"Rho":0.10,"katanin":-0.12,"spastin":0.12,
        "FtsK":-0.02,"T7_gp4":0.24,"Vps4":-0.24,"p97_VCP":-0.12,"DnaB":0.10,"SV40_LTag":-0.02,"gp16":0.0}
LOFF = {"ClpX":(0,11),"gp17":(-2,-16),"HslU":(0,12),"Rho":(-6,11),"katanin":(0,-16),
        "spastin":(-4,12),"FtsK":(10,12),"T7_gp4":(6,10),"Vps4":(0,-16),"p97_VCP":(-10,-16),
        "DnaB":(6,11),"SV40_LTag":(0,12),"gp16":(0,13)}
for r in rows:
    x = float(r["direct_gap_A"]); jam = r["jams_channel"] == "True"
    y = (1.0 if jam else 0.0) + YJIT.get(r["name"], 0.0)
    m = r["predicted_method"]
    ax.scatter(x, y, s=175, c=COL[m], edgecolor="k",
               linewidth=(2.4 if r["name"] in ANCHORS else 0.6), zorder=3)
    ax.annotate(r["name"], (x, y), textcoords="offset points", xytext=LOFF.get(r["name"], (0,11)),
                ha="center", fontsize=8.5, fontweight=("bold" if r["name"] in ANCHORS else "normal"))

ax.axvline(GAP_DIRECT, color="gray", ls="--", lw=1)
ax.text(GAP_DIRECT+0.5, 0.60, f"gap {GAP_DIRECT} Å\n(direct | longer)", fontsize=8, color="gray")
ax.axhline(0.5, color="gray", ls=":", lw=1)

# region labels in the empty mid-band (y~0.5)
ax.text(23, 0.48, "DIRECT fusion", fontsize=12, color=COL["direct"], ha="center", fontweight="bold", alpha=.9)
ax.text(60, 0.42, "'diffusion' class\n(0 validated; Rho unsupported)", fontsize=10, color=COL["diffusion"], ha="center", fontweight="bold", alpha=.8)
ax.text(57, 1.28, "CIRCULAR PERMUTATION", fontsize=12, color=COL["CP"], ha="center", fontweight="bold", alpha=.9)

# annotate the fold probes
ax.annotate("Rho-direct: apo Boltz 6/6\n(suggestive, non-discriminating)", (46.1, 0.10),
            xytext=(40, -0.40), fontsize=7.5, ha="center", color=COL["diffusion"],
            arrowprops=dict(arrowstyle="->", color=COL["diffusion"], lw=1))
ax.annotate("gp16-direct control: apo Boltz 5/5 too,\nyet method is WRONG (AF3 0/5, C-term fouls DNA)", (57.1, 1.0),
            xytext=(35, 0.80), fontsize=7.5, ha="center", color=COL["CP"],
            arrowprops=dict(arrowstyle="->", color=COL["CP"], lw=1))

ax.set_yticks([0, 1]); ax.set_yticklabels(["termini\noff-channel", "terminus fouls\nchannel (jams)"])
ax.set_xlabel("direct head-to-tail C(i)→N(i+1) gap of the motor construct (Å)")
ax.set_title("Buildability atlas: substrate-channel topology decides the single-chain method for ring ATPases\n"
             "(bold-ringed = method-validated anchors gp16/ClpX/gp17; N=13 homo-oligomeric rings)", fontsize=10.5)
ax.set_ylim(-0.55, 1.5); ax.set_xlim(15, 78)
from matplotlib.lines import Line2D
leg = [Line2D([0],[0], marker='o', color='w', markerfacecolor=COL[k], markeredgecolor='k', markersize=11, label=k)
       for k in ("direct","CP","diffusion")]
leg.append(Line2D([0],[0], marker='o', color='w', markerfacecolor='w', markeredgecolor='k', markeredgewidth=2.2, markersize=11, label="validated anchor"))
ax.legend(handles=leg, loc="upper right", fontsize=8.5, framealpha=.95)
plt.tight_layout()
out = os.path.join(BASE, "atlas_decision_boundary.png")
plt.savefig(out, dpi=155)
print("wrote", out)
