#!/usr/bin/env python3
"""Roll up outputs/tiled_fold/*.result.json into summary.csv + SUMMARY.md.
Usage: python summarize.py [outdir]   (default ../../outputs/tiled_fold)"""
import sys, os, json, glob, re, csv

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "..", "outputs", "tiled_fold")


def parse(d):
    r = d.get("ring")
    if r:
        eng, ntot, hand = r["engaged"], r["n"], r["handedness"]
    else:
        m2 = re.search(r"M2:\s*(\d+)/(\d+)", d.get("M2", ""))
        eng, ntot = (int(m2.group(1)), int(m2.group(2))) if m2 else ("", "")
        hand = ""
    g = {"hand": hand}
    for k, pat in [("radius", r"radius ([\d.]+)"), ("cv", r"radius_CV ([\d.]+)"),
                   ("planarity", r"planarity_rms ([\d.]+)"), ("compact", r"compact_ring (\w+)"),
                   ("seq", r"sequential_consistent: (\w+)")]:
        mm = re.search(pat, d.get("M1", "")); g[k] = mm.group(1) if mm else ""
    m4 = re.search(r"pLDDT mean ([\d.]+)", d.get("M4", ""))
    return eng, ntot, g, (m4.group(1) if m4 else "")


ORDER = ["cp233_WT", "cp233_E119Q_1seat", "cp233_E119Q_5seat",
         "cp285_int15_inter10", "cp297_int15_inter10",
         "rfdiff_ring_L50_d0", "rfdiff_ring_L50_d1"]
rows = {}
for f in glob.glob(os.path.join(OUT, "*.result.json")):
    rows[json.load(open(f))["name"]] = json.load(open(f))
names = [n for n in ORDER if n in rows] + [n for n in rows if n not in ORDER]

recs = []
for n in names:
    d = rows[n]
    if not d.get("ok"):
        recs.append([n, d.get("construct_len", ""), "FAILED", "", "", "", "", "", "", d.get("err", "")[:60]])
        continue
    eng, ntot, g, m4 = parse(d)
    m2pass = "PASS" if eng != "" and eng >= 4 else "fail"
    recs.append([n, d.get("construct_len", ""), f"{eng}/{ntot}", m2pass, g["hand"],
                 g["radius"], g["cv"], g["compact"], m4, f"{d.get('wall_s','')}s"])

with open(os.path.join(OUT, "summary.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["construct", "len_aa", "M2_engaged", "M2", "handedness", "M1_radius_A",
                "M1_radius_CV", "M1_compact", "M4_iface_pLDDT", "wall"])
    w.writerows(recs)

hdr = ["construct", "len", "M2(ring)", "pass", "handedness", "radius(Å)", "CV", "compact",
       "ifacepLDDT", "wall"]
md = ["| " + " | ".join(hdr) + " |", "|" + "|".join(["---"] * len(hdr)) + "|"]
for r in recs:
    md.append("| " + " | ".join(str(x) for x in r) + " |")
header = (
    "# tiled_msa_fold results — gp16 single-chain rings (Boltz-2 NIM, tiled block-diagonal MSA)\n\n"
    "**M2(ring)** = trans-R146 arginine-finger -> sequential-neighbour Walker-A < 8 Å; pass >=4/5.\n"
    "Scored handedness-robust: a near-C5 ring can wind either way in one diffusion sample, so\n"
    "engagement is counted toward the DESIGNED neighbour (k->k+1) or its MIRROR (k->k-1),\n"
    "whichever is coherent. **M1** = ring geometry (compact + planar + sequential register).\n"
    "Never global pTM. Depth-1167 monomer MSA, tiled to 5836 rows, ~55-70 s/fold.\n\n"
    "## Key findings\n"
    "- **All 7 constructs close (M2 5/5, compact planar ring).**\n"
    "- **Dead seat does NOT open the ring:** cp233 WT, E119Q_1seat and E119Q_5seat are all\n"
    "  5/5 with all R146 fingers engaged at ~5.8-6.2 Å and identical geometry (radius ~26 Å,\n"
    "  CV ~0.01). The naive direction-specific M2 shows WT 5/5-forward but E119Q 0/5-forward;\n"
    "  that is a *ring-handedness* artifact — the E119Q rings are equally closed (5/5 reverse).\n"
    "  Confirmed by 3x replicates each (handedness_replicates.json): WT wound\n"
    "  designed/designed/mirror, E119Q wound mirror/designed/designed — i.e. BOTH sequences\n"
    "  wind both ways, so handedness is a per-sample coin-flip decoupled from the mutation.\n"
    "  Across all 13 folds (7 constructs + 6 replicates) M2(ring) = 5/5 every time.\n"
    "- **Tiled MSA is required:** the RFdiffusion native-order rings score 0/5 single-seq but\n"
    "  **5/5 tiled** (d0, d1) — single-seq cannot evaluate these single-chain rings.\n"
    "- Secondary CP sites **cp285 and cp297 also close 5/5** under tiled MSA.\n\n")
open(os.path.join(OUT, "SUMMARY.md"), "w").write(header + "\n".join(md) + "\n")
print("\n".join(md))
print("\nwrote", os.path.join(OUT, "summary.csv"), "and SUMMARY.md")
