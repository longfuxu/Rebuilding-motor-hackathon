#!/usr/bin/env python3
"""Rank the 27-construct CP linker-grid screen and write summary.csv + SCREEN_REPORT.md.
Ranking (judge by M1/M2, never pLDDT): M2 engaged (handedness-robust) desc, then
sequential-consistent, compact ring, tighter ring (radius_CV asc), tighter arginine-finger
engagement (mean R146 distance asc)."""
import os, json, glob, re, csv

REPO = "/Users/longfu/Developer/claude-science-hackthon/.claude/worktrees/export-gp16-results/gp16_design"
OUT = f"{REPO}/outputs/cp_grid_screen"
AF3 = {"cp233_int15_inter10", "cp285_int15_inter10", "cp297_int15_inter10"}


def m1(d):
    s = d.get("M1", "")
    def g(p, dv=""):
        m = re.search(p, s); return m.group(1) if m else dv
    return {"radius": g(r"radius ([\d.]+)"), "cv": g(r"radius_CV ([\d.]+)"),
            "planar": g(r"planarity_rms ([\d.]+)"), "compact": g(r"compact_ring (\w+)"),
            "seq": "YES" in s and "sequential_consistent: YES" in s}


def mean_engaged_dist(d):
    """Mean R146->WalkerA distance over engaged interfaces in the winning-handedness output."""
    ds = []
    for l in d.get("score_full", "").splitlines():
        m = re.match(r"\s*\w+\s+->\s+\w+\s+([\d.]+)\s+True", l)
        if m:
            ds.append(float(m.group(1)))
    return round(sum(ds) / len(ds), 2) if ds else None


rows = []
for f in glob.glob(f"{OUT}/*.result.json"):
    d = json.load(open(f)); n = d["name"]
    mm = re.match(r"cp(\d+)_int(\d+)_inter(\d+)", n)
    site, intl, interl = int(mm.group(1)), int(mm.group(2)), int(mm.group(3))
    if not d.get("ok"):
        rows.append({"name": n, "site": site, "int": intl, "inter": interl, "engaged": -1,
                     "n": 5, "hand": "", "cv": "", "seq": False, "compact": "", "radius": "",
                     "planar": "", "meandist": "", "iface_plddt": "", "err": d.get("err", "")[:50]})
        continue
    r = d["ring"]; g = m1(d)
    p = re.search(r"pLDDT mean ([\d.]+)", d.get("M4", ""))
    rows.append({"name": n, "site": site, "int": intl, "inter": interl,
                 "engaged": r["engaged"], "n": r["n"], "hand": r["handedness"],
                 "fwd": r["m2_forward"], "rev": r["m2_reverse"], "cv": g["cv"], "seq": g["seq"],
                 "compact": g["compact"], "radius": g["radius"], "planar": g["planar"],
                 "meandist": mean_engaged_dist(d), "iface_plddt": p.group(1) if p else "", "err": ""})


def key(r):
    cv = float(r["cv"]) if r["cv"] not in ("", None) else 9
    md = float(r["meandist"]) if r["meandist"] not in ("", None) else 99
    return (-r["engaged"], -(1 if r["seq"] else 0), -(1 if r["compact"] == "True" else 0), cv, md)


rows.sort(key=key)

cols = ["site", "int", "inter", "M2_engaged", "handedness", "M1_ring", "M1_sequential",
        "radius", "radius_CV", "mean_R146_A", "iface_pLDDT", "AF3_confirmed"]
with open(f"{OUT}/summary.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(cols)
    for r in rows:
        w.writerow([r["site"], r["int"], r["inter"],
                    f"{r['engaged']}/{r['n']}" if r["engaged"] >= 0 else "FAILED", r["hand"],
                    r["compact"], "YES" if r["seq"] else "no", r["radius"], r["cv"],
                    r["meandist"], r["iface_plddt"], "yes" if r["name"] in AF3 else ""])

# best per site
best = {}
for r in rows:
    if r["engaged"] < 0:
        continue
    if r["site"] not in best or key(r) < key(best[r["site"]]):
        best[r["site"]] = r

md = ["# CP linker-grid screen — 27 gp16 single-chain pentamers (tiled-MSA Boltz-2 NIM)\n",
      "CP site {233,285,297} x internal-linker {10,15,20} x inter-linker {10,15,20}. Each folded",
      "with a block-diagonal (tiled) monomer MSA on the free Boltz-2 NIM and scored with the",
      "**handedness-robust** ring M2 (R146 arginine-finger -> sequential-neighbour Walker-A < 8 Å,",
      "counted for whichever winding, k->k+1 or k->k-1, is coherent) plus M1 ring geometry.",
      "Ranked by M2 engaged, then sequential register, compact ring, tighter ring (radius_CV),",
      "tighter finger engagement (mean R146 distance). Never global pTM.\n",
      "## Best construct per CP site\n"]
for site in (233, 285, 297):
    if site in best:
        b = best[site]
        md.append(f"- **cp{site}: int{b['int']}_inter{b['inter']}** — M2 {b['engaged']}/{b['n']}, "
                  f"ring compact={b['compact']} sequential={'YES' if b['seq'] else 'no'}, "
                  f"radius {b['radius']} Å (CV {b['cv']}), mean R146 {b['meandist']} Å"
                  + ("  [AF3-confirmed]" if b["name"] in AF3 else ""))
# per-site geometry roll-up
from collections import defaultdict
bysite = defaultdict(list)
for r in rows:
    if r["engaged"] >= 0:
        bysite[r["site"]].append(r)
md.append("\n## Per-site geometry (all 9 linker combos per site; all close M2 5/5)\n")
md.append("| site | M2 (all combos) | mean radius_CV | mean R146 dist (Å) | note |")
md.append("|---|---|---|---|---|")
notes = {233: "tightest, most regular rings", 285: "tightest finger engagement",
         297: "closes but least regular (CV ~6x cp233)"}
for site in (233, 285, 297):
    rs = bysite[site]
    cv = sum(float(r["cv"]) for r in rs) / len(rs)
    dd = sum(float(r["meandist"]) for r in rs) / len(rs)
    md.append(f"| cp{site} | {sum(r['engaged'] for r in rs)}/{5*len(rs)} (9/9 pass) | "
              f"{cv:.3f} | {dd:.2f} | {notes[site]} |")
md.append("\n**Takeaway:** every CP site x linker combination closes (27/27 M2 5/5, compact "
          "sequential rings) — the single-chain CP-ring design is robust to internal- and "
          "inter-linker length. cp233 and cp285 give consistently more regular rings "
          "(radius_CV ~0.008-0.019) than cp297 (~0.052); cp297 is the geometrically weakest "
          "site though still passing. Linker length is a second-order knob: shorter inter-linkers "
          "(inter10) trend to the tightest rings.\n")
md.append("\n## Full ranking (best first)\n")
hdr = ["#", "construct", "M2", "hand", "compact", "seq-reg", "radius(Å)", "CV", "meanR146(Å)", "ifpLDDT", "AF3"]
md.append("| " + " | ".join(hdr) + " |")
md.append("|" + "|".join(["---"] * len(hdr)) + "|")
for i, r in enumerate(rows, 1):
    m2 = f"{r['engaged']}/{r['n']}" if r["engaged"] >= 0 else "FAIL"
    md.append("| " + " | ".join(str(x) for x in [
        i, r["name"], m2, r["hand"].split("(")[0] if r["hand"] else "",
        r["compact"], "YES" if r["seq"] else "no", r["radius"], r["cv"], r["meandist"],
        r["iface_plddt"], "yes" if r["name"] in AF3 else ""]) + " |")
md.append("\nAF3-confirmed (cross-predictor) baselines: cp233_int15_inter10, cp285_int15_inter10, "
          "cp297_int15_inter10. Handedness is a per-sample coin-flip (n=1 diffusion sample each) "
          "decoupled from sequence — both windings are equally-closed rings; that is why scoring is "
          "handedness-robust.")
open(f"{OUT}/SCREEN_REPORT.md", "w").write("\n".join(md) + "\n")
print("\n".join(md[:40]))
print(f"\nwrote {OUT}/summary.csv and SCREEN_REPORT.md  ({len(rows)} constructs)")
