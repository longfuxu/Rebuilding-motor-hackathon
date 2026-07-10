#!/usr/bin/env python3
"""Buildability atlas — uniform topology descriptors for ring ATPase motors.

For each panel protein we measure, on the native ring, the two geometric
quantities that (our N=4 observation says) dictate the best single-chain method:
  (1) how far each terminus sits from the functional channel axis  -> "jams channel?"
  (2) the direct head-to-tail C(i)->N(i+1) gap between adjacent subunits.

Everything is derived from coordinates (first/last modeled CA = the native
termini); the channel axis = ring symmetry axis (least-variance normal of the
subunit centroids). Where a substrate (DNA/RNA/peptide) is modeled in the pore we
also report terminus->substrate min distance (the direct "does it touch the
translocated polymer" test, as used for gp16 res330 -> DNA 6.2 A).

Output: descriptors.csv + descriptors.json  (+ printed table).
"""
import os, json, warnings
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser
warnings.simplefilter("ignore")

BASE = "/Users/longfu/Developer/claude-science-hackthon/.claude/worktrees/export-gp16-results/gp16_design/outputs/buildability_atlas"
STRUCT = os.path.join(BASE, "structures")

AA3 = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS",
       "MET","PHE","PRO","SER","THR","TRP","TYR","VAL","MSE","SEC","PYL"}
NUC = {"DA","DT","DG","DC","DU","A","U","G","C","I"}

# ---- panel config -----------------------------------------------------------
# ring   = list of (model,chain) subunits forming the functional ring
# subA/B = substrate chains (nucleic or peptide) modeled in the channel
# family, oligo (nominal functional oligomeric state), known method (for the 4 anchors)
# cterm_motor: optional residue number to use as the C-terminus of the *motor*
#              construct when a non-motor accessory domain extends past it.
PANEL = [
 dict(name="gp16",  pdb="7JQQ.pdb1", family="ASCE packaging", oligo=5,
      ring=[(0,c) for c in "ABCDE"], substrate=[(0,"F"),(0,"G")], known="CP",
      note="phi29 packaging motor; C-term res330 contacts dsDNA"),
 dict(name="ClpX",  pdb="6PP5.pdb1", family="AAA+ unfoldase", oligo=6,
      ring=[(0,c) for c in "ABCDEF"], substrate=[(0,"S")], known="direct",
      note="ClpXdN spiral; substrate peptide chain S"),
 dict(name="gp17",  pdb="3EZK.pdb1", family="ASCE terminase", oligo=5,
      ring=[(0,c) for c in "ABCDE"], substrate=[], known="direct",
      motor_range=(10,360), note="T4 terminase; CA-only; ATPase motor = res 10-360, nuclease 360-562 appended (dropped in the single-chain motor construct)"),
 dict(name="Rho",   pdb="3ICE.pdb1", family="RecA RNA translocase", oligo=6,
      ring=[(0,c) for c in "ABCDEF"], substrate=[(0,"G")], known=None,
      motor_range=(175,414), note="E. coli Rho closed ring + polyU RNA; motor = res 175-414 (N-terminal RNA-binding domain 1-130 dropped; matches the built rho_direct construct)"),
 dict(name="FtsK",  pdb="6T8B.pdb1", family="ASCE DNA translocase", oligo=6,
      ring=[(0,c) for c in "ABCDEF"], substrate=[(0,"G"),(0,"H")], known=None,
      note="P. aeruginosa FtsK motor (alpha/beta) + dsDNA"),
 dict(name="spastin", pdb="6P07.pdb1", family="AAA+ severing (meiotic clade)", oligo=6,
      ring=[(0,c) for c in "ABCDEF"], substrate=[(0,"G")], known=None,
      note="Drosophila spastin AAA hexamer + substrate peptide"),
 dict(name="katanin", pdb="6UGD.pdb1", family="AAA+ severing (meiotic clade)", oligo=6,
      ring=[(0,c) for c in "ABCDEF"], substrate=[(0,"G")], known=None,
      note="katanin p60 AAA hexamer + substrate peptide"),
 dict(name="p97_VCP", pdb="5FTN.pdb1", family="AAA+ classic double-ring", oligo=6,
      ring=[(0,c) for c in "ABCDEF"], substrate=[], known=None,
      note="human p97; N-D1-D2, two stacked AAA rings"),
 dict(name="DnaB",  pdb="4ESV.pdb1", family="RecA replicative helicase", oligo=6,
      ring=[(0,c) for c in "ABCDEF"], substrate=[(0,"V")], known=None,
      note="Aquifex DnaB hexamer + ssDNA"),
 dict(name="SV40_LTag", pdb="1SVM.pdb1", family="SF3 AAA+ helicase", oligo=6,
      ring=[(0,c) for c in "ABCDEF"], substrate=[], known=None,
      note="SV40 large T antigen helicase domain hexamer"),
 dict(name="HslU",  pdb="1DO0.pdb1", family="AAA+ unfoldase", oligo=6,
      ring=[(0,c) for c in "ABCDEF"], substrate=[], known=None,
      note="E. coli HslU hexamer (ClpX paralog)"),
 dict(name="Vps4",  pdb="6BMF.pdb1", family="AAA+ ESCRT", oligo=6,
      ring=[(0,c) for c in "ABCDE"], substrate=[(0,"G")], known=None,
      note="Vps4 spiral; 5 of 6 subunits ordered + substrate peptide"),
 dict(name="T7_gp4", pdb="1E0J.pdb1", family="SF4 helicase", oligo=6,
      ring=[(0,"A"),(0,"B"),(0,"C"),(1,"A"),(1,"B"),(1,"C")], substrate=[], known=None,
      note="T7 gp4 helicase domain hexamer (2 models x 3 chains)"),
]

def load(path):
    p = MMCIFParser(QUIET=True) if path.endswith(".cif") else PDBParser(QUIET=True, PERMISSIVE=True)
    return p.get_structure("s", path)

def get_subunit_atoms(struct, model_id, chain_id):
    """return dict resnum -> {'CA':xyz, 'all':[xyz,...]} for protein residues."""
    model = list(struct)[model_id]
    ch = model[chain_id]
    out = {}
    for res in ch:
        nm = res.resname.strip()
        if nm in AA3 and "CA" in res:
            allc = [a.coord for a in res]
            out[res.id[1]] = {"CA": res["CA"].coord.astype(float),
                              "all": np.array(allc, float)}
    return out

def substrate_atoms(struct, subs):
    pts = []
    for (mi, cid) in subs:
        try:
            ch = list(struct)[mi][cid]
        except Exception:
            continue
        for res in ch:
            for a in res:
                pts.append(a.coord)
    return np.array(pts, float) if pts else None

def ring_axis(centroids):
    C = np.array(centroids); ctr = C.mean(0); X = C - ctr
    _, _, Vt = np.linalg.svd(X)
    return ctr, Vt[2]

def dist_to_axis(p, ctr, axis):
    v = p - ctr
    return float(np.linalg.norm(v - np.dot(v, axis) * axis))

def angular_order(labels, centroids, ctr, axis):
    X = np.array(centroids) - ctr
    Xp = X - np.outer(X @ axis, axis)
    _, _, Vt = np.linalg.svd(Xp)
    u, v = Vt[0], Vt[1]
    ang = [np.arctan2(np.dot(x, v), np.dot(x, u)) for x in Xp]
    order = list(np.argsort(ang))
    return [labels[i] for i in order]

def min_dist(a_pts, b_pts):
    if a_pts is None or b_pts is None or len(a_pts)==0 or len(b_pts)==0:
        return None
    from scipy.spatial import cKDTree
    return float(cKDTree(b_pts).query(a_pts)[0].min())

def analyze(entry):
    path = os.path.join(STRUCT, entry["pdb"])
    st = load(path)
    subs = {}
    for (mi, cid) in entry["ring"]:
        atoms = get_subunit_atoms(st, mi, cid)
        if atoms:
            subs[(mi, cid)] = atoms
    labels = list(subs.keys())
    centroids = [np.mean([subs[l][r]["CA"] for r in subs[l]], 0) for l in labels]
    ctr, axis = ring_axis(centroids)
    radii = [dist_to_axis(c, ctr, axis) for c in centroids]
    ring_radius = float(np.mean(radii))
    order = angular_order(labels, centroids, ctr, axis)

    subpts = substrate_atoms(st, entry.get("substrate", []))
    mr = entry.get("motor_range")   # (lo,hi) for multi-domain proteins; else None -> full modeled

    # per-subunit termini. "primary" = motor construct terminus (motor_range if given,
    # else first/last modeled). "full" = first/last modeled chain (secondary column).
    def term_res(l, which):
        resn = sorted(subs[l])
        if mr:
            lo, hi = mr
            if which == "N":
                return lo if lo in subs[l] else min(r for r in resn if r >= lo)
            return hi if hi in subs[l] else max(r for r in resn if r <= hi)
        return resn[0] if which == "N" else resn[-1]

    Nax, Cax, Nax_full, Cax_full = [], [], [], []
    Nsub, Csub = [], []
    first_res, last_res = [], []
    for l in labels:
        resn = sorted(subs[l])
        first_res.append(resn[0]); last_res.append(resn[-1])
        nP, cP = term_res(l, "N"), term_res(l, "C")
        Nax.append(dist_to_axis(subs[l][nP]["CA"], ctr, axis))
        Cax.append(dist_to_axis(subs[l][cP]["CA"], ctr, axis))
        Nax_full.append(dist_to_axis(subs[l][resn[0]]["CA"], ctr, axis))
        Cax_full.append(dist_to_axis(subs[l][resn[-1]]["CA"], ctr, axis))
        if subpts is not None:
            Nsub.append(min_dist(subs[l][nP]["all"], subpts))
            Csub.append(min_dist(subs[l][cP]["all"], subpts))

    # direct C(i)->N(i+1) gaps, both ring senses; keep the smaller-mean (correct) sense
    def gaps(sense):
        g = []
        for i, l in enumerate(order):
            nxt = order[(i + sense) % len(order)]
            g.append(float(np.linalg.norm(subs[l][term_res(l,"C")]["CA"] - subs[nxt][term_res(nxt,"N")]["CA"])))
        return g
    gf, gr = gaps(1), gaps(-1)
    gap = gf if np.mean(gf) <= np.mean(gr) else gr
    # full-chain gap (first/last modeled), same sense choice
    def gaps_full(sense):
        g = []
        for i, l in enumerate(order):
            nxt = order[(i + sense) % len(order)]
            g.append(float(np.linalg.norm(subs[l][sorted(subs[l])[-1]]["CA"] - subs[nxt][sorted(subs[nxt])[0]]["CA"])))
        return g
    gff, gfr = gaps_full(1), gaps_full(-1)
    gap_full = gff if np.mean(gff) <= np.mean(gfr) else gfr

    mean = lambda a: float(np.mean([x for x in a if x is not None])) if any(x is not None for x in a) else None
    res = dict(
        name=entry["name"], family=entry["family"], pdb=entry["pdb"].split(".")[0],
        oligo_nominal=entry["oligo"], n_modeled=len(labels),
        ring_radius_A=round(ring_radius, 1),
        motor_range=(f"{mr[0]}-{mr[1]}" if mr else "full-modeled"),
        Nterm_to_axis_A=round(mean(Nax), 1),
        Cterm_to_axis_A=round(mean(Cax), 1),
        Nterm_rel_radius=round(mean(Nax)/ring_radius, 2),
        Cterm_rel_radius=round(mean(Cax)/ring_radius, 2),
        Cterm_to_axis_fullchain_A=round(mean(Cax_full), 1),
        Nterm_to_substrate_A=(round(mean(Nsub),1) if Nsub else None),
        Cterm_to_substrate_A=(round(mean(Csub),1) if Csub else None),
        direct_gap_A=round(mean(gap), 1),
        direct_gap_max_A=round(float(np.max(gap)), 1),
        direct_gap_fullchain_A=round(mean(gap_full), 1),
        first_res=int(np.median(first_res)), last_res=int(np.median(last_res)),
        known_method=entry.get("known"),
        note=entry.get("note",""),
    )
    return res

def main():
    rows = []
    for e in PANEL:
        try:
            rows.append(analyze(e))
            r = rows[-1]
            print(f"{r['name']:10} {r['pdb']:6} n={r['n_modeled']}/{r['oligo_nominal']} "
                  f"R={r['ring_radius_A']:5} mtr={r['motor_range']:12} | Nrel={r['Nterm_rel_radius']:4} "
                  f"Cax={r['Cterm_to_axis_A']:5} Crel={r['Cterm_rel_radius']:4} | Csub={r['Cterm_to_substrate_A']} "
                  f"| gap={r['direct_gap_A']:5} (full {r['direct_gap_fullchain_A']}) | known={r['known_method']}")
        except Exception as ex:
            import traceback; traceback.print_exc()
            print(f"FAIL {e['name']}: {ex}")
    with open(os.path.join(BASE, "descriptors.json"), "w") as f:
        json.dump(rows, f, indent=2)
    # CSV
    import csv
    keys = list(rows[0].keys())
    with open(os.path.join(BASE, "descriptors.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
        for r in rows: w.writerow(r)
    print(f"\nwrote descriptors.csv / .json ({len(rows)} rows)")

if __name__ == "__main__":
    main()
