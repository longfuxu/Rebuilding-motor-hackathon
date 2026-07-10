#!/usr/bin/env python3
"""gp17 topology analysis on the native pentameric motor assembly (RCSB 3EZK).

Mirrors the gp16 (7JQQ) / ClpX (6PP5) analysis: the two geometric quantities that
dictate the optimal single-chain topology are
  (1) how far the termini sit from the functional (DNA) channel axis, and
  (2) the direct head-to-tail C(i)->N(i+1) gap around the ring.

3EZK is a CA-only pseudo-atomic pentamer fit into a 34 A cryo-EM map (Sun et al.
2008, Cell), so all distances are CA-based and low-precision (report ~ +/- several A).
The pore/channel axis is the 5-fold symmetry axis = normal to the plane of the five
chain centroids, through the overall centroid (no DNA is modeled in 3EZK).
"""
import numpy as np

PDB = "structures/3EZK.pdb"
# gp17 landmark residues (native P17312 numbering; grounded from UniProt + 2O0H)
N_TERM = 10          # first modeled residue
ATPASE_C = 360       # C-terminus of the ATPase domain (Sun et al; 2O0H covers 1-357)
FULL_C = 562         # last modeled residue (full-length; nuclease tip 563-610 disordered)
WALKER_A = list(range(161, 168))   # 161-167
ARG_FINGER = 162                    # R162 candidate arginine finger

def load_ca(path):
    """chain -> {resnum: xyz(CA)}"""
    ch = {}
    for L in open(path):
        if L[:4] == "ATOM" and L[12:16].strip() == "CA":
            c = L[21]; rn = int(L[22:26])
            ch.setdefault(c, {})[rn] = np.array([float(L[30:38]), float(L[38:46]), float(L[46:54])])
    return ch

def ring_axis(centroids):
    C = np.array(centroids); ctr = C.mean(0); X = C - ctr
    _, _, Vt = np.linalg.svd(X)
    axis = Vt[2]                     # least-variance direction = ring normal = 5-fold axis
    return ctr, axis

def dist_to_axis(p, ctr, axis):
    v = p - ctr
    return float(np.linalg.norm(v - np.dot(v, axis) * axis))

def domain_centroid(resmap, lo, hi):
    pts = [resmap[r] for r in resmap if lo <= r <= hi]
    return np.mean(pts, 0)

ch = load_ca(PDB)
chains = sorted(ch)
print(f"# 3EZK chains: {chains}; modeled residues {min(min(m) for m in ch.values())}-{max(max(m) for m in ch.values())}")

# ---- axis from full-chain centroids AND from ATPase-domain-only centroids ----
full_centroids = [np.mean(list(ch[c].values()), 0) for c in chains]
atp_centroids = [domain_centroid(ch[c], 10, ATPASE_C) for c in chains]
ctr_full, axis_full = ring_axis(full_centroids)
ctr_atp, axis_atp = ring_axis(atp_centroids)

# angular order of chains around the ATPase-ring axis
def angular_order(centroids, ctr, axis):
    X = np.array(centroids) - ctr
    _, _, Vt = np.linalg.svd(X - np.outer(X @ axis, axis))
    u, v = Vt[0], Vt[1]
    ang = [np.arctan2(np.dot(x, v), np.dot(x, u)) for x in X]
    return [chains[i] for i in np.argsort(ang)]

order = angular_order(atp_centroids, ctr_atp, axis_atp)
print(f"# ATPase-ring angular order (spatial neighbours): {order}")

# ring radius (ATPase-domain centroids to axis)
radii = [dist_to_axis(pc, ctr_atp, axis_atp) for pc in atp_centroids]
print(f"# ATPase-ring radius (centroid->axis): mean {np.mean(radii):.1f} A  (per-chain {[round(r,1) for r in radii]})")

print("\n## (1) TERMINI vs DNA-CHANNEL AXIS  (min over chains; radial dist to 5-fold axis)")
for label, res, axis, ctr in [
    ("N-terminus res10        ", N_TERM, axis_atp, ctr_atp),
    ("ATPase C-terminus res360 ", ATPASE_C, axis_atp, ctr_atp),
    ("full-length C-term res562", FULL_C, axis_full, ctr_full),
    ("Walker-A K166            ", 166, axis_atp, ctr_atp),
    ("arg-finger R162          ", ARG_FINGER, axis_atp, ctr_atp),
]:
    ds = [dist_to_axis(ch[c][res], ctr, axis) for c in chains if res in ch[c]]
    if ds:
        print(f"  {label}: {np.mean(ds):5.1f} A mean  ({np.min(ds):.1f}-{np.max(ds):.1f})")

# ---- (2) direct head-to-tail C(i)->N(i+1) gaps, spatially adjacent chains ----
print("\n## (2) DIRECT C->N HEAD-TO-TAIL GAPS  (spatially adjacent chains, CA-CA)")
def neighbour_gap(c_res):
    gaps = []
    for i, c in enumerate(order):
        nxt = order[(i + 1) % len(order)]
        if c_res in ch[c] and N_TERM in ch[nxt]:
            g = float(np.linalg.norm(ch[c][c_res] - ch[nxt][N_TERM]))
            gaps.append((c, nxt, g))
    return gaps

for label, c_res in [("ATPase C res360 -> next N res10", ATPASE_C),
                     ("full C res562  -> next N res10", FULL_C)]:
    gaps = neighbour_gap(c_res)
    vals = [g for _, _, g in gaps]
    print(f"  {label}: mean {np.mean(vals):.1f} A  ({np.min(vals):.0f}-{np.max(vals):.0f})  per-iface {[f'{a}->{b}:{g:.0f}' for a,b,g in gaps]}")

# also report N-term(res10) -> own ATPase C(res360) intra-subunit span (context)
intra = [float(np.linalg.norm(ch[c][N_TERM] - ch[c][ATPASE_C])) for c in chains if N_TERM in ch[c] and ATPASE_C in ch[c]]
print(f"\n# intra-subunit N(res10)..C(res360) span: mean {np.mean(intra):.1f} A")
