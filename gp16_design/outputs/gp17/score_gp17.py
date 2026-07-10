#!/usr/bin/env python3
"""gp17 M1/M2 scorer for single-chain ATPase-domain ring constructs.

Adapts gp16_design/reproduce/score_m2.py to gp17 numbering and, importantly, to gp17's
UNSETTLED coupler biology. Unlike phi29-gp16 (clear TRANS arg-finger R146), the gp17
coupling arginine is debated cis vs trans (Rao lab candidates R162/R245/R321/R406; R162
is Walker-A-proximal). So we DO NOT rest the verdict on a trans-finger assumption:

  PRIMARY closure metric  = M1 sequential_consistent (biology-agnostic register test:
      does the ring's spatial order == the designed cyclic order) + compact/planar ring.
  SUPPORTING (reported, not decisive) = arg-finger geometry under BOTH readings:
      trans R162(i)->WalkerA(i+1), cis R162(i)->WalkerA(i), alt trans R245(i)->WalkerA(i+1),
      and whether the designed neighbour is the NEAREST subunit to R162(i) (scramble flag).

gp17 landmarks (native P17312 numbering == within-copy position, copy = res1-360):
  Walker-A 161-167 (K166) | arg-finger cand R162 | alt R245 | Walker-B/E256.
engaged threshold 8 A (same as gp16 M2).
"""
import sys, os
import numpy as np
sys.path.insert(0, "/Users/longfu/Developer/claude-science-hackthon/.claude/worktrees/export-gp16-results/gp16_design/reproduce")
from score_m2 import parse_atoms, subunit_index, atoms_for, min_dist, ca_centroid, ring_geometry

GUAN = {"NE", "CZ", "NH1", "NH2"}
WALKER_A = list(range(161, 168))   # 161-167
R162, R245 = 162, 245
ENGAGED = 8.0

def main():
    struct = sys.argv[1]
    copies_arg = sys.argv[2]          # e.g. "A:1-360,A:381-740,..."
    copies = []
    for i, tok in enumerate(copies_arg.split(",")):
        ch, span = tok.split(":"); lo, hi = span.split("-")
        copies.append((f"c{i+1}", ch, int(lo), int(hi)))
    atoms = parse_atoms(struct)
    subs = subunit_index(atoms, copies, copy_start_res=1)   # within-copy res == native numbering
    labels = list(subs); n = len(labels)

    def dist(donor_sub, res_don, names, acc_sub, res_acc):
        A, _ = atoms_for(subs[donor_sub], [res_don], names)
        B, _ = atoms_for(subs[acc_sub], res_acc)
        return min_dist(A, B)

    print(f"# {os.path.basename(struct)}  copies={labels}")
    print(f"{'i':>4} {'->i+1':>6} | {'transR162':>9} {'eng':>4} | {'cisR162':>8} | {'transR245':>9} | {'nearestR162':>16}")
    trans162 = cis162 = trans245 = 0
    for i, d in enumerate(labels):
        a = labels[(i+1) % n]                       # designed sequential neighbour
        t162 = dist(d, R162, GUAN, a, WALKER_A)
        c162 = dist(d, R162, GUAN, d, [166])        # cis: R162 -> own catalytic K166 (avoid self-overlap)
        t245 = dist(d, R245, GUAN, a, WALKER_A)
        # nearest-neighbour check for R162 -> any other subunit's Walker-A (register/scramble)
        cand = [(dist(d, R162, GUAN, x, WALKER_A), x) for x in labels if x != d]
        cand = [(v, x) for v, x in cand if v is not None]
        near = min(cand) if cand else (None, None)
        flag = "" if near[1] == a else f"{near[1]}@{near[0]:.1f}!" if near[0] is not None else ""
        trans162 += (t162 is not None and t162 < ENGAGED)
        cis162 += (c162 is not None and c162 < ENGAGED)
        trans245 += (t245 is not None and t245 < ENGAGED)
        f = lambda v: f"{v:.1f}" if v is not None else "  NA"
        print(f"{d:>4} {a:>6} | {f(t162):>9} {str(t162 is not None and t162<ENGAGED):>4} | "
              f"{f(c162):>8} | {f(t245):>9} | {near[1] or '-':>10}{'':>2}{flag}")
    print(f"# M2 (trans R162->neighbour WalkerA, <{ENGAGED}A): {trans162}/{n} engaged  (floor {n-1}/{n})")
    print(f"#    cis  R162->own WalkerA engaged: {cis162}/{n} ;  alt trans R245: {trans245}/{n}")

    # M1 geometry
    cents = [ca_centroid(subs[d]) for d in labels]
    if all(c is not None for c in cents):
        g = ring_geometry(cents, designed_order=True)
        sc = g["sequential_consistent"]
        sc_txt = "YES (designed==spatial order)" if sc else "NO -> SCRAMBLED register"
        print(f"# M1: radius {g['radius']:.1f} A | radius_CV {g['radius_cv']:.2f} | "
              f"planarity_rms {g['planarity_rms']:.1f} A | compact_ring {g['compact_ring']} | "
              f"sequential_consistent {sc_txt}")

if __name__ == "__main__":
    main()
