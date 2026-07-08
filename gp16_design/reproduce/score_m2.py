#!/usr/bin/env python3
"""M2 scorer — per-interface trans-R146 arginine-finger engagement.

The load-bearing design metric (see gp16_design/SUBMISSION.md M2, DESIGN_NOTES §14):
for each subunit's R146 guanidinium (NE/CZ/NH1/NH2), the minimum heavy-atom
distance to a *neighbouring* subunit's Walker-A (res 24-31). The neighbour is the
argmin acceptor. "engaged" = < 8 Angstrom (native 7JQQ ~3.2; predicted apo ring ~6.6).

Atom selection is VALIDATED: on the local native-ring Boltz seed-1 structure this
reproduces outputs/cycle3_ring_scores.csv to 0.01 A
(A->B 6.57, B->C 6.62, C->D 6.53, D->E 6.64, E->A 6.69).

Subunits come either from separate chains (native ring, B1's 3 WT chains) or from
copy ranges within one chain (covalent constructs) via --copies. Numbering inside a
copy follows the native gp16 numbering, so R146 = copy_start+145, Walker-A = copy_start+23..30.
Also reports M4 = per-interface mean pLDDT (B-factor col) over the R146+Walker-A atoms.

Usage:
  python score_m2.py native_ring.pdb
  python score_m2.py pentamer.cif --copies "A:1-367,A:368-734,A:735-1101,A:1102-1468,A:1469-1795"
       (copy ranges are residue-number spans within the named chain; derive from the
        construct manifest / linker layout. Res 146/24-31 are taken relative to each span start.)
"""
import sys, argparse, math
import numpy as np

GUAN = {"NE", "CZ", "NH1", "NH2"}
WALKER_OFFSETS = range(23, 31)      # res 24..31 as 0-based offsets from a copy start
R146_OFFSET = 145                   # res 146
ENGAGED_A = 8.0


def parse_atoms(path):
    """Return list of (chain, resnum, atomname, xyz(np), bfactor). PDB or mmCIF atom_site."""
    out = []
    if path.lower().endswith((".cif", ".mmcif")):
        cols, in_loop, order = {}, False, []
        with open(path) as fh:
            for L in fh:
                s = L.strip()
                if s.startswith("_atom_site."):
                    order.append(s.split(".", 1)[1]); cols = {k: i for i, k in enumerate(order)}
                    in_loop = True; continue
                if in_loop and (s.startswith("ATOM") or s.startswith("HETATM")):
                    f = s.split()
                    if len(f) < len(order):
                        continue
                    try:
                        ch = f[cols.get("auth_asym_id", cols.get("label_asym_id"))]
                        rn = int(f[cols.get("auth_seq_id", cols.get("label_seq_id"))])
                        an = f[cols["label_atom_id"]].strip('"')
                        xyz = np.array([float(f[cols["Cartn_x"]]), float(f[cols["Cartn_y"]]), float(f[cols["Cartn_z"]])])
                        b = float(f[cols.get("B_iso_or_equiv", cols.get("pLDDT", -1))]) if "B_iso_or_equiv" in cols else 0.0
                    except (ValueError, KeyError, IndexError):
                        continue
                    out.append((ch, rn, an, xyz, b))
                elif in_loop and s and not s.startswith(("ATOM", "HETATM", "_")):
                    if s == "#":
                        in_loop = False
    else:
        for L in open(path):
            if L[:4] != "ATOM":
                continue
            out.append((L[21], int(L[22:26]), L[12:16].strip(),
                        np.array([float(L[30:38]), float(L[38:46]), float(L[46:54])]),
                        float(L[60:66])))
    return out


def subunit_index(atoms, copies, copy_start_res=1):
    """Map subunit label -> {gp16-residue-number: [(atomname,xyz,b)]}.

    Chain case: keys are the deposited residue numbers (native gp16 numbering).
    Copy-range case: the copy's first construct residue corresponds to gp16 residue
    `copy_start_res` (e.g. 4 for a res4-330 core), so the key is remapped to gp16 numbering.
    """
    subs = {}
    if copies:
        for label, ch, lo, hi in copies:
            d = {}
            for c, rn, an, xyz, b in atoms:
                if c == ch and lo <= rn <= hi:
                    d.setdefault(rn - lo + copy_start_res, []).append((an, xyz, b))
            subs[label] = d
    else:
        for c, rn, an, xyz, b in atoms:
            subs.setdefault(c, {}).setdefault(rn, []).append((an, xyz, b))
    return subs


def atoms_for(sub, resnums, names=None):
    P, B = [], []
    for rn in resnums:
        for an, xyz, b in sub.get(rn, []):
            if an.startswith("H"):
                continue
            if names and an not in names:
                continue
            P.append(xyz); B.append(b)
    return (np.array(P) if P else np.empty((0, 3))), B


def min_dist(A, B):
    if len(A) == 0 or len(B) == 0:
        return None
    return float(np.sqrt(((A[:, None, :] - B[None, :, :]) ** 2).sum(-1)).min())


def ca_centroid(sub):
    P = [xyz for atomlist in sub.values() for an, xyz, b in atomlist if an == "CA"]
    return np.mean(P, axis=0) if P else None


def ring_geometry(centroids, designed_order):
    """M1 — ring-closure geometry from the per-subunit CA centroids (in designed order).

    Returns radius, radius_CV (spread of centroid-to-centre distances; ~0 = on a circle),
    planarity_rms (out-of-plane RMS), compact_ring (centroids lie on a circle), and
    sequential_consistent (the angular order around the ring == the designed cyclic order,
    i.e. designed neighbours are also spatial neighbours — False = a scrambled register).
    """
    C = np.array(centroids); k = len(C)
    center = C.mean(axis=0); X = C - center
    _, _, Vt = np.linalg.svd(X)
    normal = Vt[2]
    oop = X @ normal
    planarity_rms = float(np.sqrt((oop ** 2).mean()))
    inplane = X - np.outer(oop, normal)
    radii = np.linalg.norm(inplane, axis=1)
    radius = float(radii.mean())
    radius_cv = float(radii.std() / radii.mean()) if radii.mean() > 0 else float("nan")
    ang = np.arctan2(inplane @ Vt[1], inplane @ Vt[0])
    order = list(np.argsort(ang))                       # subunit indices sorted by angle
    seq_consistent = None
    if designed_order:                                  # only meaningful when copies are in designed order
        diffs = [(order[(i + 1) % k] - order[i]) % k for i in range(k)]
        seq_consistent = all(d in (1, k - 1) for d in diffs)
    compact = radius_cv < 0.35
    return dict(radius=radius, radius_cv=radius_cv, planarity_rms=planarity_rms,
                compact_ring=compact, sequential_consistent=seq_consistent)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("structure")
    ap.add_argument("--copies", default="", help='e.g. "A:1-327,A:368-694,..." residue spans within a chain')
    ap.add_argument("--copy_start_res", type=int, default=1,
                    help="gp16 residue number of each copy's first residue (4 for a res4-330 core)")
    args = ap.parse_args()

    copies = None
    if args.copies:
        copies = []
        for i, tok in enumerate(args.copies.split(",")):
            ch, span = tok.split(":"); lo, hi = span.split("-")
            copies.append((f"{ch}{i+1}", ch, int(lo), int(hi)))

    atoms = parse_atoms(args.structure)
    subs = subunit_index(atoms, copies, args.copy_start_res)
    labels = list(subs)

    walker_res = [o + 1 for o in WALKER_OFFSETS]        # 24..31 (1-based within copy)
    r146_res = R146_OFFSET + 1                           # 146
    n = len(labels)

    def r146_to(donor, acceptor):
        r_gua, _ = atoms_for(subs[donor], [r146_res], GUAN)
        _, r146_b = atoms_for(subs[donor], [r146_res])
        wa, wa_b = atoms_for(subs[acceptor], walker_res)
        dist = min_dist(r_gua, wa)
        plddt = np.mean(r146_b + wa_b) if (r146_b or wa_b) else float("nan")
        return dist, plddt

    print(f"# {args.structure}   subunits: {labels}")
    print(f"{'donor':>7} {'->':>2} {'acceptor':>8} {'R146_min_A':>11} {'engaged':>8} {'iface_pLDDT':>12}  note")
    engaged_n = 0
    iface_plddts = []
    for i, d in enumerate(labels):
        if copies:
            # DESIGNED sequential ring neighbour (covalent order), cyclic donor -> next.
            a = labels[(i + 1) % n]
            dist, plddt = r146_to(d, a)
            # diagnostic: is R146 actually nearest to a DIFFERENT subunit? (topology scramble)
            nearest = min(((r146_to(d, x)[0], x) for x in labels if x != d and r146_to(d, x)[0] is not None),
                          default=(None, None))
            note = "" if nearest[1] == a else f"nearest={nearest[1]}@{nearest[0]:.1f}(scrambled)"
        else:
            # separate chains (symmetric ring, no designed order): nearest neighbour defines adjacency.
            cand = [(r146_to(d, x)[0], x) for x in labels if x != d and r146_to(d, x)[0] is not None]
            dist, a = min(cand)
            _, plddt = r146_to(d, a); note = ""
        if dist is None:
            continue
        eng = dist < ENGAGED_A
        engaged_n += eng
        if not (np.isnan(plddt)):
            iface_plddts.append(plddt)
        print(f"{d:>7} {'->':>2} {a:>8} {dist:>11.2f} {str(eng):>8} {plddt:>12.1f}  {note}")
    print(f"# M2: {engaged_n}/{n} DESIGNED interfaces engaged (<{ENGAGED_A} A)  ->  pass = >=4/5 for a pentamer")

    # M1 — ring-closure geometry
    centroids = [ca_centroid(subs[d]) for d in labels]
    if all(c is not None for c in centroids):
        g = ring_geometry(centroids, designed_order=bool(copies))
        sc = g["sequential_consistent"]
        sc_txt = "n/a (needs designed copies)" if sc is None else (
            "YES (designed order == spatial order)" if sc else "NO -> SCRAMBLED register")
        print(f"# M1: radius {g['radius']:.1f} A | radius_CV {g['radius_cv']:.2f} | planarity_rms "
              f"{g['planarity_rms']:.1f} A | compact_ring {g['compact_ring']} | sequential_consistent: {sc_txt}")
    # M4 — per-interface confidence roll-up
    if iface_plddts:
        print(f"# M4: designed-interface pLDDT mean {np.mean(iface_plddts):.1f} / min {np.min(iface_plddts):.1f}")


if __name__ == "__main__":
    main()
