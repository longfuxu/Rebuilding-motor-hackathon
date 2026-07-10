#!/usr/bin/env python3
"""Shared residue mapping + landmark selection for the gp16 coordination /
force-transmission analysis (cp233 single-chain design C vs native apo ring A).

Both structures are reduced to the SAME residue set for a fair comparison:
native gp16 residues 4..330 in every one of the 5 subunits (1635 CA each).

Design C (single chain, resSeq 1..1750) is a circular permutant at 233. Verified
segment layout (see build log): each subunit =
    native 234..330 (97 res) + (GGGGS)3 linker + native 4..233 (230 res) + (GGGGS)2 linker
so the 5 subunits are covalently ordered 0-1-2-3-4 along the chain. Linker residues
are DROPPED (they are not part of the folded ring and have no native counterpart).

Native A: 5 chains A-E, resSeq = native numbering directly.

Landmarks (native gp16 numbering):
  Walker-A / ATP P-loop : 24..31   (GARGIGKS) -> the ATP site of each subunit
  Y129                  : force-transmission residue (follow-on program Aim 1a)
  R146                  : trans-acting arginine finger
  DNA-contact residues  : from reproduce/functional_contacts_7jqq.json
"""
import numpy as np
import mdtraj as md

# ---- landmark residue numbers (native gp16 numbering) ----
WALKERA = list(range(24, 32))          # 24..31
Y129 = 129
R146 = 146
DNA_CONTACT = [55, 56, 57, 59, 60, 82, 83, 98, 99, 100,
               125, 126, 127, 128, 129, 130, 292, 293, 297, 330]
RES_LO, RES_HI = 4, 330                 # common comparable residue window

# ---- design subunit segment starts (design resSeq, 1-based) ----
# from segment alignment: PROT segA starts and segB starts per subunit
_SEGA_START = [1, 353, 705, 1057, 1409]     # design resSeq where native 234 sits
_SEGB_START = [113, 465, 817, 1169, 1521]   # design resSeq where native 4   sits


def design_resseq(subunit, native_res):
    """design resSeq (1-based) for a given subunit (0..4) and native residue number."""
    if 234 <= native_res <= 330:
        return _SEGA_START[subunit] + (native_res - 234)
    if 4 <= native_res <= 233:
        return _SEGB_START[subunit] + (native_res - 4)
    raise ValueError(f"native res {native_res} outside design coverage 4..330")


def load_ca(topfile):
    """Return an mdtraj Trajectory (topology/first frame) for CA selection helpers."""
    return md.load(topfile)


def ca_table_design(top):
    """For the design topology, return list of (row, subunit, native_res, atom_index)
    for every CA in native residues 4..330 across the 5 subunits, in
    (subunit, native_res) ascending order. Row = position in the reduced CA array."""
    top_ = top.topology
    # map design resSeq -> CA atom index
    res2ca = {}
    for a in top_.atoms:
        if a.name == 'CA' and a.residue.is_protein:
            res2ca[a.residue.resSeq] = a.index
    rows = []
    r = 0
    for k in range(5):
        for nres in range(RES_LO, RES_HI + 1):
            drs = design_resseq(k, nres)
            ai = res2ca.get(drs)
            if ai is None:
                raise RuntimeError(f"missing CA design resSeq {drs} (sub {k} nat {nres})")
            rows.append((r, k, nres, ai))
            r += 1
    return rows


def ca_table_native(top):
    """For the native topology (5 chains), return list of (row, subunit, native_res,
    atom_index) for CA in residues 4..330, chains 0..4, in (subunit, native_res) order."""
    top_ = top.topology
    chains = list(top_.chains)
    rows = []
    r = 0
    for k, c in enumerate(chains):
        res2ca = {}
        for res in c.residues:
            if res.is_protein:
                for a in res.atoms:
                    if a.name == 'CA':
                        res2ca[res.resSeq] = a.index
        for nres in range(RES_LO, RES_HI + 1):
            ai = res2ca.get(nres)
            if ai is None:
                raise RuntimeError(f"missing CA native chain {k} res {nres}")
            rows.append((r, k, nres, ai))
            r += 1
    return rows


def ca_table(top, kind):
    return ca_table_design(top) if kind == 'design' else ca_table_native(top)


def landmark_rows(rows):
    """Given a ca_table, return dicts of row indices for each landmark set, per subunit.
    Returns: walkerA[k]=list rows, y129[k]=row, r146[k]=row, dna[k]=list rows,
    plus per-row lookup (subunit,native_res)->row."""
    lut = {(k, nres): r for (r, k, nres, ai) in rows}
    walkerA = {k: [lut[(k, n)] for n in WALKERA] for k in range(5)}
    y129 = {k: lut[(k, Y129)] for k in range(5)}
    r146 = {k: lut[(k, R146)] for k in range(5)}
    dna = {k: [lut[(k, n)] for n in DNA_CONTACT] for k in range(5)}
    return dict(walkerA=walkerA, y129=y129, r146=r146, dna=dna, lut=lut)


def subunit_ids(rows):
    return np.array([k for (_, k, _, _) in rows], dtype=int)


def native_res_ids(rows):
    return np.array([n for (_, _, n, _) in rows], dtype=int)


def atom_indices(rows):
    return np.array([ai for (_, _, _, ai) in rows], dtype=int)


def ring_neighbor_order(top, rows, xyz=None):
    """Determine the spatial ring order of the 5 subunits from Walker-A centroids.
    Returns list of subunit ids in ring order and the neighbor map k -> next_k."""
    if xyz is None:
        xyz = top.xyz[0]
    lm = landmark_rows(rows)
    ai = atom_indices(rows)
    cents = []
    for k in range(5):
        idx = ai[lm['walkerA'][k]]
        cents.append(xyz[idx].mean(axis=0))
    cents = np.array(cents)
    c0 = cents.mean(axis=0)
    v = cents - c0
    # project onto best-fit plane (ring plane) via SVD, take angle
    u, s, vt = np.linalg.svd(v - v.mean(0))
    basis = vt[:2]
    ang = np.arctan2(v @ basis[1], v @ basis[0])
    order = list(np.argsort(ang))
    nxt = {order[i]: order[(i + 1) % 5] for i in range(5)}
    prv = {order[i]: order[(i - 1) % 5] for i in range(5)}
    return order, nxt, prv, cents


if __name__ == '__main__':
    for kind, top in [('design', 'gp16_design/md/openmm_validation/trajectories/C/C_start.pdb'),
                      ('native', 'gp16_design/md/openmm_validation/trajectories/A/A_start.pdb')]:
        t = load_ca(top)
        rows = ca_table(t, kind)
        print(f"{kind}: {len(rows)} CA in reduced set (expect 1635)")
        order, nxt, prv, cents = ring_neighbor_order(t, rows)
        print(f"  spatial ring order (by Walker-A centroid angle): {order}")
        print(f"  neighbor map (k -> k+1 around ring): {nxt}")
        lm = landmark_rows(rows)
        print(f"  subunit0 Y129 row={lm['y129'][0]}, R146 row={lm['r146'][0]}, "
              f"WalkerA rows={lm['walkerA'][0][:3]}...")
