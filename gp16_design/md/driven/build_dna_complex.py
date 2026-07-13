#!/usr/bin/env python3
"""Track A / step S0 of md/STEERED_TARGETED_MD_PLAN.md — thread 7JQQ dsDNA through the
design ring channel to build the SMD (in-silico optical-tweezers) starting complex.

Method (no per-residue correspondence needed): compute the ring FRAME (center, axis,
azimuth) of the design ring and of the 7JQQ protein ring from their 5 subunit CA
centroids, then apply the rigid transform that maps the 7JQQ frame onto the design
frame to the 7JQQ dsDNA (chains F,G). The DNA lands along the design channel axis.

Outputs (md/driven/inputs/):
  C_plus_dna.pdb   design protein (from C_start.pdb) + threaded dsDNA chains F,G
  build_report.json  min protein-DNA distance (clash check), DNA COM offset, axis info
"""
import argparse, os, json
import numpy as np
import mdtraj as md

COPY_LOS = [1, 353, 705, 1057, 1409]
COPY_LEN = 342
KEEP_DNA = {'F', 'G'}   # 7JQQ dsDNA


def ring_frame(centroids):
    """centroids: (5,3) A. Return (center, axis, u, w) orthonormal frame."""
    C = np.asarray(centroids, float)
    ctr = C.mean(0)
    X = C - ctr
    _, _, Vt = np.linalg.svd(X)
    axis = Vt[2]
    if axis[2] < 0:            # consistent orientation (point +z)
        axis = -axis
    u = X[0] - (X[0] @ axis) * axis     # in-plane direction to subunit 0
    u = u / np.linalg.norm(u)
    w = np.cross(axis, u)
    return ctr, axis, u, w


def design_centroids(pdb):
    t = md.load(pdb); xyz = t.xyz[0] * 10.0
    ch = list(t.topology.chains)[0]
    resmap = {r.resSeq: r for r in ch.residues}
    cents = []
    for lo in COPY_LOS:
        hi = lo + COPY_LEN - 1
        cas = [a.index for rn in range(lo, hi + 1) if rn in resmap
               for a in resmap[rn].atoms if a.name == 'CA']
        cents.append(xyz[cas].mean(0))
    return np.array(cents), t


def parse_cif_atoms(cif):
    """Return list of dicts for ATOM records: chain, resn, resi, name, elem, xyz(A)."""
    order, cols, inloop, out = [], {}, False, []
    for L in open(cif):
        st = L.strip()
        if st.startswith('_atom_site.'):
            order.append(st.split('.', 1)[1]); cols = {k: i for i, k in enumerate(order)}
            inloop = True; continue
        if inloop and (st.startswith('ATOM') or st.startswith('HETATM')):
            f = st.split()
            if len(f) < len(order):
                continue
            def g(*names):
                for n in names:
                    if n in cols:
                        return f[cols[n]]
                return None
            out.append(dict(
                group=f[cols['group_PDB']],
                chain=g('auth_asym_id', 'label_asym_id'),
                resn=f[cols['label_comp_id']],
                resi=int(g('auth_seq_id', 'label_seq_id')),
                name=f[cols['label_atom_id']].strip('"'),
                elem=(g('type_symbol') or f[cols['label_atom_id']][0]).strip(),
                alt=g('label_alt_id') or '.',
                xyz=np.array([float(f[cols['Cartn_x']]), float(f[cols['Cartn_y']]),
                              float(f[cols['Cartn_z']])]))
            )
        elif inloop and st == '#':
            inloop = False
    return out


def jqq_protein_centroids(atoms):
    """5 chain CA centroids (chains A-E) for the 7JQQ ring frame."""
    chains = {}
    for a in atoms:
        if a['group'] == 'ATOM' and a['chain'] in set('ABCDE') and a['name'] == 'CA' \
                and a['alt'] in ('.', 'A'):
            chains.setdefault(a['chain'], []).append(a['xyz'])
    cents = [np.mean(chains[c], 0) for c in sorted(chains)]
    return np.array(cents)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--design', default='../openmm_validation/trajectories/C/C_start.pdb')
    ap.add_argument('--cif', default='../../data/raw/7JQQ.cif')
    ap.add_argument('--out', default='inputs/C_plus_dna.pdb')
    ap.add_argument('--report', default='inputs/build_report.json')
    args = ap.parse_args()

    d_cents, d_traj = design_centroids(args.design)
    d_ctr, d_axis, d_u, d_w = ring_frame(d_cents)
    atoms = parse_cif_atoms(args.cif)
    q_cents = jqq_protein_centroids(atoms)
    q_ctr, q_axis, q_u, q_w = ring_frame(q_cents)

    # rigid transform: q-frame -> d-frame.  p' = R (p - q_ctr) + d_ctr
    Fd = np.stack([d_u, d_w, d_axis], 0)   # rows = d-frame basis (world->d coords)
    Fq = np.stack([q_u, q_w, q_axis], 0)
    R = Fd.T @ Fq                           # maps q-world vec -> d-world vec

    # transform DNA (chains F,G) into the design frame
    dna = [a for a in atoms if a['chain'] in KEEP_DNA and a['alt'] in ('.', 'A')]
    dna_xyz0 = np.array([R @ (a['xyz'] - q_ctr) + d_ctr for a in dna])

    # protein heavy atoms (subsampled for the scan)
    prot_xyz = d_traj.xyz[0] * 10.0
    prot_heavy = np.array([prot_xyz[a.index] for a in d_traj.topology.atoms
                           if a.element is not None and a.element.symbol != 'H'])
    sub = prot_heavy[::4]

    def rot_about(vecs, ax, theta):
        c, s = np.cos(theta), np.sin(theta)
        return vecs * c + np.cross(ax, vecs) * s + np.outer(vecs @ ax, ax) * (1 - c)

    # 1) align the DNA helical (principal) axis to the ring axis
    com0 = dna_xyz0.mean(0)
    Xd = dna_xyz0 - com0
    _, _, Vt = np.linalg.svd(Xd)
    dna_axis = Vt[0]
    if dna_axis @ d_axis < 0:
        dna_axis = -dna_axis
    v = np.cross(dna_axis, d_axis); s = np.linalg.norm(v); cth = dna_axis @ d_axis
    if s > 1e-6:
        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        Ralign = np.eye(3) + vx + vx @ vx * ((1 - cth) / (s * s))
    else:
        Ralign = np.eye(3)
    # 2) put DNA COM exactly on the ring axis at the ring center (radial+axial centered)
    dna_aligned = (Ralign @ Xd.T).T + d_ctr

    def score(xyz):
        com = xyz.mean(0)
        axial = (com - d_ctr) @ d_axis
        radial = np.linalg.norm((com - d_ctr) - axial * d_axis)
        dmin = np.min(np.linalg.norm(sub[:, None, :] - xyz[None, :, :], axis=2))
        return dmin, axial, radial

    # 3) scan azimuth about the ring axis x small axial slide; maximize clearance
    best = None
    for theta in np.linspace(0, 2 * np.pi, 72, endpoint=False):
        rotated = rot_about(dna_aligned - d_ctr, d_axis, theta) + d_ctr
        for ds in np.linspace(-8, 8, 9):
            xyz = rotated + d_axis * ds
            dmin, axial, radial = score(xyz)
            if best is None or dmin > best[0]:
                best = (dmin, theta, ds, axial, radial, xyz)
    dmin, theta, slide, axial_off, radial_off, dna_xyz = best
    for a, p in zip(dna, dna_xyz):
        a['xyz'] = p
    # TRUE clearance on ALL protein heavy atoms (the scan used prot_heavy[::4], which can miss the
    # worst clash by ~1.4 A -- that subsample once masked a 0.427 A overlap as "1.79 A").
    dmin_full = float(np.min(np.linalg.norm(prot_heavy[:, None, :] - dna_xyz[None, :, :], axis=2)))
    usable = dmin_full >= 2.2
    print(f"[place] axis-aligned; best azimuth={np.degrees(theta):.0f}deg slide={slide:+.1f}A -> "
          f"scan-clearance={dmin:.2f}A  FULL-ATOM clearance={dmin_full:.2f}A  "
          f"axial_off={axial_off:.2f}A radial_off={radial_off:.2f}A")
    if not usable:
        print(f"[place] *** WARNING: full-atom clearance {dmin_full:.2f} A < 2.2 A = STERIC OVERLAP. "
              f"Rigid threading of 7JQQ DNA into the (differently-shaped) design pore cannot yield a "
              f"usable S0. The complex MUST be relaxed (restrained energy minimization / soft-core "
              f"insertion) before any SMD. See relax_complex.py. ***")

    # ---- write combined PDB: design protein (chain A) + DNA (chains F,G) ----
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    lines, serial = [], 0
    # design protein via mdtraj topology + xyz
    top = d_traj.topology
    for atom in top.atoms:
        serial += 1
        r = atom.residue
        x, y, z = prot_xyz[atom.index]
        nm = atom.name
        nm4 = nm if len(nm) == 4 else ' ' + nm
        elem = atom.element.symbol if atom.element is not None else nm[0]
        lines.append(f"ATOM  {serial:>5d} {nm4:<4s}{'':1s}{r.name:>3s} A{r.resSeq:>4d}{'':1s}   "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {elem:>2s}")
    lines.append('TER')
    prev = None
    STRIP_5P = {'P', 'OP1', 'OP2', 'OP3'}   # 5'-terminal dangling phosphate -> make 5'-OH (matches *5 template)
    chain_min = {}
    for a in dna:
        chain_min[a['chain']] = min(chain_min.get(a['chain'], a['resi']), a['resi'])
    for a in dna:
        if a['resi'] == chain_min[a['chain']] and a['name'] in STRIP_5P:
            continue   # drop 5'-phosphate on the 5'-terminal residue
        serial += 1
        if prev is not None and a['chain'] != prev:
            lines.append('TER')
        x, y, z = a['xyz']
        nm = a['name']; nm4 = nm if len(nm) == 4 else ' ' + nm
        lines.append(f"ATOM  {serial:>5d} {nm4:<4s}{'':1s}{a['resn']:>3s} {a['chain']:1s}{a['resi']:>4d}{'':1s}   "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {a['elem']:>2s}")
        prev = a['chain']
    lines += ['TER', 'END']
    open(args.out, 'w').write('\n'.join(lines) + '\n')

    report = dict(
        design=args.design, cif=args.cif, out=args.out,
        n_dna_atoms=len(dna), n_protein_atoms=int(top.n_atoms),
        design_ring_center_A=[round(float(v), 2) for v in d_ctr],
        design_axis=[round(float(v), 3) for v in d_axis],
        dna_com_axial_offset_A=round(float(axial_off), 2),
        dna_com_radial_offset_A=round(float(radial_off), 2),
        scan_clearance_subsampled_A=round(float(dmin), 2),
        min_protein_dna_dist_full_atom_A=round(float(dmin_full), 2),
        usable_without_relaxation=bool(usable),
        note=("FULL-ATOM min protein-DNA distance is the real clearance. Rigid threading tops out "
              "near ~1 A (steric overlap) because the design pore differs from 7JQQ's -> the raw "
              "complex is a DRAFT that MUST be relaxed (relax_complex.py: restrained minimization) "
              "before SMD. usable_without_relaxation=false means: relax first."))
    json.dump(report, open(args.report, 'w'), indent=2)
    print(json.dumps(report, indent=2))
    print(f"\nwrote {args.out} ({serial} atoms) and {args.report}")


if __name__ == '__main__':
    main()
