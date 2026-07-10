#!/usr/bin/env python3
"""Single-trajectory MM-GBSA interface energy for gp16 rings (design vs native).

Orthogonal, NON-deep-learning validation signal for the cp233 single-chain design.
For each adjacent subunit pair (the M2 coupler interface: donor R146 -> acceptor
Walker-A) compute the interface binding free energy

    dG_bind = G(complex) - G(subunitA) - G(subunitB)

with amber14 (ff14SB) + OBC2 implicit GBSA, in the single-trajectory approximation:
the complex is protonated and energy-minimized ONCE, and the two isolated subunits are
scored at those SAME coordinates (bonded + intramolecular terms cancel exactly, so the
result is the intermolecular vdW+elec plus the GB/SA desolvation of forming the interface).

Decomposition (all = complex - A - B):
    dE_MM   = intermolecular van der Waals + Coulomb (NonbondedForce, NoCutoff)
    dG_solv = GB polar + ACE nonpolar desolvation (OBC2 CustomGBForce)
    dG_bind = dE_MM + dG_solv

Subunits:
  native ring  -> separate chains (A_apo.pdb, chains A-E), adjacency by R146->Walker-A nearest neighbour
  design ring  -> copy ranges within one chain (C_design.cif), designed cyclic order k->k+1;
                  inter-copy linkers are dropped and each copy is treated as an independent
                  (zwitterionic-terminus) chain, exactly as the native chains are.

Usage:
  python mmgbsa_interface.py --native <A_apo.pdb> --out <dir> [--smoke]
  python mmgbsa_interface.py --design <C_design.cif> --copies "A:1-342,..." \
         --r146_incopy 255 --walker_incopy 133-140 --out <dir> [--smoke]
  (run once per structure; writes <label>_interface_energy.csv)
"""
import os, sys, argparse, json, csv, math, time
import numpy as np
from openmm import app, unit, CustomGBForce, NonbondedForce, Platform, LangevinMiddleIntegrator
from openmm.app import PDBFile, PDBxFile, ForceField, Modeller, NoCutoff, CutoffNonPeriodic, HBonds
from pdbfixer import PDBFixer

kcal = unit.kilocalorie_per_mole
KJ2KCAL = 1.0 / 4.184

GUAN = {"NE", "CZ", "NH1", "NH2"}
WALKER_OFFSETS = range(23, 31)      # gp16 res 24..31 (0-based offsets from copy start)
R146_OFFSET = 145
ENGAGED_A = 8.0

# ---------------------------------------------------------------- structure parsing
def parse_atoms(path):
    """Ordered list of dict(chain,resnum,resname,atom,elem,xyz). PDB or mmCIF _atom_site."""
    out = []
    if path.lower().endswith((".cif", ".mmcif")):
        order, cols, in_loop = [], {}, False
        for L in open(path):
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
                    rname = f[cols["label_comp_id"]]
                    an = f[cols["label_atom_id"]].strip('"')
                    el = f[cols["type_symbol"]]
                    xyz = (float(f[cols["Cartn_x"]]), float(f[cols["Cartn_y"]]), float(f[cols["Cartn_z"]]))
                except (ValueError, KeyError, IndexError):
                    continue
                out.append(dict(chain=ch, resnum=rn, resname=rname, atom=an, elem=el, xyz=xyz))
            elif in_loop and s and not s.startswith(("ATOM", "HETATM", "_")) and s == "#":
                in_loop = False
    else:
        for L in open(path):
            if L[:4] != "ATOM":
                continue
            an = L[12:16].strip()
            el = L[76:78].strip() or "".join(c for c in an if c.isalpha())[:1]
            out.append(dict(chain=L[21], resnum=int(L[22:26]), resname=L[17:20].strip(),
                            atom=an, elem=el,
                            xyz=(float(L[30:38]), float(L[38:46]), float(L[46:54]))))
    return out


def group_residues(atoms):
    """[(resnum,resname,[atomdict...])] preserving first-seen residue order."""
    res, seen = [], {}
    for a in atoms:
        key = a["resnum"]
        if key not in seen:
            seen[key] = len(res); res.append((a["resnum"], a["resname"], []))
        res[seen[key]][2].append(a)
    return res


def subunits_native(atoms):
    """chain -> ordered residue list."""
    by = {}
    for a in atoms:
        by.setdefault(a["chain"], []).append(a)
    return {c: group_residues(v) for c, v in by.items()}


def subunits_design(atoms, copies):
    """copy label -> ordered residue list (gp16-renumbered within copy: 1..N)."""
    subs = {}
    for label, ch, lo, hi in copies:
        sel = [a for a in atoms if a["chain"] == ch and lo <= a["resnum"] <= hi]
        # renumber to 1..N within copy so termini/numbering are handled like native
        res = group_residues(sel)
        subs[label] = [(rn - lo + 1, rname, ats) for (rn, rname, ats) in res]
    return subs


# ---------------------------------------------------------------- adjacency (M2 seam)
def guan_xyz(reslist, r146_res):
    for rn, rname, ats in reslist:
        if rn == r146_res:
            return np.array([a["xyz"] for a in ats if a["atom"] in GUAN and not a["atom"].startswith("H")])
    return np.empty((0, 3))


def walker_xyz(reslist, walker_res):
    P = []
    for rn, rname, ats in reslist:
        if rn in walker_res:
            P += [a["xyz"] for a in ats if not a["atom"].startswith("H")]
    return np.array(P) if P else np.empty((0, 3))


def min_dist(A, B):
    if len(A) == 0 or len(B) == 0:
        return math.inf
    return float(np.sqrt(((A[:, None, :] - B[None, :, :]) ** 2).sum(-1)).min())


def m2_adjacency(subs, r146_res, walker_res, designed_order=None):
    """Return ordered list of (donor,acceptor,R146_min_A). Designed order forces k->k+1."""
    labels = list(subs)
    pairs = []
    if designed_order:
        labels = designed_order
        n = len(labels)
        for i, d in enumerate(labels):
            a = labels[(i + 1) % n]
            dist = min_dist(guan_xyz(subs[d], r146_res), walker_xyz(subs[a], walker_res))
            pairs.append((d, a, dist))
    else:
        for d in labels:
            cand = [(min_dist(guan_xyz(subs[d], r146_res), walker_xyz(subs[x], walker_res)), x)
                    for x in labels if x != d]
            dist, a = min(cand)
            pairs.append((d, a, dist))
    return pairs


# ---------------------------------------------------------------- PDB writer
def write_pair_pdb(path, subA, subB):
    """Two subunits -> clean 2-chain PDB (chain A / chain B), residues renumbered 1..N each.
    Only heavy atoms (drop any H present); PDBFixer re-adds H consistently for both systems."""
    def emit(fh, reslist, chain, serial):
        rr = 0
        for (_, rname, ats) in reslist:
            rr += 1
            for a in ats:
                if a["atom"].startswith("H") or a["elem"] == "H":
                    continue
                x, y, z = a["xyz"]
                name = a["atom"]
                nm = (" " + name) if len(name) < 4 else name
                fh.write(f"ATOM  {serial:>5} {nm:<4} {rname:>3} {chain}{rr:>4}    "
                         f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00          {a['elem']:>2}\n")
                serial += 1
        return serial
    with open(path, "w") as fh:
        s = emit(fh, subA, "A", 1)
        fh.write("TER\n")
        emit(fh, subB, "B", s)
        fh.write("END\n")


# ---------------------------------------------------------------- energetics
# Minimize with a 2.0 nm cutoff (the project's validated OpenMM-MD protocol: 7x faster,
# geometry-preserving; interface contacts are all << 2 nm). Score the 3 final single points
# at NoCutoff (only 3 evals/interface) for an accurate direct-Coulomb + full-GB interface energy.
MIN_CUTOFF_NM = 2.0

def build_system(topology, ff, cutoff_nm=None):
    if cutoff_nm:
        return ff.createSystem(topology, nonbondedMethod=CutoffNonPeriodic,
                               nonbondedCutoff=cutoff_nm * unit.nanometer,
                               constraints=None, rigidWater=False)
    return ff.createSystem(topology, nonbondedMethod=NoCutoff, constraints=None, rigidWater=False)


def tag_force_groups(system):
    """NonbondedForce -> group 1 (MM vdW+elec); CustomGBForce -> group 2 (GB+SA)."""
    for f in system.getForces():
        if isinstance(f, NonbondedForce):
            f.setForceGroup(1)
        elif isinstance(f, CustomGBForce):
            f.setForceGroup(2)
        else:
            f.setForceGroup(0)


def energies(topology, positions, ff, platform):
    """Return dict of MM (group1) and GB (group2) potential energies in kcal/mol."""
    system = build_system(topology, ff)
    tag_force_groups(system)
    integ = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.001 * unit.picoseconds)
    ctx = app.Simulation(topology, system, integ, platform).context
    ctx.setPositions(positions)
    mm = ctx.getState(getEnergy=True, groups={1}).getPotentialEnergy().value_in_unit(kcal)
    gb = ctx.getState(getEnergy=True, groups={2}).getPotentialEnergy().value_in_unit(kcal)
    del ctx, integ, system
    return dict(MM=mm, GB=gb)


def split_modeller(modeller, keep_chain_index):
    """Copy modeller, delete all chains except the one with keep_chain_index -> (topology, positions)."""
    m = Modeller(modeller.topology, modeller.positions)
    chains = list(m.topology.chains())
    to_del = [c for i, c in enumerate(chains) if i != keep_chain_index]
    m.delete(to_del)
    return m.topology, m.positions


def mmgbsa_pair(pair_pdb, ff, platform, min_iters, min_tol_kj):
    """Prep + minimize the 2-chain complex; single-trajectory MM-GBSA. Returns dict of dG terms."""
    fixer = PDBFixer(filename=pair_pdb)
    fixer.findMissingResidues(); fixer.missingResidues = {}      # do NOT model gaps between chains
    fixer.findMissingAtoms(); fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    modeller = Modeller(fixer.topology, fixer.positions)

    system = build_system(modeller.topology, ff, cutoff_nm=MIN_CUTOFF_NM)
    integ = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    sim = app.Simulation(modeller.topology, system, integ, platform)
    sim.context.setPositions(modeller.positions)
    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kcal)
    sim.minimizeEnergy(tolerance=min_tol_kj * unit.kilojoule_per_mole / unit.nanometer, maxIterations=min_iters)
    st = sim.context.getState(getPositions=True, getEnergy=True)
    e1 = st.getPotentialEnergy().value_in_unit(kcal)
    minpos = st.getPositions()
    modeller = Modeller(modeller.topology, minpos)             # minimized complex coords
    del sim, integ, system

    n_chains = modeller.topology.getNumChains()
    assert n_chains == 2, f"expected 2 chains, got {n_chains}"
    e_cx = energies(modeller.topology, modeller.positions, ff, platform)
    topA, posA = split_modeller(modeller, 0)
    topB, posB = split_modeller(modeller, 1)
    e_A = energies(topA, posA, ff, platform)
    e_B = energies(topB, posB, ff, platform)

    dMM = e_cx["MM"] - e_A["MM"] - e_B["MM"]
    dGB = e_cx["GB"] - e_A["GB"] - e_B["GB"]
    return dict(dG_bind=dMM + dGB, dE_MM=dMM, dG_solv=dGB,
                E_min_start=e0, E_min_final=e1,
                nres_A=topA.getNumResidues(), nres_B=topB.getNumResidues(),
                natom_cx=modeller.topology.getNumAtoms())


# ---------------------------------------------------------------- driver
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--native", help="native ring PDB (chains = subunits)")
    ap.add_argument("--design", help="design ring mmCIF (copies within one chain)")
    ap.add_argument("--copies", default="", help='e.g. "A:1-342,A:353-694,..."')
    ap.add_argument("--r146_incopy", type=int, default=None)
    ap.add_argument("--walker_incopy", default=None)
    ap.add_argument("--label", default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min_iters", type=int, default=2000)
    ap.add_argument("--min_tol_kj", type=float, default=10.0)
    ap.add_argument("--platform", default="CPU")
    ap.add_argument("--smoke", action="store_true", help="only the first interface")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    ff = ForceField('amber14-all.xml', 'implicit/obc2.xml')
    try:
        platform = Platform.getPlatformByName(args.platform)
    except Exception:
        platform = Platform.getPlatformByName("Reference")
    if args.platform == "CPU":
        try: platform.setPropertyDefaultValue("Threads", str(os.cpu_count() or 4))
        except Exception: pass

    if args.native:
        path, label = args.native, args.label or "native"
        atoms = parse_atoms(path)
        subs = subunits_native(atoms)
        walker_res = [o + 1 for o in WALKER_OFFSETS]; r146_res = R146_OFFSET + 1
        pairs = m2_adjacency(subs, r146_res, walker_res, designed_order=None)
        get = lambda lbl: subs[lbl]
    else:
        path, label = args.design, args.label or "design"
        copies = []
        for i, tok in enumerate(args.copies.split(",")):
            ch, span = tok.split(":"); lo, hi = span.split("-")
            copies.append((f"{ch}{i+1}", ch, int(lo), int(hi)))
        atoms = parse_atoms(path)
        subs = subunits_design(atoms, copies)
        r146_res = args.r146_incopy; lo, hi = (int(x) for x in args.walker_incopy.split("-"))
        walker_res = list(range(lo, hi + 1))
        order = [c[0] for c in copies]
        pairs = m2_adjacency(subs, r146_res, walker_res, designed_order=order)
        get = lambda lbl: subs[lbl]

    if args.smoke:
        pairs = pairs[:1]

    rows = []
    print(f"# {label}: {path}")
    print(f"# {'donor':>6} {'acc':>6} {'R146_A':>7} {'dG_bind':>9} {'dE_MM':>9} {'dG_solv':>9}  (kcal/mol)")
    for d, a, r146d in pairs:
        pdb = os.path.join(args.out, f"_pair_{label}_{d}_{a}.pdb")
        write_pair_pdb(pdb, get(d), get(a))
        t0 = time.time()
        res = mmgbsa_pair(pdb, ff, platform, args.min_iters, args.min_tol_kj)
        res.update(system=label, donor=str(d), acceptor=str(a), R146_min_A=round(r146d, 2),
                   seconds=round(time.time() - t0, 1))
        rows.append(res)
        print(f"  {str(d):>6} {str(a):>6} {r146d:>7.2f} {res['dG_bind']:>9.1f} "
              f"{res['dE_MM']:>9.1f} {res['dG_solv']:>9.1f}  ({res['seconds']}s)")
        try: os.remove(pdb)
        except OSError: pass

    dg = [r["dG_bind"] for r in rows]
    print(f"# {label} mean dG_bind = {np.mean(dg):.1f} +/- {np.std(dg):.1f} kcal/mol  (n={len(dg)})")

    csv_path = os.path.join(args.out, f"{label}_interface_energy.csv")
    cols = ["system", "donor", "acceptor", "R146_min_A", "dG_bind", "dE_MM", "dG_solv",
            "E_min_start", "E_min_final", "nres_A", "nres_B", "natom_cx", "seconds"]
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in cols})
    json.dump(dict(system=label, structure=path, mean_dG_bind=float(np.mean(dg)),
                   std_dG_bind=float(np.std(dg)), rows=rows),
              open(os.path.join(args.out, f"{label}_interface_energy.json"), "w"), indent=2)
    print(f"# wrote {csv_path}")


if __name__ == "__main__":
    main()
