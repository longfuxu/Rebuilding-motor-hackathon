#!/usr/bin/env python3
"""PyRosetta InterfaceAnalyzer interface energy for gp16 rings (design vs native).

Second, independent semi-physical signal (Rosetta ref2015 energy units, REU) — orthogonal
both to the deep-learning predictors AND to the OpenMM MM-GBSA (different energy function,
different solvation model). For each adjacent subunit pair (the M2 coupler interface) we:
  1. build the 2-chain complex (chain A = donor, chain B = acceptor; H/side-chains rebuilt),
  2. FastRelax with coordinate constraints to the model (relieve clashes, no big drift),
  3. InterfaceAnalyzerMover with pack_separated -> dG_separated (REU), buried SASA, n_interface_res.
Reports the per-interface dG_separated and the mean for design vs native.

Subunit definition + adjacency (M2 seam) are shared with the MM-GBSA script.

Usage: same --native / --design flags as mmgbsa_interface.py.
"""
import os, sys, argparse, json, csv, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from mmgbsa_interface import (parse_atoms, subunits_native, subunits_design,
                              m2_adjacency, write_pair_pdb, WALKER_OFFSETS, R146_OFFSET)

import pyrosetta
from pyrosetta import pose_from_pdb, get_fa_scorefxn
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.relax import FastRelax


def analyze_pair(pair_pdb, sfxn, relax_cycles):
    pose = pose_from_pdb(pair_pdb)
    if relax_cycles > 0:
        fr = FastRelax(sfxn, relax_cycles)
        fr.constrain_relax_to_start_coords(True)      # keep the modelled interface, just de-clash
        fr.apply(pose)
    ia = InterfaceAnalyzerMover("A_B")
    ia.set_scorefunction(sfxn)
    ia.set_compute_packstat(True)
    ia.set_pack_separated(True)                        # repack the separated state -> proper dG_separated
    ia.set_compute_interface_sc(True)
    ia.apply(pose)
    return dict(
        dG_separated=float(ia.get_separated_interface_energy()),
        dG_dsasa_ratio=float(ia.get_interface_dG() / ia.get_interface_delta_sasa()
                             if ia.get_interface_delta_sasa() else float("nan")),
        interface_dG=float(ia.get_interface_dG()),
        dSASA=float(ia.get_interface_delta_sasa()),
        n_interface_res=int(ia.get_num_interface_residues()),
        total_score=float(sfxn(pose)),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--native"); ap.add_argument("--design")
    ap.add_argument("--copies", default="")
    ap.add_argument("--r146_incopy", type=int, default=None)
    ap.add_argument("--walker_incopy", default=None)
    ap.add_argument("--label", default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--relax_cycles", type=int, default=1)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    pyrosetta.init("-mute all -ignore_unrecognized_res -ignore_zero_occupancy false "
                   "-detect_disulf false -no_optH false")
    sfxn = get_fa_scorefxn()   # ref2015

    if args.native:
        path, label = args.native, args.label or "native"
        subs = subunits_native(parse_atoms(path))
        walker_res = [o + 1 for o in WALKER_OFFSETS]; r146_res = R146_OFFSET + 1
        pairs = m2_adjacency(subs, r146_res, walker_res, designed_order=None)
        order = None
    else:
        path, label = args.design, args.label or "design"
        copies = []
        for i, tok in enumerate(args.copies.split(",")):
            ch, span = tok.split(":"); lo, hi = span.split("-")
            copies.append((f"{ch}{i+1}", ch, int(lo), int(hi)))
        subs = subunits_design(parse_atoms(path), copies)
        r146_res = args.r146_incopy
        lo, hi = (int(x) for x in args.walker_incopy.split("-"))
        walker_res = list(range(lo, hi + 1))
        order = [c[0] for c in copies]
        pairs = m2_adjacency(subs, r146_res, walker_res, designed_order=order)

    if args.smoke:
        pairs = pairs[:1]

    rows = []
    print(f"# {label}: {path}   (Rosetta ref2015 InterfaceAnalyzer, relax_cycles={args.relax_cycles})")
    print(f"# {'donor':>6} {'acc':>6} {'R146_A':>7} {'dG_sep(REU)':>12} {'dSASA':>8} {'n_res':>6}")
    for d, a, r146d in pairs:
        pdb = os.path.join(args.out, f"_ros_pair_{label}_{d}_{a}.pdb")
        write_pair_pdb(pdb, subs[d], subs[a])
        t0 = time.time()
        res = analyze_pair(pdb, sfxn, args.relax_cycles)
        res.update(system=label, donor=str(d), acceptor=str(a),
                   R146_min_A=round(r146d, 2), seconds=round(time.time() - t0, 1))
        rows.append(res)
        print(f"  {str(d):>6} {str(a):>6} {r146d:>7.2f} {res['dG_separated']:>12.1f} "
              f"{res['dSASA']:>8.0f} {res['n_interface_res']:>6}  ({res['seconds']}s)")
        try: os.remove(pdb)
        except OSError: pass

    dg = [r["dG_separated"] for r in rows]
    print(f"# {label} mean dG_separated = {np.mean(dg):.1f} +/- {np.std(dg):.1f} REU  (n={len(dg)})")

    cols = ["system", "donor", "acceptor", "R146_min_A", "dG_separated", "interface_dG",
            "dSASA", "dG_dsasa_ratio", "n_interface_res", "total_score", "seconds"]
    with open(os.path.join(args.out, f"{label}_rosetta_interface.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in cols})
    json.dump(dict(system=label, structure=path, method="pyrosetta_InterfaceAnalyzer_ref2015",
                   relax_cycles=args.relax_cycles, mean_dG_separated=float(np.mean(dg)),
                   std_dG_separated=float(np.std(dg)), rows=rows),
              open(os.path.join(args.out, f"{label}_rosetta_interface.json"), "w"), indent=2)
    print(f"# wrote {label}_rosetta_interface.csv")


if __name__ == "__main__":
    main()
