#!/usr/bin/env python3
"""Inspect downloaded ring-motor structures: per (model,chain) list protein vs
nucleic content, residue range, CA count. Prints a compact report so we can pick
the correct file + ring chains + substrate chains for each panel protein.

Handles both legacy PDB biological-assembly files (.pdb1, symmetry copies in
separate MODEL records) and deposited mmCIF (.cif). Uses only CA atoms (protein)
and P atoms (nucleic) so it is robust to CA-only entries (e.g. 3EZK)."""
import sys, os, warnings, glob
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser
warnings.simplefilter("ignore")

AA3 = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS",
       "MET","PHE","PRO","SER","THR","TRP","TYR","VAL","MSE","SEC","PYL"}
NUC = {"DA","DT","DG","DC","DU","A","U","G","C","I","N"}

def load(path):
    if path.endswith(".cif") or "assembly1.cif" in path:
        p = MMCIFParser(QUIET=True)
    else:
        p = PDBParser(QUIET=True, PERMISSIVE=True)
    return p.get_structure(os.path.basename(path), path)

def inspect(path):
    try:
        s = load(path)
    except Exception as e:
        print(f"  !! parse fail {os.path.basename(path)}: {e}")
        return
    n_models = len(list(s))
    rows = []
    for mi, model in enumerate(s):
        for ch in model:
            prot = []   # (resnum, ca_xyz)
            nuc = 0
            for res in ch:
                rn = res.id[1]
                nm = res.resname.strip()
                if nm in AA3 and "CA" in res:
                    prot.append(rn)
                elif nm in NUC:
                    nuc += 1
            if prot:
                rows.append((mi, ch.id, "prot", len(prot), min(prot), max(prot)))
            elif nuc >= 3:
                rows.append((mi, ch.id, "nuc", nuc, 0, 0))
    # summarize
    prot_rows = [r for r in rows if r[2] == "prot"]
    nuc_rows  = [r for r in rows if r[2] == "nuc"]
    print(f"\n=== {os.path.basename(path)}  models={n_models}  protein_subunits={len(prot_rows)}  nucleic_chains={len(nuc_rows)}")
    # cluster protein subunit lengths
    lens = sorted(set(r[3] for r in prot_rows))
    for r in prot_rows[:20]:
        print(f"    m{r[0]} ch{r[1]:<3} prot n={r[3]:<5} res {r[4]}-{r[5]}")
    if len(prot_rows) > 20:
        print(f"    ... (+{len(prot_rows)-20} more protein subunits)")
    if nuc_rows:
        print(f"    nucleic: " + ", ".join(f"m{r[0]}ch{r[1]}(n={r[3]})" for r in nuc_rows[:12]))

if __name__ == "__main__":
    d = "/Users/longfu/Developer/claude-science-hackthon/.claude/worktrees/export-gp16-results/gp16_design/outputs/buildability_atlas/structures"
    args = sys.argv[1:]
    if args:
        files = [os.path.join(d, a) for a in args]
    else:
        files = sorted(glob.glob(os.path.join(d, "*.pdb1"))) + sorted(glob.glob(os.path.join(d, "*.cif")))
    for f in files:
        inspect(f)
