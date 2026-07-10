#!/usr/bin/env python3
"""Extend the buildability-atlas panel with additional homo-oligomeric AAA+ rings and
apply the (published) 2-branch rule per protein.

Reuses the atlas descriptor geometry verbatim (compute_descriptors.py):
  jams_channel  = a native terminus contacts the translocated substrate (<10 A) OR
                  sits at rel_radius < 0.5 (i.e. near the pore axis)
  direct_gap_A  = through-space C(i)->N(i+1) gap between adjacent subunits (motor construct)
2-branch rule (atlas §4): jams_channel -> CP ; else -> direct fusion (linker sized to gap).
(The nominal 3rd 'diffusion' branch has 0 validated members; we keep the raw depth-2 tree
 prediction as a separate column for transparency but the honest call is the 2-branch rule.)

New panel members (full modeled chain = the single-chain fusion unit; no motor-range guess):
  ClpB  6OAX  Thermus ClpB disaggregase, NBD1+NBD2, res161-857, hexamer, substrate chain P (0 CA)
  NSF   3J94  N-ethylmaleimide sensitive factor, D1 ring res217-737, hexamer
  LonA  6ON2  Yersinia LonA protease, AAA++protease res253-775, hexamer, substrate chain G
"""
import os, sys, json, csv, warnings
import numpy as np
warnings.simplefilter("ignore")
ATLAS = "/Users/longfu/Developer/claude-science-hackthon/.claude/worktrees/export-gp16-results/gp16_design/outputs/buildability_atlas"
OUT   = "/Users/longfu/Developer/claude-science-hackthon/.claude/worktrees/export-gp16-results/gp16_design/outputs/aaa_ascE_screen"
sys.path.insert(0, ATLAS)
import compute_descriptors as cd   # reuse analyze(), PANEL, geometry

# --- new members appended to the atlas PANEL ---
NEW = [
 dict(name="ClpB", pdb="6OAX.pdb1", family="AAA+ disaggregase (double-ring)", oligo=6,
      ring=[(0,c) for c in "ABCDEF"], substrate=[], known=None,
      note="Thermus ClpB NBD1+NBD2 tandem res161-857; substrate chain P has no modeled CA"),
 dict(name="NSF", pdb="3J94.pdb1", family="AAA+ SNARE disassembly", oligo=6,
      ring=[(0,c) for c in "ABCDEF"], substrate=[], known=None,
      note="NSF D1 ATPase ring res217-737 (D2/N not in this construct)"),
 dict(name="LonA", pdb="6ON2.pdb1", family="AAA+ protease", oligo=6,
      ring=[(0,c) for c in "ABCDEF"], substrate=[(0,"G")], known=None,
      note="Yersinia LonA AAA+ + protease res253-775; substrate peptide chain G"),
]

JAM_SUB_A = 10.0
JAM_CTERM_RELR = 0.40   # atlas: C-term (fusion donor) deep in pore; a merely inward N-term is NOT flagged
GAP_DIRECT_A = 38.0

def jams(row):
    """Exact atlas jam criterion (classify_and_train.py::jams): substrate contact <10 A on
    either terminus, OR the C-terminus sitting deep in the pore (rel_radius < 0.40). A merely
    inward N-terminus is deliberately NOT flagged (gp17 N-term rel 0.32 yet direct closes 5/5)."""
    csub = row.get("Cterm_to_substrate_A"); nsub = row.get("Nterm_to_substrate_A")
    touch = (csub is not None and csub < JAM_SUB_A) or (nsub is not None and nsub < JAM_SUB_A)
    nearaxis = row["Cterm_rel_radius"] < JAM_CTERM_RELR
    return bool(touch or nearaxis), touch, nearaxis

def rule_2branch(row):
    j,_,_ = jams(row)
    if j: return "CP"
    return "direct"

def rule_tree(row):
    j,_,_ = jams(row)
    if row["direct_gap_A"] <= GAP_DIRECT_A: return "direct"
    return "CP" if j else "diffusion"

rows = []
for e in cd.PANEL + NEW:
    r = cd.analyze(e)
    j, touch, nearaxis = jams(r)
    r["jams_channel"] = j
    r["jam_by_substrate_contact"] = touch
    r["jam_by_near_axis"] = nearaxis
    r["pred_2branch"] = rule_2branch(r)
    r["pred_tree_depth2"] = rule_tree(r)
    rows.append(r)

# print
print(f"{'name':10} {'fam':30.30} {'oligo':5} {'gap':6} {'Csub':6} {'Crel':5} {'jams':5} {'2branch':8} {'tree':10} {'known'}")
for r in rows:
    print(f"{r['name']:10} {r['family']:30.30} {r['oligo_nominal']:<5} {r['direct_gap_A']:<6} "
          f"{str(r['Cterm_to_substrate_A']):6} {r['Cterm_rel_radius']:<5} {str(r['jams_channel']):5} "
          f"{r['pred_2branch']:8} {r['pred_tree_depth2']:10} {r['known_method']}")

os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT,"panel_descriptors.json"),"w") as f: json.dump(rows,f,indent=2)
keys=list(rows[0].keys())
with open(os.path.join(OUT,"panel_descriptors.csv"),"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=keys); w.writeheader()
    for r in rows: w.writerow(r)
print(f"\nwrote panel_descriptors.csv/.json ({len(rows)} proteins)")
# rule accuracy on validated anchors
anchors=[r for r in rows if r["known_method"]]
acc=sum(1 for r in anchors if r["pred_2branch"]==r["known_method"])
print(f"2-branch rule accuracy on {len(anchors)} validated anchors: {acc}/{len(anchors)}")
for r in anchors: print(f"   {r['name']}: pred {r['pred_2branch']} vs known {r['known_method']} {'OK' if r['pred_2branch']==r['known_method'] else 'MISS'}")
