#!/usr/bin/env python3
"""Structural credibility of the cp233 single-chain design: superimpose it on the native
gp16 ring and quantify fold + ring fidelity. Answers "is cp233 actually a good design or
is it off from the apo structure?"

Inputs (from the MD validation set):
  A_apo.pdb      native apo gp16 ring (5 chains x 332)   -- reference
  C_design.cif   cp233_int15_inter10 single-chain design (1 chain, 1750 CA = 5 CP copies)
  7JQQ.cif       EXPERIMENTAL ATP-helical gp16 (predictor-independent anchor)

Method: US-align with -cp (circular-permutation-aware) for the subunit fold; then apply the
single-subunit transform to the whole design and check that all 5 copies fall onto the 5
native subunits (tests the C5 ring geometry, not just the fold).

Requires: USalign (compiled from zhanggroup.org/US-align), biotite, matplotlib, numpy.
Results (2026-07-09):
  subunit fold vs native apo : RMSD 1.80 A, TM 0.94, 330/332 aligned, seq_id 0.99
  subunit fold vs EXP 7JQQ   : RMSD 3.72 A, TM 0.79 (larger b/c 7JQQ is the ATP-helical conformer)
  ring geometry              : outer radius 44 vs 42 A; all 5 design copies -> 5 distinct native subunits
  => cp233 recapitulates the native fold AND the pentameric ring. The only non-native density is
     the CP linker + relocated termini (by design, to make it single-chain).
"""
# (the run used biotite for parsing + a US-align -m matrix applied to the full design;
#  see usalign_cp233_vs_nativeA.txt for the raw subunit alignment. Kept as a record.)
