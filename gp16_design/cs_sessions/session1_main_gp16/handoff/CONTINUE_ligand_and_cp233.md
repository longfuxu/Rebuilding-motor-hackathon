# Continue — AF3 ligand states + the cp233 CP winner

## 1. AF3 ligand-state control (apo/ADP/ATP, template OFF) — done, closed-book
Structures in `handoff/ligand_states/`. Result (scored by Claude Code): apo/ADP/ATP are ALL closed planar
(planarity 0.01–0.08 Å; R146 ~7–8.6 symmetric) — **even ATP template-off does NOT open the ring**. So AF3
fundamentally can't model the ATP-driven opening; 7JQQ (3-ATP) stays the experimental open-state ground truth.
Don't over-invest here — just fold this into the ligand-state section of RUN_SUMMARY/folding_results as a clean
3-state predictor-limit control. Detail: `handoff/af3_ligand_states.md`.

## 2. cp233 is the CP winner — next fold task = E119Q addressability on it
Your CP task nailed it: **cp233_int15_inter10 closes the ring under both Boltz (3/3 seeds 5/5) and OpenFold3
(5/5)** — the storyline's CP rung succeeded. The user is running AF3 on cp233_int15_inter10 (the decisive 3rd
predictor); when those come back Claude Code will score them.

**Your next fold task (highest value):** run the **E119Q addressability experiment on the cp233 winner** — the
payoff on the best single-chain construct. Put the catalytically-dead Walker-B E119Q at 1–2 DEFINED copies of the
cp233 single chain (note: CP shuffles the sequence, so E119=res119 sits inside part2 = the res4..233 half of each
copy — compute its construct position from the manifest), fold Boltz-2 (≥3 seeds) + OpenFold3, score sequential
M2 + M1 in ring context. Question: does the ring still close with the dead seat in place (assembly tolerance),
and is any perturbation local to the known seat? Compare to WT cp233 and to B1_L40_E119Q. Free NIM, sequential
M2 not global pTM, save incrementally.
