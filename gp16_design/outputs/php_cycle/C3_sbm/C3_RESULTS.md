# C3 SBM (coarse-grained Gō) — result: INCONCLUSIVE (2026-07-11)

Consolidated in `../../PHP_CYCLE_PROGRESS_REPORT.md`. Model: hand-rolled dual-basin Cα-Gō in OpenMM
(`pipelines/php_cycle/c3_sbm.py`), 1635 matched beads (P=A_apo, H=7JQQ; P→H Cα-RMSD 43 Å), run on GCP A100.

## Two attempts, both inconclusive
- **v1** — union dual-basin (P+H contacts), CV = fraction of ALL H-contacts. Q stuck at **0.68** because P
  and H share ~60 % of contacts (intra-domain) → the CV is insensitive to the inter-subunit rearrangement,
  and the union potential is frustrated (competing P/H contacts).
- **v2** — single-basin-on-P + CV = fraction of **H-UNIQUE** contacts (2063 of them) + steered by a global
  RMSD-to-H moving restraint. Q_Huniq stayed at **~0.25** as r0 was ramped 4.3→2.25 nm: the global RMSD CV
  drags the ring toward H **without local rearrangement** (H-unique contacts don't form). **Spot-preempted at
  ~win 30 and the VM auto-deleted → run lost (no checkpoint).**

## Why (converges with C2)
C2 found P↔H is a **collective multi-mode** motion (best single mode ~33 %). A **single global CV (RMSD or
contact-fraction) + simple steer cannot induce it** — confirmed here. Track B (all-atom) reached the design
only ~⅓ of the way with a monotonic staircase CV, for the same reason.

## Fixes for the next run
1. **CV**: use the per-subunit **axial-staircase target** CV (the one that actually drove Track B), not a
   global RMSD/contact CV; or a multi-mode CV built from the C2 mode set.
2. **Sampling**: umbrella / metadynamics, not a single monotonic steer.
3. **Preemption-safe**: checkpoint + `scp` the incremental `steer_refined.json` during the run (Spot deletes
   the VM on preemption).
4. **DNA**: protein-only is incomplete (P→H is DNA/ATP-gated) — add DNA (AWSEM+3SPN.2C, or a DNA-grip SBM).

Inputs kept: `P_matched.pdb`, `H_matched.pdb` (1635 matched Cα, same bead order).
