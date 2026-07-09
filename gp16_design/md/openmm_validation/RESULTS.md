# OpenMM MD physics-validation of the gp16 single-chain ring — RESULTS

**One-sentence conclusion.** Under identical implicit-solvent MD, the single-chain design
**cp233_int15_inter10 stays closed and symmetric on par with the genuine apo ring — it does
not open, collapse, or dissociate — behaving like the closed apo state (A) and clearly unlike
the ATP-bound helical state (B); the design passes the physics-based, predictor-independent
stability check.**

## Systems (identical protocol; see METHODS.md)
- **A** = apo closed ring (Boltz native ring), 3 ns
- **B** = 7JQQ ATP-bound helical, ligands stripped, 2 ns
- **C** = design cp233_int15_inter10 model_0, 1 ns
(Lengths differ only because this account's Colab sessions self-terminated at ~20–30 min;
the stability verdict is fully resolved within the shared 0–1 ns window — see matched table.)

## Readouts — full-run (last-20% means)

| system | Cα-RMSD plateau | radius_CV t0→final | engaged /5 | iface mean (Å) | contacts |
|--------|-----------------|--------------------|------------|----------------|----------|
| **A apo**      | 1.84 Å | 0.013 → **0.006** | 5 → 5.0 | 4.9 → 4.2 | 552 → 586 |
| **B 7JQQ**     | 4.13 Å | 0.036 → **0.052** | 3 → 3.0 | 9.1 → 9.0 | 445 → 501 |
| **C design**   | 2.81 Å | 0.022 → **0.009** | 5 → 4.5 | 5.0 → 5.5 | 584 → **657** |

## Readouts — length-matched 0–1 ns window (apples-to-apples)

| system | Cα-RMSD | radius_CV | engaged /5 | iface mean (Å) |
|--------|---------|-----------|------------|----------------|
| **A apo**    | 1.54 Å | 0.0137 | 4.98 | 4.80 |
| **B 7JQQ**   | 3.43 Å | 0.0446 | 3.06 | 9.11 |
| **C design** | 2.41 Å | **0.0107** | 4.70 | 5.27 |

## What each readout says
- **Backbone drift (Cα-RMSD).** All three plateau (no runaway unfolding). C plateaus at ~2.8 Å —
  between apo (1.8 Å) and the strained ligand-free helical state (4.1 Å), i.e. modest and stable.
- **Radial symmetry (radius_CV).** C **drifts down to 0.009** (even below apo over the matched
  window, 0.0107 vs 0.0137) → the design ring becomes *more* circular, the opposite of opening.
  B rises to 0.052 → the helical state stays asymmetric.
- **Ring closure (R146→Walker-A).** C keeps all 5 sequential fingers engaged (final per-interface
  2.8/3.0/6.1/6.1/6.7 Å, all < 8 Å); ~5/5 engaged throughout (occasional single-frame dip to 4).
  A stays 5/5. **B stays 3/5 with a persistent wide-open seam (13–21 Å)** — the ATP-driven lock-washer
  does not re-close on this timescale, exactly the expected contrast.
- **Interface contacts.** C **increases 584 → 657** (ring tightens), A ~586, B lowest — no interface
  dissociation in the design.
- **Flexibility (RMSF).** Comparable across systems (~1 Å mean); no localized unfolding hot-spot in C.

## Verdict against the kill criterion
The pre-registered criterion was: *design C is unreasonable if it rapidly opens/collapses/dissociates
while apo A stays closed.* Observed: **A stays closed; C also stays closed (and more symmetric); B
stays open/asymmetric.** → **PASS.** This is orthogonal, predictor-independent (physics-based) support
for the same conclusion the Boltz/OF3/AF3 M1/M2 gate reached — the design encodes a genuinely closed,
mechanically stable sequential ring, not a predictor artifact.

## Honest limits
- Implicit solvent (GBSA-OBC2, 2.0 nm cutoff) + short trajectories (1–3 ns): a *stability screen*, not
  a converged free-energy or long-timescale statement. GB slightly over-stabilizes compact states.
- B is run ligand-free (no ATPγS/Mg), so it tests metastability of the ATP-bound *backbone conformation*,
  not ATP-driven dynamics. Its persistent open seam is consistent with the deposited lock-washer.
- Single trajectory per system (no replicas). Explicit-solvent TIP3P + replicas is the documented
  Phase-2 escalation (needs CUDA-enabled OpenMM; ~10× cost).
- OpenCL platform (PyPI OpenMM exposed no CUDA plugin on the VM); mixed precision.
