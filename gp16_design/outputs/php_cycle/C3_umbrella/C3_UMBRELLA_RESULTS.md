# Umbrella-sampling — gp16 single-chain ring, P→H landscape (FORWARD pass)

**Run:** A100 Spot, OpenCL, implicit GBSA-OBC2, 27,770-atom single-chain ring (cp233, chain A only — NOT the
native pentamer). Staircase axial CV, `--passes both` launched; **forward 24/24 windows completed**, reverse
**stopped at the quality gate** (see below). `k_eff=4248 kJ/mol/nm²`, kT=0.596 kcal/mol. Started 07:33 UTC,
forward done ~15:57 UTC. Portable `state_fwd_end.xml` saved → reverse resumable anytime.

## ⚠️ Gate verdict: the absolute free energy is NOT reliable (windows too stiff → poor overlap)
- WHAM gave ΔG(H−P) = **+81 kcal/mol**, barrier ≈ 62 kcal/mol — but the targeted-MD **fast-pull work was only
  35–53 kcal/mol**, and a reversible ΔG cannot exceed the dissipative work bound. The numbers are inconsistent.
- **Root cause:** `kcv=30000` is too stiff → each window's ξ distribution is σ≈0.023 while adjacent window
  centres are spaced 0.050 apart → **overlap ratio 1.07 (>1 = insufficient overlap)**. WHAM is under-determined
  between windows, so the absolute ΔG/barrier is a sampling artifact, not physics.
- **Why the reverse was stopped:** it uses the same kcv/spacing → the same overlap gap → it cannot produce a
  trustworthy hysteresis either. Continuing it would spend ~8 h/$10 to re-confirm this caveat. Decision: stop.
- **Fix for a trustworthy G(ξ):** re-run with **softer kcv (~3–4× lower, e.g. 8000)** and/or ~2× more windows so
  adjacent windows overlap, forward **and** reverse. (~8 h/$10, one A100 run.) `tmd_umbrella.py --passes both`.

## ✅ What IS reliable (per-window ensemble averages — independent of WHAM overlap)
### (a) qualitative shape
Planar is a deep basin; helical is strongly uphill — the **apo design ring does not spontaneously go helical**
(needs the ATP/occupancy drive). Directionally consistent with C1 (predictors all planar) and the targeted-MD.
### (b) M2 coupling SURVIVES the whole P→H drive  ★ core result
n_engaged 5.0 → min **4.79** → 4.99 across all 24 windows. The R146→neighbour-WalkerA coupling stays engaged as
the ring is driven planar→helical; the single chain does not shed its coupling to move toward the helical basin.
### (c) pathway = graded, UNEQUAL staircase  ★ ties to the MINFLUX prediction
per-subunit axial travel: **−2.64 / −1.53 / −0.05 / +1.60 / +2.62 Å** (subunit 3 ≈ hinge). Monotonic and
**unequal-tread** — the graded, position-dependent per-subunit motion the dwell-burst model predicts (not five
identical steps). half-travel-ξ spread = 0.20 (some sequential character, not fully concerted).

## Figures / data
- `pmf.png`, `pmf_forward.csv` — the forward WHAM PMF **(shape only; absolute value unreliable — see caveat)**
- `coupling_vs_xi.png` — M2 engagement vs ξ (reliable)
- `per_subunit_z.png` — per-subunit axial z vs ξ (reliable; the graded staircase)
- `umb/window_data.json`, `umb/state_fwd_end.xml` — raw windows + reverse resume seed

_Caveat stack: idealized monotonic 1-D staircase CV (C2: the real transition is multi-mode); implicit solvent;
7JQQ is a 3-ATP partial intermediate (~4.8 Å), not the full ~35 Å stroke; absolute ΔG under-sampled (overlap).
This is a mechanism/pathway result, NOT an absolute rate or a converged barrier._
