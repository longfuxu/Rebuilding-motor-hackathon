# Elastic-network normal-mode analysis (NMA) of the gp16 ring + overlap with the real transition

**Tool:** Anisotropic Network Model (ANM) — springs between CA within 13 Å; diagonalize; the softest modes
are the ring's cheapest collective motions. Run on the apo closed-planar native ring (Boltz), 1635 CA (res 4–330).
Code: `reproduce/anm.py`; overlap `reproduce/overlap.py` (numpy/scipy, seconds, free).

## Raw modes (character of the softest motions)
The two softest non-trivial modes (7, 8) are **out-of-plane / helical** (helical fraction 0.88, 0.90).
**But this alone is weak evidence** — out-of-plane bending is the softest motion of *any* flat ring (a drumhead),
so "softest mode is out-of-plane" is partly generic, not gp16-specific.

## The rigorous test — overlap of the soft modes with the ACTUAL apo→7JQQ(helical) transition
Aligned the apo ring to the experimental helical state 7JQQ (best cyclic chain assignment; ring RMSD **6.6 Å** —
they genuinely differ = the planar↔helical change). Projected that difference vector onto the ANM modes:

| softest modes used | cumulative overlap | fraction of the transition captured |
|---|---|---|
| 2 | 0.11 | **1%** |
| 5 | 0.44 | 19% |
| 10 | 0.57 | 33% |
| 20 | 0.64 | **41%** |

Best single mode = mode 10 (overlap 0.33); no single mode dominates. (Random overlap for 20 modes in this space
≈ 0.06 ≈ 0.4%, so 0.64 is a strong enrichment.)

## Honest conclusion (this corrects the raw-NMA first impression)
- The **very softest modes are NOT the transition** (softest 2 capture only 1%) — they are the generic flat-ring bending.
- The **planar→helical opening IS enriched in the low-energy modes** (~41% captured by the 20 softest, vs ~0.4%
  random) but is **moderately collective, not a single clean soft mode**.
- So: the transition lies **partly along the ring's low-energy collective directions** (the ring is *somewhat*
  predisposed to open), but it is not "one soft mode = the opening."
- **Caveat:** the apo structure is a *prediction*, so part of the unexplained 59% is apo-prediction-vs-cryoEM
  noise (not real high-frequency motion) — the true functional-transition overlap is likely somewhat higher.

## Why this tool, and what's next
Static predictors (AF3/Boltz/OF3) give a closed ring for ATP and cannot show the opening; plain all-atom MD
cannot reach the transition timescale; NMA reads the intrinsic soft directions and (with the overlap) *quantifies*
how soft the real transition is — free, minutes. For the actual dynamics/path/barrier of the transition, the
next rung is a **Cα structure-based (Gō) dual-basin model** (apo + 7JQQ as the two minima) on GPU (~1 h) — see
the MD escalation ladder in HANDOFF.
