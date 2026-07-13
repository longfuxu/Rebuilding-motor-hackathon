# Integrated helical↔planar + descent-grip — the honest field-model test

_2026-07-12, A100. `md/driven/tmd_integrated.py`, 2 modes × 2 seeds, ring+dsDNA implicit, grip on REAL DNA-contact
residues (110/226/286/127/115 heavy atoms/subunit), active only during the H→P descent, catch-and-carry._

## What this is (and why it matters)
The earlier grip ratchet was a **hand-over-hand firing** model — NOT the field's mechanism. This run puts the grip
into the actual **helical↔planar** model: the whole ring forms the helical lock-washer (P→H) then collapses to
planar (H→P), and the grip engages the DNA **only during the descent** (the power stroke), on the real
DNA-contacting residues. This is the closest in-silico test here to the φ29 spiral↔planar hypothesis, and it
addresses reviewer concerns (ii, replicates) and (iv, real grip residues).

## Result
| mode | DNA net (4 cycles, mean±sd, n=2) | per-subunit z(t) | coupling |
|---|---|---|---|
| concerted  | **−2.46 ± 0.20 Å** | symmetric staircase, all subunits in-phase | 5→4 (min 3) |
| sequential | **−3.18 ± 0.71 Å** | staggered travelling wave | 5→5 (min 4) |

Both are **small (~0.25 bp/cycle)** and **concerted ≈ sequential** (difference within/near the error bars).
Contrast: the hand-over-hand ratchet gave **−16 Å**. Figure `integrated_compare.png`.

## Interpretation — the key honest finding
**The symmetric helical↔planar cycle does NOT efficiently translocate DNA, even with a physical descent-grip.**
Why: the helical staircase is symmetric (subunits split ±2.4 Å about planar). During the H→P collapse, the grip
pulls DNA in BOTH axial directions at once (subunits above pull it down, subunits below pull it up) → near
cancellation → poor net. The hand-over-hand ratchet worked (−16 Å) precisely because it BROKE that symmetry
(all subunits stroke the same way, always ≥1 gripping).

**⇒ Efficient translocation needs GRIP + a SYMMETRY-BREAKING / DIRECTIONAL element.** The helical↔planar cycle
alone (even gripped) lacks it. This is exactly where the **special subunit** comes in: the per-position dead-seat
scan (A) showed the special subunit is *dynamic/timing*, not structural — a natural candidate for the
symmetry-breaker (a subunit that stays anchored, or that starts/phases a directional wave).

## The whole session's mechanistic arc (honest)
1. helical↔planar, NO grip (campaign/B): DNA net ~0.
2. helical↔planar + grip, SYMMETRIC (this run): still small (~2.5–3 Å); concerted ≈ sequential.
3. hand-over-hand, SYMMETRY-BROKEN (ratchet): −16 Å.
→ **Translocation = grip × broken symmetry.** The symmetric collective conformational cycle is not enough.

## Sharpened prediction for the lab
MINFLUX should test not just concerted-vs-sequential timing but whether the per-subunit motion is **symmetric**
(a poor translocator) or carries a **directional/broken-symmetry signature** (efficient) — and whether the
special subunit is the one that breaks it. The addressable single chain lets you put the special/anchor subunit
at a chosen seat and read whether translocation efficiency tracks its position.

## Caveats (do not over-read — still a model)
Grip + clock are imposed; the −2.5/−3.2 Å magnitudes are parameter-dependent (kgrip, cycle_ps) and small/noisy
(n=2); implicit solvent; no ATP·Mg; the grip is a spring, not emergent electrostatics. The **symmetry argument is
robust**, the exact magnitudes are not. Next: an emergent grip (explicit solvent + Mg + real contacts) and an
anchored-special-subunit variant to test the symmetry-breaker hypothesis directly.
