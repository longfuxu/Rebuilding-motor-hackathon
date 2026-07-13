# C — mechanochemical ratchet: does grip + an ATP-clock translocate DNA, and concerted vs sequential?

_2026-07-12, A100. Model driver `md/driven/tmd_ratchet.py`; comparison `md/driven/compare_ratchet.py`._

## What this is
A **minimal mechanochemical model** built to answer the questions the plain driven-MD could not:
(1) can the ring's conformational cycle be **coupled to processive DNA translocation** by adding a **grip**;
(2) how do **concerted** (all subunits fire together) vs **sequential/hand-over-hand** (fire one-by-one) differ?
We PUT IN two ingredients the earlier runs lacked: a per-subunit **power stroke on a clock** (`--mode`), and a
**catch-and-carry grip** (a spring that couples the dsDNA to a subunit only during its down-stroke, then releases
on reset). Ring+dsDNA, implicit, ~31k atoms, 4 cycles.

## Result
| mode | DNA net (4 cycles) | per-subunit z(t) signature | M2 coupling |
|---|---|---|---|
| **concerted** | **−1.6 Å** (barely moves) | in-phase / noisy, ~0 amplitude | 5→5 (min 4, flickers) |
| **sequential** | **−16.1 Å** (~1.2 bp/cycle) | **travelling wave** (each subunit peaks in turn) | 5→5 (min 4, steady) |

Figure: `ratchet_compare.png`.

## Interpretation
- **Grip is the load-bearing element.** With grip, DNA translocates (unlike the plain staircase drive, net ~0).
- **Efficient translocation requires SEQUENTIAL (hand-over-hand) firing.** Sequential always keeps ≥1 subunit
  gripping → the DNA is continuously held and cannot slip back → it ratchets steadily. Concerted releases ALL
  grips simultaneously each cycle → the DNA springs back → almost no net motion. **This is the physical reason
  processive motors use hand-over-hand: continuous grip.**
- **The 3D-MINFLUX discriminator** falls straight out: read the phase relationship of the 5 per-subunit z(t)
  traces — a **travelling wave** (staggered) ⇒ sequential + DNA substeps; **in-phase** ⇒ concerted + one big
  step or no net translocation.

## Ties to A and B
- **A** (per-position dead-seat: position-INDEPENDENT) says the special subunit is a **dynamic/timing** feature,
  not a static-structure one. **C** says efficient translocation needs a **sequential clock**. Together: the
  special subunit is a candidate **phase-setter / clock-starter** for the sequential wave — a concrete, testable
  hypothesis for the mixed-ring / MINFLUX experiments.
- **B** (continuous 8 ns stroke, DNA net +3.5 Å vs ~0 fast) says time/rate also matters — but grip + continuous
  hold (C) is what gives robust, directional stepping.

## Limitations (do not over-read)
- This is a **MODEL**: the grip and the clock are IMPOSED, not emergent. It demonstrates the **logic** (continuous
  grip ⇒ translocation; sequential ⇒ continuous grip) and the **signature to look for**, NOT that φ29 is sequential.
- The concerted mode's poor result is partly because 5 simultaneous catch-and-carry grips on one DNA fight each
  other; a different concerted grip scheme could differ. The robust, model-independent point is *continuous grip*.
- Implicit solvent; coarse (whole-subunit) grip; short (4 cycles); the −16 Å slightly exceeds the ideal 13.6 Å
  (grip over-carry) — the magnitude is illustrative, the direction/mechanism is the message.

## Next
- Replicates + a concerted grip variant that doesn't self-conflict (fairer concerted test).
- Grip on real DNA-binding residues (not whole-subunit) + explicit solvent for a physical grip.
- Render the sequential travelling-wave + DNA-translocation animation (payoff visual).
