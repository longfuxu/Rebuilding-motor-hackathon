# Elastic-network normal-mode analysis (NMA) of the gp16 ring

**Tool:** Anisotropic Network Model (ANM) — springs between CA atoms within 13 Å; diagonalize; the
lowest-frequency (softest) modes are the ring's cheapest collective motions. Run on the apo closed-planar
native ring (Boltz), 1660 CA. Code: `reproduce/anm.py` (numpy/scipy, ~seconds, free).

**Result:** the two softest non-trivial modes (7, 8) — the ring's very lowest-energy motions — are
**out-of-plane / helical** (helical fraction 0.88, 0.90); modes 10, 11 are strongly helical too (≈1.00).
Radial "breathing" appears only at higher modes (12, 15).

**Interpretation:** the **planar → helical (lock-washer) opening lies along a soft (low-energy) direction**
of the ring. This is consistent with the motor cycling planar↔helical, and with the idea that ATP does not
force a high-energy new conformation — it **biases the ring along a pre-existing soft mode**.

**Honest caveat:** out-of-plane bending is *generically* the softest motion of any flat ring-shaped object
(like a drumhead), so this is **suggestive, not proof** of a gp16-specific designed opening. The rigorous
next step is to compute the **overlap** between these low modes and the actual apo→7JQQ(helical) difference
vector — high overlap would show the observed transition IS a soft mode. (Needs 7JQQ coordinates aligned.)

**Why this is the right tool:** the static structure predictors (AF3/Boltz/OF3) gave a closed ring for ATP
and can't show the opening; plain all-atom MD can't reach the transition timescale; NMA directly reads the
intrinsic soft collective motions. For the actual dynamics of the transition, coarse-grained MD is next.
