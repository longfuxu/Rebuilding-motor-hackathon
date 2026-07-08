# Topology Failure Taxonomy — φ29 gp16 Single-Chain Design
_Generated 2026-07-07 · from 7JQQ geometry (helical state) + idealized planar C5 model_

This is the map of *why* designs fail — it prevents blind synthesis and defines the negative controls.

## F1 — Native-terminus fusion fouls the DNA channel (topology failure)
The native C-terminus (res 330) sits **6.2 Å from the dsDNA** — it is a DNA-channel contact. Direct
head-to-tail fusion routes a linker across the translocation channel. **All `nativefus_*` constructs are
expected to fail and serve as negative controls** (composite 0.24–0.36).

## F2 — Native-terminus span too long (topology failure)
Native C_i→N_{i+1} spans **54–66 Å** (max at the A→E seam). No practical short linker bridges this without
either strain or a long floppy insert that mislocalizes. Circular permutation at res 280/297 cuts the required
junction to **24–27 Å**, inside GS-linker reach.

## F3 — Linker splints the planar↔helical transition (functional failure)
The decisive state-dependent failure. 7JQQ is the **open helical** state; the resting ring is hypothesized
**planar**. A linker sized to the planar junction is **too short to permit opening**; one sized rigidly to the
helical state resists closure. The per-junction **stroke** is the tolerance a linker must absorb:
native termini 13 Å, CP@297 6.2 Å, **CP@280 3.0 Å (best)**. Design rule: size for the helical maximum + flexible
excess so the linker is compliant, never load-bearing.

## F4 — CP cut in a functional element (topology + functional failure)
Cutting in Walker A (24–31), Walker B (115–119), the trans arginine finger R146/loop 104–109, ATP pocket, DNA
channel, or pRNA interface destroys catalysis or substrate binding. cut@60 is rejected on this basis (5.5 Å to DNA).

## F5 — Global fold recovered but trans-contacts lost (functional failure)
A predicted structure can show a plausible ring RMSD while the *inter-subunit* catalytic registry (the trans R146
into the neighbor's pocket) is broken. **Do not rank on global fold alone** — score the arginine-finger geometry,
ATP-pocket integrity, and pRNA/DNA contact preservation separately (Aim-3 folding stage).

## F6 — Special-interface not addressable without changing state (scientific outcome, not a bug)
If the A→E "special" interface can only be made single-chain-addressable by locking the ring into one
conformational state, that is itself the finding — the special-subunit role may be dynamic rather than fixed.

## Where each failure is caught
| Stage | Catches |
|---|---|
| Static topology (this analysis, CPU) | F1, F2, F4, and the F3 stroke estimate |
| Structure prediction (Boltz-2/OpenFold3, GPU) | F5, refined F3, fold feasibility of long chains |
| Wet-lab (expression → assembly → motor assay) | true F3/F5, expression-length limits |
