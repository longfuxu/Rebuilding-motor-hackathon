# C2 agent prompt — is planar↔helical an intrinsic SOFT MODE of the gp16 ring? (ENM/ANM)

_Paste to a fresh agent. Self-contained. Part of Goal 4 (computationally test the φ29 P→H→P cycle);
full plan: `gp16_design/COMPUTATIONAL_AIMS_PHP_CYCLE.md`. This is the cheap, CPU-only, runs-on-a-laptop
tier — no GPU. It does NOT need C1/C3._

## The question
The P→H→P model says the ring cycles between a **planar** ring (P) and a **helical staircase** (H).
**Is the P↔H transition a low-frequency (soft) collective mode of the ring architecture?** — a necessary
(not sufficient) condition for an easy, ATP-tippable transition. And: **does the single-chain DESIGN retain
this mode**, or did the covalent linkers stiffen it out? (This connects to the "don't make a rigid brick"
concern and to whether C3's driven transition is well-motivated.)

Context worth knowing: both AF3 and Boltz predict the ring **planar for apo/ATP/ADP** (they don't open with
ATP — a predictor-basin bias, see `outputs/php_cycle/`). So the accessibility of the H state must be probed
by physics, not prediction. ENM is the cheapest such probe.

## Inputs (already in the repo)
- **P (planar) reference:** `gp16_design/md/openmm_validation/inputs/A_apo.pdb` (native apo ring, ~1660 Cα).
- **H (helical) reference:** `gp16_design/md/openmm_validation/inputs/B_7jqq_helical.pdb` (7JQQ, ~1635 Cα).
  Calibration already measured: 7JQQ subunit-centroid axial spread ≈ **4.8 Å** (planar ≈ 0.1 Å).
- **Design:** `gp16_design/md/openmm_validation/trajectories/C/C_start.pdb` (cp233 single chain, 1750 Cα).
- Existing tooling to reuse or check: `gp16_design/reproduce/anm.py`, `overlap.py`, `enm_openclose.py`.

## Method
1. **Build the P→H difference vector (native):** map A_apo ↔ B_7jqq by residue (both are native gp16
   sequence, different conformations), superpose on Cα (Kabsch), Δ = H − P (per-Cα displacement). Normalize.
   Report the per-subunit axial component of Δ (should be the ~planar→staircase motion; sanity vs the 4.8 Å
   centroid spread).
2. **ANM on the planar state** (A_apo, Cα-only, standard cutoff ~13–15 Å; use **ProDy** `calcANM` if simplest,
   else the repo's `anm.py`). Get the lowest ~20 non-trivial modes.
3. **Overlap:** compute |mode_i · Δ̂| for each low mode + the **cumulative overlap** of the first k modes
   (ProDy `calcOverlap`, `calcCumulOverlap`). Report: which mode(s) carry P→H, and the cumulative overlap of
   the lowest 1/3/5/10 modes. High cumulative overlap in the lowest few modes ⇒ P↔H is an intrinsic soft mode.
4. **Design vs native:** ANM on the design `C_start.pdb`; assess whether it has an analogous low-frequency
   opening / planar↔helical mode (e.g., project the design's low modes onto a design-equivalent axial-spread
   / ring-opening vector, or compare the softness/character of its lowest modes to native's). The design is a
   different topology (covalent single chain) so an exact Δ mapping isn't possible — instead compare **mode
   character** (is there a low-frequency ring open/close mode?) and the **effective stiffness** of that mode.
5. (Optional) `enm_openclose.py` may already implement an open↔close ENM analysis — check and reuse it.

## Deliverables (write to `gp16_design/outputs/php_cycle/C2_enm/`)
- `C2_REPORT.md`: the overlap table (mode index, frequency/eigenvalue, overlap, cumulative), **which mode =
  P→H and its rank/softness**, native-vs-design comparison, and an **honest** read: ENM tells you the motion
  is *geometrically accessible / low-energy*, NOT its driven direction or absolute energetics; a high overlap
  supports that P↔H is built into the fold; a low overlap would be the surprising (and important) result.
- A figure: overlap-per-mode bar + cumulative-overlap curve (native), and the design's low-mode spectrum.
- The Δ vector and mode arrays saved (npz) for reuse by C3/C4.

## Tooling / compute
- Prefer **ProDy** (`pip install prody`): `parsePDB`, `calcANM`, `calcOverlap`, `calcCumulOverlap`,
  `superpose`. Fall back to the repo's `reproduce/anm.py` + `overlap.py`.
- **CPU-only, ~minutes, runs on the MacBook** (no GPU, no GCP). ~1660–1750 Cα ANM is small.

## Honest scope
ENM/soft-mode overlap is a **necessary-condition** test (is the motion accessible + how soft), not proof of
the mechanism. It feeds C3 (the enhanced-sampling / coarse-grained free-energy landscape) and C4 (concerted
vs sequential). Report cleanly whichever way it comes out.
