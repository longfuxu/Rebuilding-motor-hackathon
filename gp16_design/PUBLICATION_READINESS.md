# Publication readiness — gp16 single-chain in-silico P→H→P mechanism

_Honest self-assessment for longfuxu. Companion to `SESSION_LOG.md` (the full computational trail). Bottom line:
there is a **real, publishable paper here** — but it is a **design + controlled-computational-platform + falsifiable
predictions** paper, NOT a "we proved φ29's mechanism by MD" paper. Frame it as the former and it stands; frame it
as the latter and it will (correctly) be rejected._

## 1. What is solid vs what is not
| tier | claim | status |
|---|---|---|
| 🟢 **publishable now** | the **addressable single-chain** design closes the correct ring across 3 predictors (Boltz-2/OF3/AF3) and encodes the designed sequential R146→WalkerA coupling | solid (prior work, predictor-independent) |
| 🟢 | the single chain lets you do things the native oligomer can't: **defined asymmetric occupancy, a positionable special/dead subunit, monomer-only ensemble** | solid (design fact) |
| 🟢 | **feasibility/plausibility**: the ring can be driven P→H with the M2 coupling retained; the per-subunit lift is a **graded, unequal, reproducible** staircase | solid AS a driven-MD plausibility + geometry result (honestly caveated) |
| 🟡 | conformational cycling alone does NOT translocate DNA; **time/rate matters partly** | suggestive (n small, driven, implicit) |
| 🟡 | a single dead seat is **position-independent** by static coordination → special subunit is dynamic/timing | suggestive (static predictor + geometric proxy) |
| 🔴 **NOT a standalone result** | the **grip ratchet** numbers (sequential −16 Å vs concerted −1.6 Å) | illustrative model — shows a known PRINCIPLE + a MINFLUX signature, not a discovery about φ29 |

## 2. The four reviewer objections — status + what we did + what's still needed
**(i) Driven-CV artificiality (steered work ≠ ΔG; a 1-D CV imposes concerted / a chosen path).**
- Reality: steered/targeted MD gives a *dissipative upper bound*, not ΔG; the staircase CV forces a monotonic path.
- Done here: we state W is an upper bound everywhere; C2 (ENM) showed the transition is multi-mode so a 1-D CV
  mis-estimates barriers; we do NOT report a ΔG (the umbrella's absolute number was flagged unreliable).
- Still needed: converged **enhanced sampling** (well-tempered **metadynamics** / string method on a difference-RMSD
  or 2-D CV) for any thermodynamic claim; OpenMM ships metadynamics (no PLUMED). Until then: no rate/ΔG claims.

**(ii) Insufficient sampling / no replicates (umbrella window overlap 1.07; ratchet n=1/mode).**
- Done here: the cycle campaign used **3+3 seeds** (treads ±0.08 Å reproducible; the DNA-null survived replication
  and killed a single-seed false positive). The **integrated run uses 2 modes × 2 seeds** (`tmd_integrated.py`).
- Still needed: soften umbrella `kcv≈8000` (or 2× windows) for adequate overlap + WHAM error bars; ≥3–5 seeds and
  convergence checks on every mechanism figure; report all replicates, never a single trajectory.

**(iii) Grip is imposed, not emergent.**
- Reality: the grip (and the ATP clock) are hand-coded restraints — they encode a hypothesis, they don't discover it.
- Done here: we label every grip result 🔴 illustrative; the integrated run at least anchors the grip on the **real
  DNA-contacting residues** (not whole subunits) and activates it only during the H→P descent (physically motivated).
- Still needed: an **emergent** grip = explicit protein–DNA electrostatics/H-bonds (explicit solvent + ions) or a
  DNA-coupled coarse-grained model with a nucleotide-state-dependent contact map — where the DNA moves because of
  the force field, not a spring.

**(iv) Implicit solvent, no ATP·Mg, no real grip residues.**
- Done here: grip now on **real DNA-contact residues** (integrated run) — the "no real grip residues" part is addressed.
- Still needed: **explicit TIP3P + Mg²⁺/K⁺** (implicit GB screens the backbone electrostatics that ARE the grip);
  **ATP·Mg parameters** (GAFF/AM1-BCC) so the drive is nucleotide-state-driven, not a CV; the earlier explicit box
  NaN'd at 1.16 M atoms → **trim DNA to ~24 bp** + staged equilibration + a CUDA build first.

## 3. The defensible paper (what to actually write)
1. **The addressable single-chain gp16** — design, cross-predictor ring closure, designed coupling. (the anchor)
2. **A controlled in-silico platform** the single chain uniquely enables: driven P→H with retained coupling; the
   reproducible **graded, unequal per-subunit treads** (geometry); conformational-cycling-alone-does-not-translocate.
3. **A mechanistic framework + falsifiable predictions** (clearly labelled as model/prediction, not proof):
   - 3D-MINFLUX: per-subunit z(t) is a **travelling wave (sequential) vs in-phase (concerted)**; DNA **substeps vs
     one step**; the **special subunit is a dynamic phase-setter** (positionable on the single chain).
   - Optical tweezers: **unequal substep sizes**; a grip-weakening mutation shrinks substeps without changing the
     conformational cycle (structure-unchanged / translocation-changed double readout).
   - smFRET: R146→WalkerA stays engaged through the whole cycle in the single chain.
4. **Let the experiments arbitrate.** The paper's value is the platform + the sharp, testable predictions.

## 4. To-do before any MECHANISM claim (ranked)
1. Integrated **helical↔planar + emergent grip** in **explicit solvent + Mg·ATP** (trim DNA ~24 bp; CUDA).
2. Converged **metadynamics/umbrella** (soft kcv, difference-RMSD/string CV) for the landscape — replicates + error bars.
3. **AWSEM + 3SPN.2C** DNA-coupled model with a nucleotide-state contact map → emergent per-bp register + concerted-
   vs-sequential from committor/unbiased dynamics (the CV can't decide this).
4. Experimental cross-checks: MINFLUX per-subunit z(t), OT substep sizes, smFRET coupling — the real arbiters.
