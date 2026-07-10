# BioEmu pilot — cp233 1750aa ring vs gp16 domain (2026-07-10)

_First BioEmu run for this project. BioEmu (Microsoft, Science 2025) samples the equilibrium
conformational ensemble of a protein MONOMER from sequence. Setup/how-to: `../../gcp_pipeline/BIOEMU_SETUP.md`._

## What was run
- **cp233 (1750 aa single-chain ring)** — the main design, fed our tiled block-diagonal a3m
  (`outputs/directed_evolution/folds/cp233_WTrep.tiled.a3m`), 10 samples.
- **gp16 monomer (327 aa single domain)** — control, `pipelines/tiled_msa_fold/gp16_core.a3m`, 10 samples.
- Compute: Spot **A100-40GB** (us-central1-a), `bioemu[cuda]`, GPU embedding, **~13 min / 10 samples, no OOM**.
  Credit acct `0161E2-E9BA4A-8F5A6A`, self-pay $0. `--filter_samples=False` (see gotcha below).

## Result — BioEmu folds the gp16 DOMAIN but NOT the 1750 aa artificial ring

**Visual + interactive write-up (artifact):** https://claude.ai/code/artifact/79463dd3-f82a-4948-aa6c-d563ab726bad
Reproduce metrics/figures: `analyze_ensemble.py` (needs mdtraj + matplotlib).

Judged by FOUR geometric measures against the **actual folded design** (`md/openmm_validation/trajectories/C/C_start.pdb`),
not Rg alone:

| measure | folded design (ref) | BioEmu cp233 samples | reading |
|---|---|---|---|
| total Rg | **3.68 nm** (native apo 3.54; 7JQQ helical 3.84) | **5.85–8.64 nm** | expanded beyond even the open helical state |
| whole-ring Cα-RMSD | 0 | **45–76 Å** | ring assembly destroyed |
| per-subunit Cα-RMSD | 0 | **9–28 Å** (a correct ~342-aa domain is <5 Å) | individual domains mis-folded too |
| per-subunit Rg | ~2.57 nm | 2.3–2.9 nm (a few 3.5–4.5) | most domains stay as compact *blobs*, but blob ≠ folded |
| gp16 monomer control (327 aa) | — | **Rg 2.29 nm, 0 breaks** | **folds fine ✓** |

- **Rg-reference correction (important):** an initial pass mis-compared Rg to a compact-*globule* estimate
  (~3.3 nm) — wrong reference for a ring. Recomputed against the real folded structures, a folded cp233 ring
  is only **3.68 nm** (native 3.54; even the open 7JQQ helical state 3.84). The samples (5.85–8.64 nm) sit far
  past all of them → not "large because it's a ring", and not a helical state — genuinely expanded.
- The five subunits are **strung apart, never assembled into a ring**, and the domains themselves are mis-folded
  (per-subunit RMSD 9–28 Å). BioEmu's samples are **not the folded design in any respect.**
- This **contradicts every other line of evidence** (AF3 + Boltz + OF3 all 5/5 closed; unbiased MD stays closed).
  Given that weight, the honest reading is **BioEmu fails to fold this construct — not that the design is
  unfolded.** The intact, compact gp16-**domain** control (327 aa → Rg 2.29 nm) confirms the pipeline/MSA are
  sound; the failure is specific to the 1750 aa 5-repeat ring.

## Why (mechanistic)
The ring's fold is set by the **inter-subunit interfaces**, whose signal is **absent from the block-diagonal
tiled MSA** (each repeat's co-evolution is independent; there is no cross-repeat coupling). AF3/Boltz still
assemble the ring via their structure module, but **BioEmu's diffusion prior — trained on natural monomers —
pulls a 1750 aa artificial 5-repeat chain toward generic expanded states.** The construct is outside its
training distribution.

## Implication (for paper / next steps)
- **Sequence-only equilibrium emulators (BioEmu) do NOT reliably fold mega artificial single-chain rings.**
  A clean, honest **methods limitation + future direction** (do NOT headline it in the hackathon).
- The ensemble/dynamics route for these constructs needs **structure-conditioned / steered methods**
  (start from the known folded structure, don't generate de novo) — exactly the hybrid architecture in
  `../../FUTURE_DIRECTIONS_ENSEMBLE_AND_NONEQ.md` (ligand-aware states → steered/targeted MD → MSM cycle).
- The gp16-**domain** ensemble (`gp16_monomer/`) IS usable — BioEmu captures single-domain gp16 dynamics.
- Untried levers if revisiting BioEmu on the ring: `--denoiser_config …/physical_steering.yaml`;
  smaller sub-constructs (2-repeat); or conditioning on the AF3/Boltz ring as a template (if supported).

## Files
- `cp233_apo/{topology.pdb,samples.xtc}` — 10-frame ensemble of the 1750 aa ring (expanded; the negative result).
- `gp16_monomer/{topology.pdb,samples.xtc}` — 10-frame ensemble of the 327 aa domain (folded; usable).
- Load: `mdtraj.load('samples.xtc', top='topology.pdb')`. Backbone-only.
