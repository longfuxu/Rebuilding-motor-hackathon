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

| construct | Rg (mean) | compact ref | chain breaks | verdict |
|---|---|---|---|---|
| gp16 monomer (327 aa) | **2.29 nm** (2.23–2.37) | ~1.86 nm | 0.0% | **folded, compact ✓** |
| cp233 ring (1750 aa) | **6–9 nm** (5.85–8.64) | ~3.3 nm | 0.8% | **expanded, NOT a closed ring ✗** |

- cp233's BioEmu samples are **expanded/extended (Rg 2–3× too large), with some chain breaks + clashes** —
  they are **not** the compact closed ring.
- This **contradicts every other line of evidence** (AF3 + Boltz + OF3 all 5/5 closed; MD stays closed).
  Given that weight, the honest reading is **BioEmu fails to fold this construct — not that the design is
  actually expanded.** The intact, compact gp16-domain control confirms the pipeline/MSA are fine; the
  failure is specific to the 1750 aa 5-repeat ring.

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
