# BioEmu run log — cp233 1750aa ring + gp16 domain (2026-07-10)

_A neutral record of what was run, what was measured, and the numbers. Interpretation is left open
(the divergence between these samples and the other structure/MD evidence is noted, not resolved).
Setup/how-to: `../../gcp_pipeline/BIOEMU_SETUP.md`. Visual: the artifact linked below._

## What was run
BioEmu (Microsoft, *Science* 2025) samples an equilibrium conformational ensemble of a protein
**monomer** from sequence.

| system | input | samples | compute |
|---|---|---|---|
| cp233 (1750 aa single chain) | tiled block-diagonal a3m `outputs/directed_evolution/folds/cp233_WTrep.tiled.a3m` | 10 | Spot A100-40GB, `bioemu[cuda]`, GPU embedding |
| gp16 monomer (327 aa) | `pipelines/tiled_msa_fold/gp16_core.a3m` | 10 | same VM |

- Settings: `--batch_size_100 350` (avoids the long-sequence batch underflow), `--filter_samples=False`
  (BioEmu's built-in physical filter removed **all** cp233 samples and then crashed on its own
  unphysical-save path — so filtering was turned off to obtain the raw samples).
- Version set (verified working): clean venv + `bioemu` + `numpy<2` + `tensorflow-cpu==2.18.0`.
- Wall time ~13 min / 10 samples; no GPU OOM. Credit acct `0161E2-E9BA4A-8F5A6A`, self-pay $0.

## What was measured
Per sampled frame (backbone/Cα), computed with `analyze_ensemble.py` (mdtraj):
- **total radius of gyration** (Rg);
- **whole-ring Cα-RMSD** to a folded reference (Kabsch);
- **per-subunit Cα-RMSD** and **per-subunit Rg** (the 1750 aa split into its 5 covalent copies:
  residues 1–342, 353–694, 705–1046, 1057–1398, 1409–1750);
- reference structures for comparison: `md/openmm_validation/trajectories/C/C_start.pdb` (folded design,
  used for MD), `.../inputs/A_apo.pdb` (native apo), `.../inputs/B_7jqq_helical.pdb` (7JQQ helical).

## The numbers

| quantity | folded design (ref) | native apo | 7JQQ helical | cp233 BioEmu (10 samples) | gp16 monomer (10) |
|---|---|---|---|---|---|
| total Rg | 3.68 nm | 3.54 nm | 3.84 nm | **5.85–8.64 nm** | 2.23–2.37 nm |
| whole-ring Cα-RMSD vs folded design | 0 | — | — | **45–76 Å** | — |
| per-subunit Cα-RMSD vs folded design | 0 | — | — | **9–28 Å** | — |
| per-subunit Rg | ~2.57 nm | — | — | 2.3–2.9 nm (a few 3.5–4.5) | — |
| adjacent-Cα breaks | — | — | — | 0.8% | 0.0% |

Per-frame table + 3D structure renders (folded reference vs most-compact vs most-expanded sample,
coloured by subunit) + Rg distribution are in the **artifact:**
https://claude.ai/code/artifact/79463dd3-f82a-4948-aa6c-d563ab726bad

### Note on the Rg reference
An initial pass compared Rg to a compact-globule scaling estimate (~3.3 nm); that is not the right
reference for a ring topology. The values above use the **actual folded structures** instead (folded
cp233 ring 3.68 nm; native apo 3.54; 7JQQ helical 3.84).

## Observations (descriptive)
- The **gp16 single-domain** samples are compact (Rg 2.29 nm) with no chain breaks.
- The **cp233** samples have larger total Rg (5.85–8.64 nm) than the three reference structures, high
  whole-ring RMSD (45–76 Å) to the folded design, and per-subunit RMSD of 9–28 Å; per-subunit Rg is
  mostly 2.3–2.9 nm with a few larger. In the 3D renders the five subunits appear separated rather than
  packed as in the folded reference.
- These cp233 samples **differ from** the fold-in-context predictions (AF3 + Boltz + OF3, previously
  scored closed 5/5) and from the unbiased OpenMM MD (which stayed closed). This divergence is recorded
  as an observation; it is **not resolved here.**

## Things NOT yet done / open questions (why the divergence is unresolved)
- Only 10 samples were drawn; convergence/coverage of the equilibrium distribution not assessed.
- BioEmu's `physical_steering` denoiser config was not used.
- Not tried: a 2-repeat sub-construct; conditioning/initializing on the AF3/Boltz folded ring; a
  cross-MSA (paired rather than block-diagonal) variant; longer sampling or a different seed set.
- Whether the block-diagonal tiled MSA supplies enough inter-subunit signal for BioEmu was not tested.
- No wet-lab or independent ensemble method has been run to adjudicate.

## Files
- `cp233_apo/{topology.pdb,samples.xtc}` — 10-frame ensemble of the 1750 aa ring.
- `gp16_monomer/{topology.pdb,samples.xtc}` — 10-frame ensemble of the 327 aa domain.
- `analyze_ensemble.py` — metrics + figures (needs mdtraj + matplotlib).
- `cp233_ensemble_analysis.html` — the artifact page.
- Load: `mdtraj.load('samples.xtc', top='topology.pdb')`. Backbone-only.
