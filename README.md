# gp16-single-chain-design — engineering repo & build history

Design and computational evaluation of a **single-chain, genetically-addressable φ29 gp16 DNA-packaging ring** — the pipelines, fold/MD scripts, result data, and the full commit history behind the work.

*Longfu Xu · Bustamante Lab, UC Berkeley / HHMI · Built with Claude Science + Claude Code*

> **Companion showcase repo →** [**`Rebuilding-motor-hackathon`**](https://github.com/longfuxu/Rebuilding-motor-hackathon)
> That repo is the narrated story (README + full report + figures). **This** repo is the engine room: the code that produced every result, the raw outputs, and the day-by-day commit history of how it was built.

---

## What this repo is

The project rebuilds nature's strongest molecular motor — a ring of five *identical* subunits that packs a virus's genome into its capsid — as **one covalent single chain** (`cp233`, 1,750 aa), so every subunit becomes individually addressable and the whole ring becomes computable. It then evaluates that design two ways: **static** (does it keep native's DNA grip, channel, and catalytic coupling?) and **dynamic** (does it run the packaging cycle, and what makes it translocate DNA?).

The scientific narrative, headline results, and figures live in the [companion repo](https://github.com/longfuxu/Rebuilding-motor-hackathon). This repo holds the machinery.

## Layout

| Path | What's in it |
|---|---|
| `gp16_design/pipelines/` | Design + analysis pipelines — circular-permutation construct generation, dead-seat / addressability variants, tiled-MSA folding of the 1,750-aa chain, P→H→P cycle drivers |
| `gp16_design/fold_scripts/` | Structure-prediction driver scripts (Boltz-2 / OpenFold3 / AlphaFold3 harnesses, NIM calls) |
| `gp16_design/md/` | Molecular dynamics — OpenMM predictor-independent stability check (`openmm_validation/`), driven / steered-targeted MD for the packaging cycle |
| `gp16_design/reproduce/` | Function-first metric scorers (grip on DNA · channel · catalytic coupling) and grounded functional-residue definitions |
| `gp16_design/gcp_pipeline/` | Ensemble-generation setup (BioEmu) on cloud GPU |
| `gp16_design/outputs/` | Result data — per-construct confidence summaries, coordination / dead-seat scans, cycle campaign, buildability atlas across 13 ring motors, ClpX retrospective validation |
| `gp16_design/cs_sessions/` | Claude Science session working outputs (folds, scores, handoffs) |

## Headline results (see the companion repo for the full write-up)

- `cp233` folds into the native ring: **3 independent predictors agree 5/5**, subunit **RMSD 1.80 Å / TM 0.94**, and it **threads dsDNA** through its channel like native.
- A catalytically-dead seat (E119Q) can be placed at any **chosen** position; a per-position scan shows all five seats are **statically equivalent** → the "special" subunit is a **timing** feature, not a structural one.
- Driven all-atom MD runs the full **P→H→P cycle reversibly** with coupling held; the conformational drive **alone does not translocate DNA** (honest 3-seed null).
- A **mechanochemical-ratchet model** shows sequential / hand-over-hand firing translocates DNA (~1.2 bp/cycle) while concerted barely moves — a clean **3D-MINFLUX discriminator**.
- The method **generalizes**: a two-branch buildability rule across **13 ring motors**, retrospectively validated on **ClpX** (passes the working motor 6/6, flags dead-coupler mutants 0/6).

## Built with Claude

- **Claude Code** orchestrated the compute — GPU jobs (Colab / GCP A100 / NVIDIA NIM) through hard-timeout harnesses, plus a multi-agent *score → cross-check → synthesize* workflow that verified every fold and surfaced the honest nulls.
- **Claude Science** ran the systematic folding / validation campaigns and structure rendering.
- Total compute **< $50 GPU**; the design/scoring logic is packaged as a reusable `protein-design-toolkit`.

## Reproduce

- Metric scorers + grounded functional residues: `gp16_design/reproduce/`
- Folding pipelines: `gp16_design/pipelines/tiled_msa_fold/` · `gp16_design/fold_scripts/`
- Predictor-independent MD check: `gp16_design/md/openmm_validation/`
- Large regenerable artifacts (trajectories, MSAs, predicted structures) are intentionally excluded from git; the scripts regenerate them.

---

*Computational design study. Contact the author before reusing results in a publication.*
