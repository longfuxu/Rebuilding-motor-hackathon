# Rebuilding Molecular Motors towards Programmable Viral DNA Packaging

**A single-chain, genetically-addressable φ29 gp16 DNA-packaging ring — designed *and* evaluated end-to-end, entirely in silico.**

*Longfu Xu, Ph.D · Bustamante Lab, UC Berkeley / HHMI · Built with Claude Science + Claude Code*
*"Built with Claude: Life Sciences" (Anthropic × Gladstone) — Research track.*

This README is the guide to the whole submission: where the idea came from, what we did, what we found, and where every material lives. Watch the 3-minute video for the story; browse this repo for the storyline; read the report if you want the rigorous detail.

---

## ▶ Start here
| | | |
|---|---|---|
| **Demo video** (3 min, narrated) | the whole story, start to finish | [**▶ Watch on YouTube**](https://youtu.be/ZxMc6qudB5Q) · in repo: [`gp16_demo.mp4`](gp16_demo.mp4) |
| **Deck** (10 slides, animated) | the same story as slides | [`gp16_deck.pptx`](gp16_deck.pptx) |
| **Full report** (13 pp, rigorous) | every claim, with figures and honest bounds | [`report.pdf`](report.pdf) |

---

## Where this came from

I am a wet-lab single-molecule biophysicist. I have spent my career at the bench measuring molecular motors one molecule at a time, and I had never done computational protein design — I had never written simulation code. This project is my first entirely in-silico work, and it happened because Claude Science and Claude Code let me design a motor, fold it, simulate it, and cross-check it, end to end, without leaving the reasoning to guesswork. The science below is the result; the fact that a bench scientist could do it at all is the point.

## What we did (the story)

Nature's strongest motor is a **ring of five *identical* subunits** — the φ29 gp16 ATPase — that packs a virus's genome into its shell against tens of piconewtons. That sameness is a wall: you **can't tune one subunit**, you **can't read the firing order or mechanism**, and the oligomer **isn't even computable**. We did two things, both in silico:

1. **DESIGN** — rebuild the ring as **one covalent single chain** (`cp233`, 1,750 aa), so every seat is individually addressable and the whole ring becomes computable.
2. **EVALUATE** it two ways — **static** (does it keep native's machinery: grip on DNA, channel, catalytic coupling?) and **dynamic** (does it run the packaging cycle, and what makes it translocate DNA?).

Along the way we did far more than a single design: we screened **20 circular-permutation cut sites**, tried **generative diffusion (RFdiffusion) linkers**, ran **in-silico directed evolution**, and tested the recipe across **13 ring motors** plus a blind retrospective on **ClpX**. The report tells all of it; the highlights are below.

## What we found (headline results)

- `cp233` folds into the native ring: **3 independent predictors agree 5/5** (Boltz-2 + OpenFold3 + AlphaFold3), subunit **RMSD 1.80 Å / TM 0.94**, and it **threads dsDNA** through its channel at the native contacts — designs that merely *close* the ring cannot.
- Why circular permutation: native termini are **~53 Å apart** and swing **~13 Å** per stroke, so direct fusion scrambles the register; permutation moves the seam to a short, low-stroke loop. A systematic screen found **residue 233** as the unique cut that spares the DNA-gripping C-domain (rivals close the ring but collapse the channel).
- A **catalytically-dead seat (E119Q)** can be placed at any **chosen** position; a per-position scan shows all five seats are **statically equivalent** → the "special" subunit is a **timing** feature, not a structural one.
- Driven all-atom MD: the ring runs the full **P→H→P cycle reversibly**, coupling held, with a graded unequal per-subunit staircase. The conformational drive **alone does not translocate DNA** (honest 3-seed null).
- A **mechanochemical-ratchet model** shows **sequential / hand-over-hand firing translocates DNA (~1.2 bp/cycle) while concerted barely moves** — a clean **3D-MINFLUX discriminator** (a travelling wave across the five subunits = sequential).
- The method **generalizes**: a two-branch buildability rule across **13 ring motors**, retrospectively validated on **ClpX** (passes the working motor 6/6, flags the dead-coupler mutants 0/6).

---

## Full manifest — everything in this submission

| Item | What it is | Where |
|---|---|---|
| **Demo video** | 3-min narrated walkthrough | [▶ YouTube](https://youtu.be/ZxMc6qudB5Q) · [`gp16_demo.mp4`](gp16_demo.mp4) |
| **Deck** | 10-slide animated deck | [`gp16_deck.pptx`](gp16_deck.pptx) |
| **Report** | 13-page rigorous report (figures + honest bounds) | [`report.pdf`](report.pdf) |
| **Figures** | the 8 figures used across deck and report | [`figs/`](figs/) |
| **Result records** | cycle campaign · ratchet · per-position scan · session log | [`results/`](results/) |
| **Code + full history** | every pipeline, fold/MD script, scorer, and raw output — plus the complete commit history of how this was built | [`gp16_design/`](gp16_design/) |

## Reading the code (`gp16_design/`)

| Path | What's in it |
|---|---|
| `gp16_design/pipelines/` | design + analysis pipelines — CP construct generation, addressable dead-seat variants, tiled-MSA folding of the 1,750-aa chain, P→H→P drivers |
| `gp16_design/fold_scripts/` | structure-prediction drivers (Boltz-2 / OpenFold3 / AlphaFold3 harnesses) |
| `gp16_design/md/` | molecular dynamics — OpenMM predictor-independent stability check, driven MD for the packaging cycle |
| `gp16_design/reproduce/` | function-first metric scorers + grounded 7JQQ functional residues |
| `gp16_design/outputs/` | result data — confidence summaries, coordination / dead-seat scans, cycle campaign, buildability atlas, ClpX validation |

The commit history is the day-by-day record of the build; browse the *Commits* tab to see it.

## Built with Claude

- **Claude Code** orchestrated the compute — GPU jobs (Colab / cloud A100 / NVIDIA NIM) through hard-timeout harnesses, plus a multi-agent *score → cross-check → synthesize* workflow that verified every fold and surfaced the honest nulls (DNA translocation, absolute ΔG).
- **Claude Science** ran the systematic folding / validation campaigns (Boltz-2 / OpenFold3 / AlphaFold3) and structure rendering.
- Total compute **< $50 GPU**; the design and scoring logic is packaged as a reusable protein-design toolkit.

## Honest bounds

Everything here is **computational** — cross-checked by ≥2 structure predictors + physics-based MD, and reported with explicit reliability. **Nothing is wet-lab-validated yet** — that is the stated next step. We steer by interface / channel **geometry, never global pTM**. The dynamics are driven models, not rates; the plain-cycle DNA translocation is an honest null; the coupling gain is gp16-specific. Full detail and every caveat are in [`report.pdf`](report.pdf). The single-chain design is a defined scaffold for the bench: single-molecule (optical tweezers · 3D-MINFLUX), ATPase biochemistry, and cryo-EM.
