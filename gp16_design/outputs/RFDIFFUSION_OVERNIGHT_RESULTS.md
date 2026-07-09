# RFdiffusion overnight run — de-novo connector (Track G / Mode B)

**Date:** 2026-07-09 (overnight). **Compute:** Colab Pro A100-40GB driven headless via the `colab` CLI
(first time this route was actually exercised end-to-end — prior sessions could not run RFdiffusion for
lack of a GPU). **Downstream scorer:** `reproduce/score_m2.py` (M1 ring geometry / M2 trans-R146 / M4).

## Objective (new work only)
Design a **de-novo connector** that fuses two ADJACENT gp16 subunits (A→B) into one chain while holding
both subunits — and therefore the R146(A)→Walker-A(B) trans interface (M2) — fixed. This is the generative
comparator to the circular-permutation lead **cp233**. Nothing already done was repeated; this is the
RFdiffusion backbone generation that the plan flagged as blocked-on-GPU.

## The design problem (measured from the motif)
- Input motif: `fold_inputs/generative/rfdiffusion/motif_AB_adjacent_subunits.pdb` — chains A, B, res 1–332
  each, carved from the closed native ring so their relative geometry IS the trans interface.
- **The connector must span a 53 Å gap** (A:330 Cα → B:4 Cα = 53.2 Å; A:330→B:1 = 50.5 Å). This large gap
  is exactly why the native N→C fusion cannot close and why cp233's permutation was invented. A rigid
  de-novo strut of ~30–50 residues is required; a floppy Gly-Ser linker would leave the register ambiguous.
- The M2 interface is present in the motif: R146(A) Cα → nearest Walker-A(B) Cα = 6.6 Å.

## Method
- ColabDesign RFdiffusion (`sokrypton/ColabDesign`, `rf/examples/diffusion.ipynb`), installed on the A100
  (dgl 2.4.0+cu124 works under Colab's newer torch). Contig `A4-330/<L>/B4-330` (fix both subunits,
  diffuse an L-residue connector into one chain), T=50, motif held fixed (verified RMSD 0.17 Å in a smoke run).
- Backbone layout is deterministic: subunit A = out res 1–327 (gp16 4–330), connector = middle,
  subunit B = last 327; connector length = N − 654. → `score_m2` copies `A:1-327,A:{N-326}-{N}` `--copy_start_res 4`.
- Geometric pre-filter (`reproduce/`-style, Cα-only since backbones are poly-Gly): connector compactness
  (Rg/residue), backbone clash (min connector↔motif Cα), and interface preservation (R146(A)→WalkerA(B) Cα).

## Backbones generated
Connector lengths sampled at 32 / 38 / 44 / 50 residues (spanning the plausible range for a 53 Å gap).

4 backbones generated and captured (1 per length), files in `outputs/rfdiffusion_modeB/backbones/`:

| design | connector L | N res | connector Rg/res | min clash (Å) | **R146(A)→WalkerA(B) Cα (Å)** | helix-frac | verdict |
|---|---|---|---|---|---|---|---|
| **sal_L50** | **50** | 704 | **0.448** | 3.76 | **6.96** | **0.49** | best — compact, ~50% structured rigid strut |
| sal_L44 | 44 | 698 | 0.539 | 3.75 | 6.95 | 0.27 | good |
| sal_L38 | 38 | 692 | 0.597 | 3.69 | 6.97 | 0.31 | ok |
| sal_L32 | 32 | 686 | 0.711 | 3.74 | 6.95 | 0.00 | strained extended coil (too short for 53 Å) |

**What this shows (real, preliminary):**
1. **RFdiffusion held the trans interface exactly** — R146(A)→Walker-A(B) Cα = 6.95–6.97 Å across all designs,
   vs 6.6 Å in the input motif (Δ ≈ 0.3 Å). The catalytic register that M2 measures is preserved by construction.
2. **No backbone clashes** (min connector↔motif Cα ≥ 3.69 Å) in any design.
3. **Clear length trend:** the 50-residue connector folds into the most **compact** (Rg/res 0.45) and most
   **structured** (helix-frac 0.49) strut, i.e. the most rigid candidate; the 32-residue connector is forced
   into a strained extended coil (Rg/res 0.71, 0% SS) because 32 residues cannot rigidly span 53 Å.
   → **sal_L50 is the lead de-novo connector backbone.**

Caveat: these are backbone-only (poly-Gly) geometry metrics. They establish the connector is *buildable* and
preserves the interface; they do NOT yet establish that a real sequence encodes it (that is the MPNN→fold→M2
gate below). `helix-frac` is a crude Cα(i)–Cα(i+3) proxy, not DSSP.

## Downstream — ProteinMPNN + Boltz-2 fold (DONE 2026-07-09)
**ProteinMPNN (local, soluble v_48_020, T=0.1):** threaded the native gp16 sequence onto each backbone's
motif, then designed **only the connector** with the entire motif (both subunits, incl. Walker-A/B, R146,
ATP pocket) frozen. Verified: catalytic residues R146/K30/E119 native in both subunits of every design,
motif byte-identical to native, connector fully redesigned. 16 designs (4/backbone) in
`mpnn_designs/connector_designs.fasta`.

**Boltz-2 NIM fold + `score_m2` (15 designs folded):** every design folds (pLDDT 0.44–0.53) but the isolated
**A→connector→B dimer gives M2 = 0/2** across all connector lengths and designs (`mpnn_designs/dimer_fold_score.json`).

**This is expected, not a failure of the connector.** The project already established (HANDOFF §5, B1 work)
that **the trans interface is a ≥3-subunit property** — even the *native-sequence* isolated 2-mer does not
reform it (Boltz 0/5); it needs both neighbours. So an isolated A-connector-B **dimer is the wrong test unit**
for M2. The dimer fold confirms the designed connector sequences are foldable; it cannot report the trans-interface metric.

**Correct next test (deferred, ~1700 aa like cp233):** tile the lead connector (sal_L50) 5× into a single-chain
native-order ring, fold that ring, and score sequential M2. Only the tiled ring can show whether the de-novo
connector closes with correct register — the true generative comparator to cp233. Scripts (`mpnn_pipeline.py`,
`nim_fold_score.py`) and the lead backbone/sequences are in this directory.

## Tiled ring — the decisive test (DONE 2026-07-09)
Tiled the lead connector **sal_L50** 5× into a single-chain native-order ring:
`5×[gp16 4-330] + 4×[designed L50 connector] = 1835 aa`. Folded all 4 L50 connector designs with
Boltz-2 NIM (single-seq) and scored sequential M2 (copies `A:1-327,A:378-704,A:755-1081,A:1132-1458,A:1509-1835`
`--copy_start_res 4`). Structures + scores in `rfdiffusion_modeB/tiled_ring/`.

| design | ring radius (Å) | M1 sequential | compact ring | **M2** | pLDDT |
|---|---|---|---|---|---|
| sal_L50_d0 | 27.4 | **YES** | True | **0/5** | 0.45 |
| sal_L50_d1 | 30.3 | **YES** | True | 0/5 | 0.45 |
| sal_L50_d2 | 29.6 | **YES** | True | 0/5 | 0.48 |
| sal_L50_d3 | 26.7 | **YES** | True | 0/5 | 0.44 |
| _cp233 (ref, tiled MSA)_ | _25.4_ | _YES_ | _True_ | _**5/5**_ | _—_ |

**Result — a real, mixed outcome:**
- **Topological success:** the de-novo connector, tiled 5×, folds into a **compact pentameric ring with the
  correct designed sequential order** (M1 sequential_consistent = YES, compact_ring = True, all 4 designs). This
  is exactly what native-order *passive* fusion could NOT do (it scrambled or failed to close) — the generative
  connector **locks the ring register**. That is a genuine gain from generative design.
- **Catalytic shortfall:** but every ring is **~2–5 Å too loose radially** (27–30 Å vs cp233's 25.4 Å), so the
  R146→Walker-A trans interface never engages — **M2 = 0/5 across all four designs**. cp233 (circular permutation)
  gets **M2 5/5**. So on the load-bearing catalytic metric, **cp233 remains the lead; the de-novo connector does
  not (yet) beat it.**

**Honest conclusion (answers the plan's Milestone-G question):** generative de-novo design bought *register-locking
+ ring compactness* but **not catalytic-interface engagement beyond cp233** — consistent with the plan's expectation
that cp233 is the winner and the generative track is an optional comparator, now with concrete negative-on-M2 evidence.
Caveats: single-seq Boltz (cp233's 5/5 used a **tiled block-diagonal MSA** — a fair re-test needs the same tiled MSA,
which could tighten the ring; deferred), and the connector was optimised for one 2-subunit geometry so tiling
accumulates small errors. Next refinements if this track is pursued: tiled-MSA refold; a shorter/tighter connector
window; or partial-diffusion tightening (Mode A).

## Honest notes on the compute run
- The `colab` CLI route works but is **flaky for long unattended runs**: `colab exec`/`download` websockets
  intermittently hang past their timeout, and Colab sessions self-terminate on idle. Two A100 sessions were
  lost mid-run before a **local Python orchestrator with hard per-call timeouts + one-tarball capture +
  immediate `colab stop`** made capture reliable. Operational details recorded in the memory
  `colab-cli-overnight-gotchas.md`.
- This is a small, honestly-scoped pilot (backbone generation + geometry). It does NOT claim a validated
  ring; two-predictor M2 agreement on MPNN-designed sequences is the next gate (deferred, scripts ready).
