# C1 result — ligand-conditioned structure prediction does NOT open the ring (2026-07-10)

_Aim C1 of `../../COMPUTATIONAL_AIMS_PHP_CYCLE.md`: does adding ATP to the planar apo ring shift it toward
the 7JQQ helical (H) state? Tested on two independent predictors (Boltz-2 NIM + AF3) and two constructs
(native pentamer + single-chain design), under apo / ATP / ADP / AGS._

## Readout
Per-subunit axial **spread** of the 5 subunit centroids (best-fit-plane normal): **planar ≈ 0.1 Å; the real
7JQQ helical staircase ≈ 4.8 Å** (calibrated from `md/openmm_validation/inputs/{A_apo,B_7jqq_helical}.pdb`).

| predictor · construct · state | axial spread (Å) | radius (Å) | verdict |
|---|---|---|---|
| **ref** native apo (A_apo) | 0.1 | 26.9 | planar |
| **ref** 7JQQ helical (B) | **4.8** | 29.1 | helical |
| AF3 native apo | 0.2 | 27.2 | planar |
| AF3 native **ATP** | 0.0 | 27.9 | planar |
| AF3 native ADP | 0.2 | 28.0 | planar |
| Boltz native apo | 0.25 | 27.3 | planar |
| Boltz native **ATP** | 0.08 | 27.4 | planar |
| Boltz native **AGS** (ATPγS) | 0.07 | 27.5 | planar |
| Boltz single-chain design apo | 0.23 | 25.5 | planar (design closes ✓) |
| Boltz single-chain design **ATP** | 0.02 | 25.8 | planar |
| Boltz native apo **+DNA** | 0.81 | 27.8 | planar (DNA barely perturbs) |
| Boltz native ATP **+DNA** | 0.16 | 27.4 | planar |
| Boltz native **AGS +DNA** (= 7JQQ conditions) | 0.22 | 27.4 | planar (still no helix) |
| Boltz single-chain ATP **+DNA** | 0.04 | 25.8 | planar |

## Finding
**Every structure-prediction condition returns a planar ring** (spread ~0.02–0.81 Å) — regardless of
nucleotide (apo/ATP/ADP/AGS), construct (native/single-chain), predictor (AF3/Boltz), **or DNA in the
channel**. Crucially, **AGS + DNA (the exact 7JQQ conditions) still predicts planar (0.22 Å)** — it does
not reproduce the helical H state (4.8 Å). Adding DNA at most nudges the apo ring to 0.81 Å (still planar).
**None reproduces the helical state, and none responds to ATP occupancy.**

Because two independent predictors agree, this is best read as a **predictor-basin bias** — AF3/Boltz are
trained toward the most stable/common conformation (the planar closed ring) and do not capture the
ATP-driven allosteric opening — **not** as evidence that ATP fails to open the ring. This is the
inconclusive-C1 branch anticipated in the aims doc.

## Consequence
- **Structure prediction cannot test the P→H→P transition.** The ligand-driven opening must be probed by
  **physics** → escalate to Aim C3, now framed as **coarse-grained / structure-based (dual-basin) modeling**
  (the P↔H transition is ~10–100 ms and ~35 Å — beyond all-atom MD; see the revised C3).
- Useful positive by-products: the **single-chain design closes to a planar ring** (spread 0.23 Å, radius
  25.5 Å — slightly tighter than native 27.3) with and without ATP, consistent with the closed-ring evidence.
- The planar predictions (`outputs/php_cycle/C1_boltz/*.cif`, `AF3/*.cif`) are usable as the **P endpoint /
  starting structures** for the C3 model.

## Files
- `C1_boltz/*.cif` + `C1_results.json` — Boltz-2 folds (native apo/atp/ags, single-chain apo/atp).
- `AF3/*.cif` — the user's AF3 native pentamer apo/atp/adp folds (scored here).
- `../../pipelines/tiled_msa_fold/` + the C1 driver `pipelines/php_cycle/c1_ligand_fold.py` (Boltz-2 NIM with
  ligands via CCD codes; scorer = per-subunit centroid axial spread / planarity / radius).
- Cross-predictor design note: run AF3 jobs in `AF3_JOBS.md` if more AF3 states are wanted.
