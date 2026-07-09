# AlphaFold3 linker sweep — the third-predictor verdict (2026-07-08)

**Question:** Claude Science found L30/L34 close under Boltz-2 + OpenFold3 (the two NIM predictors). Does
AlphaFold3 — the strictest predictor, which gave 1/5 on L40 — confirm it? User ran the full AF3 sweep
(L20–L50, 5 models each). Scored here with the validated designed-sequential M2 (reproduce/score_m2.py).

## AlphaFold3 result — scrambles at EVERY linker length

M2 = designed-sequential interfaces engaged / 5, per AF3 model (all `sequential_consistent: NO` unless noted):

| linker | model0 | model1 | model2 | model3 | model4 | AF3 verdict |
|---|---|---|---|---|---|---|
| L20 | **5/5 (correct)** | 1/5 | 2/5 | 2/5 | 2/5 | 1 of 5 models correct — not robust |
| L25 | 1/5 | 1/5 | 1/5 | 1/5 | 1/5 | scrambled |
| **L30** | 0/5 | 0/5 | 0/5 | 0/5 | 0/5 | **scrambled (worst)** |
| **L34** | 2/5 | 2/5 | 2/5 | 2/5 | 2/5 | **scrambled** |
| L36 | 2/5 | 2/5 | 2/5 | 2/5 | 2/5 | scrambled |
| L50 | 2/5 | 2/5 | 2/5 | 2/5 | 2/5 | scrambled |

**AlphaFold3 never robustly forms the designed ring at any GS-linker length — including the L30/L34 that
Boltz-2 + OpenFold3 both closed.** Every AF3 model is a compact ring (radius_CV 0.00, planar) with the
copies threaded in the wrong order.

## Combined three-predictor picture (single-chain pentamer, M2 designed-sequential)

| linker | Boltz-2 (tiled MSA) | OpenFold3 (tiled MSA) | AlphaFold3 (own MSA) | 3-predictor robust? |
|---|---|---|---|---|
| L30 | 5/5 (5 seeds) | 5/5 | **0/5** | **NO** |
| L34 | 5/5 | 5/5 | **2/5** | **NO** |
| L40 | 5/5 (4/5 seeds) | 2/5 | 1/5 | no |
| L50 | 3/3 | 1/5 | 2/5 | no |

## Honest interpretation

1. **The fully single-chain GS-linker pentamer is NOT robust across all three predictors at any tested
   length.** The "linker fix" (L30/L34) was a *two-predictor* result (Boltz + OpenFold3); AlphaFold3 does
   not confirm it. We cannot declare the single chain solved.
2. **Caveat (fair to state):** AlphaFold3-server builds its **own** MSA on the full tandem-repeat sequence —
   it cannot be given our tiled core MSA (the AF3 server does not accept custom MSAs). Boltz/OF3's L30/L34
   closure used the tiled MSA. So part of AF3's pessimism *may* be the messier tandem MSA.
   **But the predictor-dependence is real, not purely MSA:** OpenFold3 with the *same tiled MSA as Boltz*
   still scrambled at L40 (2/5). So even matched-MSA, predictors disagree.
3. **Consequence:** to get a fully-single-chain construct that is robust under all predictors (a stronger,
   more determined design), escalate beyond a passive linker → **generative (RFdiffusion directional
   connector)** and/or **domain-recombination (chimera, swap the CTD)**. The conservation ladder's cheap
   rung (linker) did not survive the third predictor.
4. **B1 (multi-chain covalent dimer + 3 WT) remains the safest robust addressable construct** — it avoids
   the single-chain tandem-threading ambiguity entirely. Worth an AF3 check to make it three-predictor-clean.

## Ligand-state run — native ring + 5 ATP + 5 Mg (a clean win)

AF3 native ring with ATP+Mg loaded at all 5 sites → **all 5 R146 arginine fingers engage tightly and
symmetrically: 3.14–3.23 Å, interface pLDDT ~81** (vs the apo predicted ring's symmetric ~6.6 Å).

- **ATP tightens every finger** from the apo ~6.6 Å to ~3.2 Å — the arginine finger senses the neighbour's
  nucleotide (consistent with the ClpX central-coupler mechanism, Sosa et al. 2026).
- **Symmetric, not a lock-washer:** a *fully-loaded* ring (all 5 ATP) is symmetric. This **resolves the
  7JQQ mixed-occupancy caveat**: 7JQQ's asymmetry (tight at some interfaces, open seam) comes from *partial*
  ATP occupancy; fully-loaded = symmetric tight, apo = symmetric loose.

This is the ligand-state (nucleotide) validation depth, and it is clean.

Structures: `outputs/structures/af3_sweep/{L20,L25,L30,L34,L36,L50,native_ring_ATP_Mg}/` (5 models + confidences
each; 29 MB PAE full_data left in ~/Downloads). Cross-checked by three predictors, not validated; apo except the ATP/Mg run.
