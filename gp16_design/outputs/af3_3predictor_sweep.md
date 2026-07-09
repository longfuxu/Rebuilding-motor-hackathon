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

## Ligand-state run — native ring + 5 ATP + 5 Mg (CORRECTED — an unreliable prediction, not a win)

**Retraction of an earlier over-claim.** I first reported this as "ATP tightens the fingers, resolves the
7JQQ caveat." That was wrong — I only measured the R146→neighbour distance, not the ring's open/closed
(planar vs helical) geometry. Measuring the actual ring geometry:

| state | R146 | radius | planarity_rms | ring shape |
|---|---|---|---|---|
| apo (Boltz, **no template**) | ~6.6 Å symmetric | 26.9 Å | **0.05 Å** | closed **planar** |
| ATP/Mg (AF3, **template ON**) | ~3.2 Å symmetric | 27.9 Å | **0.01 Å** | closed **planar** |

**Both are flat, symmetric, closed rings** — the ATP run did **not** produce the open/helical (lock-washer)
state expected for the ATP-bound motor (the state captured experimentally by 7JQQ). The only change AF3 made
was tightening the fingers within a still-closed planar ring. Why AF3 is unreliable here:
1. **The ATP job used `useStructureTemplate: true`** → biased toward a closed symmetric template.
2. **More fundamentally, Boltz/OF3/AF3 give a single static structure and do not model the ATP-driven
   allosteric planar↔helical transition** — they dock ATP into whatever ring they already build.

**Honest conclusion:** the ligand-state / switch-position conformational question is **beyond what these
structure predictors can answer**. The experimental answer already exists — **7JQQ is the ATP-analog-bound
helical/asymmetric state**. The predictors failing to reproduce it is a *predictor limitation*, not biology.
Correct note on apo: the apo ring IS closed planar (planarity 0.05 Å) — the ~6.6 Å is just the loosely-engaged
finger distance in the closed state, **not** an open ring. Verifying on AF3 would require `useStructureTemplate:
false` and an apo/ADP/ATP comparison, but expect the predictors to still fail to open the ring.

Structures: `outputs/structures/af3_sweep/{L20,L25,L30,L34,L36,L50,native_ring_ATP_Mg}/` (5 models + confidences
each; 29 MB PAE full_data left in ~/Downloads). Cross-checked by three predictors, not validated; apo except the ATP/Mg run.
