# AF3 tie-breaker — single-chain (GGGGS)×8 gp16 pentamer

**Question:** does the fully single-chain pentamer form the DESIGNED cyclic ring? Boltz-2 said yes (5/5),
OpenFold3 said no (2/5, scrambled). AlphaFold3 (5 models, own MSA, seed 1) is the independent third vote.

**Verdict: AF3 breaks the tie AGAINST Boltz — the single chain scrambles.**

## Three-predictor M2 (designed-sequential trans-R146 → neighbour Walker-A, <8 Å; scored by reproduce/score_m2.py)

| Predictor | global pTM | M2 (designed ring) | topology |
|---|---|---|---|
| Boltz-2 (tiled MSA) | 0.48 | **5/5** @ 4.9–5.5 Å | designed cyclic ring closes |
| OpenFold3 (tiled MSA) | 0.58 | **2/5** (only 3→4, 5→1) | compact but scrambled |
| AlphaFold3 (5 models) | 0.49–0.51 | **1/5** (only 1→2), identical across all 5 | compact but scrambled; ring threads 1→2→5→4→3 |

**Consensus: 2 of 3 predictors (OpenFold3 + AlphaFold3) fail the designed ring; AF3 is unanimous across all
5 models. Boltz-2 is the lone optimist.** The (GGGGS)×8 single-chain pentamer does **not** robustly close
into the designed cyclic topology — the flexible linkers let the chain thread into a scrambled register
that buries every R146 against a *wrong* partner (~7 Å) while leaving the *designed* sequential interfaces
open at 53–56 Å.

## Why this matters (two load-bearing points)

1. **It hardens the metric-discipline argument.** All three predictors sit at global pTM 0.48–0.58 — a band
   that cannot tell the correct ring (Boltz) from the scrambled ones (OF3, AF3). Only the interface-resolved,
   **designed-sequential** M2 separates them. This is the cleanest possible evidence that global pTM is the
   wrong objective. (AF3 even reports `has_clash 0.0`, `fraction_disordered 0.09` — it "looks" clean.)
2. **It fires the generative-escalation trigger cleanly, not speculatively.** The passive-linker single chain
   is now majority-rejected across independent predictors → the RFdiffusion + LigandMPNN/ProteinMPNN
   "self-repair" roadmap (5-day goal) is warranted: replace the non-directional GS linker with a designed,
   directional connector that makes the correct register the only low-energy solution.

## Consequence for the paper / project

- **B1 (multi-chain covalent dimer + 3 WT) remains the robust position-addressable construct** (Boltz 15/15,
  Chai 47/50). Lead the paper with B1; present the fully single-chain pentamer as the harder case where
  predictors diverge → generative design is the principled next step (a *result*, not a failure).
- Cross-checked by three predictors, not validated; apo (no ATP/pRNA/DNA); geometric seat.

Structures + confidences: `outputs/structures/ladder/af3/` (5 models; 29 MB PAE full_data left in Downloads).
Boltz-2 pTM 0.48 / OF3 pTM 0.58 from the ladder-MSA run; AF3 pTM 0.49–0.51 from `summary_confidences_*.json`.
