# cp233 is THREE-predictor-robust — the circular-permutation single chain works (2026-07-09)

The decisive third-predictor check. cp233 (circular permutation at res 233) single-chain gp16 pentamer,
scored by designed-sequential M2 + M1 (reproduce/score_m2.py, CP-aware --r146_incopy 255 --walker_incopy 133-140):

| predictor | cp233_int15_inter10 | cp233_int15_inter20 |
|---|---|---|
| Boltz-2 | 3/3 seeds, 5/5 | 3/3 seeds, 5/5 |
| OpenFold3 | 5/5 | 5/5 |
| **AlphaFold3** | **5/5 all 5 models, sequential YES** | **5/5 all 5 models, sequential YES** |

**All three independent predictors close cp233 into the correct designed cyclic ring (M2 5/5, M1
sequential_consistent = YES), unanimously.** This is the opposite of the N-C direct-linker single chain,
which AF3 rejected at every linker length (scrambled). It confirms the escalation-ladder conclusion:
direct fusion fails all three predictors → **circular permutation (cp233) succeeds under all three** →
diffusion not needed. cp233_int15_inter10 is the definitive single-chain lead (with B1_L40_E119Q the
multi-chain addressable lead). Cross-checked by three predictors, not validated (no wet lab); apo; geometric seat.
Structures: outputs/structures/af3_sweep/cp233/{inter10,inter20}/.
