# cp233 is THREE-predictor-robust — the circular-permutation single chain works (2026-07-09)

The decisive third-predictor check. cp233 (circular permutation at res 233) single-chain gp16 pentamer,
scored by designed-sequential M2 + M1 (reproduce/score_m2.py, CP-aware --r146_incopy 255 --walker_incopy 133-140):

| predictor | cp233_int15_inter10 | cp233_int15_inter15 | cp233_int15_inter20 |
|---|---|---|---|
| Boltz-2 | 3/3 seeds, 5/5 | (NIM, prior) | 3/3 seeds, 5/5 |
| OpenFold3 | 5/5 | (NIM, prior) | 5/5 |
| **AlphaFold3** | **5/5 all 5 models, seq YES** | **4/5 models 5/5, 1 model 3/5; all 5 seq YES** | **5/5 all 5 models, seq YES** |

**All three independent predictors close cp233 into the correct designed cyclic ring (M2 5/5, M1
sequential_consistent = YES), unanimously.** This is the opposite of the N-C direct-linker single chain,
which AF3 rejected at every linker length (scrambled). It confirms the escalation-ladder conclusion:
direct fusion fails all three predictors → **circular permutation (cp233) succeeds under all three** →
diffusion not needed. cp233_int15_inter10 is the definitive single-chain lead (with B1_L40_E119Q the
multi-chain addressable lead). Cross-checked by three predictors, not validated (no wet lab); apo; geometric seat.
Structures: outputs/structures/af3_sweep/cp233/{inter10,inter15,inter20}/.

### AF3 inter15 added 2026-07-09 (completes the AF3 inter-linker sweep 10/15/20)
AlphaFold3 on `cp233_int15_inter15` (1770 aa; copies `A:1-342,A:358-699,A:715-1056,A:1072-1413,A:1429-1770`,
`--copy_start_res 1 --r146_incopy 255 --walker_incopy 133-140`): **models 1–4 close 5/5** (R146→Walker-A
7.2–7.6 Å), **model_0 is marginal at 3/5** (all five fingers 7.4–8.2 Å, two just over the 8 Å cutoff), and
**all 5 models are M1 sequential_consistent = YES** (compact planar ring, radius_CV 0.00). So inter15 also
passes the AF3 gate (majority 5/5, correct register everywhere) but is slightly softer than inter10/inter20
(whose fingers sit tighter). **inter10 remains the tightest/committed lead; inter-linker length 10 ≥ 15 ≈ 20
for AF3 robustness.** Global pTM 0.39 (not used — sequential M2 is the metric). Structures + summaries in
`outputs/structures/af3_sweep/cp233/inter15/` (140 MB `full_data` JSONs left in `~/Downloads`, not committed).
