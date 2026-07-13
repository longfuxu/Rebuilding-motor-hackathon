# Per-position dead-seat scan — cp233 single chain, R146A at each of the 5 seats

One dead trans-finger (R146A) placed at each ring position, folded on Boltz-2 NIM (3 reps),
M2 = trans-R146→neighbour Walker-A <8 Å, handedness-robust (engaged_max of 3 reps). WT baseline = 5.

| construct | engaged_max | engaged (3 reps) | radius_CV | sequential | interpretation |
|---|---|---|---|---|---|
| WT | **5** | [5, 5, 5] | 0.00 | True | baseline (all fingers) |
| R146A_pos1 | **4** | [4, 4, 4] | 0.00 | True | one finger lost (local, 5→4 as expected) |
| R146A_pos2 | **4** | [4, 4, 4] | 0.00 | True | one finger lost (local, 5→4 as expected) |
| R146A_pos3 | **4** | [4, 4, 4] | 0.00 | True | one finger lost (local, 5→4 as expected) |
| R146A_pos4 | **4** | [4, 4, 4] | 0.00 | True | one finger lost (local, 5→4 as expected) |
| R146A_pos5 | **4** | [4, 4, 4] | 0.00 | True | one finger lost (local, 5→4 as expected) |

**Position dependence:** engaged_max across the 5 positions = {'R146A_pos1': 4, 'R146A_pos2': 4, 'R146A_pos3': 4, 'R146A_pos4': 4, 'R146A_pos5': 4}. Spread = **0**.
- spread 0 ⇒ position-INDEPENDENT (every seat loses exactly one finger, 5→4) — no 'special' seat by this readout.
- spread ≥1 ⇒ position-DEPENDENT — some seat's dead finger costs more/less coordination (a candidate special seat).

_Caveat: Boltz-2 static fold + M2 geometric proxy; predictor basin-bias (all cp233 fold planar);
3 reps only. This is a design-prioritisation prediction for the mixed-ring experiment, not proof._