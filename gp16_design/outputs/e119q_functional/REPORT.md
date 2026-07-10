# E119Q dead-seat AF3 folds — assembly tolerance confirmed (3rd predictor, ± ligand) 2026-07-09

The user ran the E119Q dead-seat series on AlphaFold3 (the 3rd independent predictor), apo AND with ATP·Mg. Scored handedness-robust M2 (cp233 landmarks R146=255, Walker-A=133-140, copies of 342).

| construct | best M2 | verdict |
|---|---|---|
| cp233_WT (apo) | **5/5** | baseline closed |
| E119Q_1seat (apo) | **5/5** | dead seat tolerated |
| E119Q_5seat (apo) | **5/5** | all-5 dead tolerated |
| E119Q_1seat + ATP·Mg | **5/5** | tolerated with ligand |
| E119Q_5seat + ATP·Mg | **5/5** | tolerated with ligand |
| rfdiff_ring d1 | 2/5 (mirror) | consistent w/ d0/d2 — rfdiff ring doesn't close on AF3 |

**Result:** the genetically-addressable E119Q dead seat (1 or all 5 seats) keeps the ring fully closed and sequential — **5/5 across all three predictors (Boltz-2 tiled, OpenFold3, AlphaFold3) and both ligand states.** This is the assembly-tolerance leg of the addressability claim, now 3-predictor + ligand-confirmed.

**Pending (agent died on an API error, re-run when stable):** the deeper *catalytic-locality* analysis — from the +ATP·Mg folds, per-seat Mg/attacking-water geometry to test whether E119Q locally disrupts catalysis *at its own seat* while neighbours stay intact (turning "assembly tolerance" into "position-specific catalytic knockout"). The folds are in `af3_2026_07_09/*_atp_mg/`.
