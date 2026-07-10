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

## Catalytic-locality (ATP·Mg folds, per-seat) — 2026-07-09
AF3 placed **5 ATP + 5 Mg** (one per pocket) in both ligand folds. Per-seat distance from residue-119 side-chain (E: OE1/OE2; Q→OE1/NE2) to nearest Mg:

| construct | dead seat(s) Q→Mg | wt seat(s) E→Mg |
|---|---|---|
| E119Q_1seat +ATP·Mg | copy1 (Q) **4.2 Å** | copies 2–5 (E) 4.3–4.4 Å |
| E119Q_5seat +ATP·Mg | all (Q) **4.2 Å** | — |

**Finding (honest, and a cleaner story):** E119Q does **NOT** structurally displace the Mg or the pocket — the ground-state geometry is essentially identical to WT (Δ ≤ 0.2 Å). So the addressable dead seat is a **structurally-silent, position-specific catalytic knockout**: it perturbs neither assembly (M2 5/5), nor the neighbours, nor the visible ATP·Mg pocket. The catalytic loss is **chemical** — E119 is the base that activates the lytic water; Q cannot — which a *static* structure predictor cannot capture. This is exactly what you want in an addressable perturbation (minimal structural side-effects); to *see* the catalytic effect requires QM/MM or the wet ATPase/single-molecule assay. It also revises our earlier hypothesis (we expected local Mg displacement) — honestly, there is none.

## WT + ATP·Mg baseline landed (2026-07-09) — confirms structurally-silent
cp233_WT + ATP·Mg (the missing baseline) folded on AF3: **M2 5/5**, 5 ATP + 5 Mg placed, per-seat **E119→Mg 4.1–4.2 Å** at all five WT pockets. This equals the E119Q Q→Mg distance (4.2 Å) — so **WT and E119Q hold Mg identically**; the dead seat does not displace it. Definitive: E119Q is a structurally-silent, position-addressable catalytic knockout (the loss is chemical, not structural).
