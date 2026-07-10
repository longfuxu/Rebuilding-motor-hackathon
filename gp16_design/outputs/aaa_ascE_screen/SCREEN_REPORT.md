# AAA+/ASCE single-chain ring screen — topology rule, three-method closure, and a generative winner

**Scope.** Extend the buildability atlas to a 16-protein AAA+/ASCE panel; for the foldable
subset run all three single-chain topologies (N→C direct fusion / circular permutation / RFdiffusion
generative connector); score each with the handedness-robust M1/M2 coupler metrics under a **tiled
MSA** and gate with **M3** (translocation channel clear + functional residues undisturbed); cross-check
survivors on a second predictor (OpenFold3). Goal: find which method closes each ring, and in
particular **≥1 RFdiffusion generative ring that closes** for wet-lab.

All folds used the free Boltz-2 / OpenFold3 NIMs with tiled (block-diagonal) MSAs. Global pTM/pLDDT
were never used as a closure criterion. Nothing here is called "validated" — these are *in-silico*
predictions across two independent structure predictors.

---

## 1. Extended panel + topology descriptors (16 proteins)

`panel_descriptors.csv` — descriptors computed from PDB biological assemblies (residues/termini read
from structure, never guessed). Three new verified homohexamers added to the atlas: **ClpB (6OAX),
NSF (3J94), LonA (6ON2)**.

**Two-branch rule** (reproduces 3/3 on the validated anchors gp16→CP, ClpX→direct, gp17→direct):

```
if C-terminus jams the channel        → circular permutation (CP)
elif head-to-tail C→N gap ≤ 38 Å      → direct N→C fusion
else                                  → RFdiffusion generative connector
```
Jam criterion (from the atlas classifier): substrate contact < 10 Å on either terminus **OR**
C-terminus relative radius < 0.40 (a merely inward N-terminus is *not* flagged).

**Result:** only **gp16** routes to CP — its C-terminus (res330) sits 4 Å from the translocated
substrate and 6.2 Å from the DNA channel, so a naive fusion would foul the pore. All 15 other panel
members have peripheral C-termini and route to **direct** under the two-branch rule. (A deeper
gap-only tree would split the 9 members with gap > 38 Å into a "diffusion" bin, but the validated
determinant of CP is channel-jamming, not raw gap — gp16 is the only jammer.)

---

## 2. Three-method closure on the foldable subset

NIM single-chain fold ceiling ≈ 1900 aa. Full hexameric single chains exceed it for most motor-domain
rings (ClpX 2378, FtsK 2598, p97 4692, …); the descriptor work covers all 16, the **fold** comparison
is tractable on **Rho (1590–1650 aa)** and **T7 gp4 (1878 aa)**, plus the prior gp16/gp17/ClpX anchors.

`screen_3method_results.csv` — master per-protein × three-method table. Which method closes:

| protein | oligomer | closing method | notes |
|---|---|---|---|
| *gp16* | pentamer | **CP** (cp233) | direct fouls DNA channel; CP off-channel cut; AF3 5/5 (3-predictor). RFdiffusion salami = documented negative (Boltz 0/5, AF3 split). |
| *gp17* | pentamer | **direct** | known-good anchor. |
| ClpX | hexamer | **direct** | geometry ruling (gap ~20 Å, termini off-pore; CP strictly worse). Full 6-copy chain exceeds NIM ceiling. |
| **Rho** | hexamer | **direct** + **RFdiffusion** | direct M2 6/6, M3-clean, OF3-confirmed. Two generative rings also close on both predictors (§3). CP390 closes M2 but fails M3. |
| **T7 gp4** | hexamer | **direct** | M2 6/6 (mirror winding — chirality to confirm); CP490 disrupts the coupler (1/6). |

### Rho three methods (`rho_3method_scores.csv`)
- **direct** (6× motor res175-414 + (GGGGS)×8): Boltz M2 **6/6**, M1 compact/sequential, M3 clean
  (pore 7.0 Å admits ssRNA; R366 within 0.6 Å of native). OpenFold3 **6/6** → **cross-predictor consistent.**
- **CP390**: Boltz M2 6/6 **but M3 FLAG** — arginine finger R366 pulled ~20 Å inward from its native
  radial position. This is exactly the failure mode the M3 gate exists to catch: an apo fold can satisfy
  the coupler contact yet distort the active site. **Rejected.**
- **CP330**: M2 0/6 — coupler not engaged. Fail.

### T7 gp4 (`t7gp4_scores.csv`) — a third topology class (SF4 helicase)
- **direct** (6× helicase res262-549 + (GGGGS)×6): M2 **6/6**, compact sequential ring — the coupler
  closes, matching the two-branch prediction (direct). Handedness is *mirror* (reverse winding); the
  physiological chirality should be confirmed before build.
- **CP490**: M2 1/6 — the permutation breaks the interface. Fail.

---

## 3. Generative ring — the headline result

`generative_summary.json`, `generative_winners_wetlab.json`.

Prior state: **no RFdiffusion ring had ever closed.** The gp16 "salami" backbones (adjacent-subunit
motifs frozen, connector diffused) folded to M2 0/5 even under tiled MSA — arginine fingers 22–57 Å
off, and AF3 split — a confirmed negative.

This campaign ran the full loop on the **Rho** salami backbones (10 RFdiffusion backbones already in
the repo): ProteinMPNN connector-only sequence design (motor motifs frozen, model v_48_020) → tile
6× native motor + 5× designed connector → tiled-MSA Boltz-2 fold → M1/M2/M3 → OpenFold3 cross-check.

**Outcome — first generative rings to close and survive the full gate:**

| backbone | Boltz M2 | OF3 M2 | M3 channel | verdict |
|---|---|---|---|---|
| **rhosal_L30_40_0_d0** | 6/6 designed | **6/6 designed** | clean | **CLOSES (cross-predictor confirmed)** |
| **rhosal_L40_52_2_d0** | 6/6 designed | **6/6 designed** | clean | **CLOSES (cross-predictor confirmed)** |
| rhosal_L30_40_1_d0 | 6/6 designed | 3/6 mirror | clean | Boltz-only (predictor-split) |
| rhosal_L30_40_2_d0 | 6/6 designed | 3/6 mirror | clean | Boltz-only (predictor-split) |
| rhosal_L40_52_0_d0 | 6/6 mirror | — | clean | M2-only (wrong winding) |
| rhosal_L40_52_1_d0 | 3/6 scrambled | — | — | fail |

All M3-clean closers keep the pore open (radius 6.3–7.4 Å, admits ssRNA), place Walker-A K184 and
arginine-finger R366 at native-like radial positions, and hold the designed connector outside the
lumen (radial 18.8–24.6 Å — no channel intrusion).

**Two generative Rho rings (`rhosal_L30_40_0_d0`, `rhosal_L40_52_2_d0`) close on both Boltz-2 and
OpenFold3 with the channel clean** — the first generative designs to pass M1 + M2 + M3 + cross-predictor.
Two further Boltz closers are predictor-split (the same pattern the AF3 index records for the gp16
generative rings), i.e. they should be treated as unconfirmed.

---

## 4. Does the topology rule extend, or is it falsified?

**Extends.** Across every folded system the closing method matched the two-branch prediction:
- gp16 (only jammer) → CP closes, direct fails on the channel. ✓
- gp17 / ClpX / Rho / T7 gp4 (peripheral C-termini) → direct closes. ✓
- Rho additionally admits a generative solution; T7 gp4 (SF4 helicase, a third topology class beyond
  the ASCE gp16/gp17 and RecA-like Rho) obeys the rule (direct closes, CP fails).

**The M3 gate is load-bearing.** Rho CP390 shows an apo M2 6/6 ring that the M3 channel gate rejects
(R366 displaced 20 Å). Apo single-chain closure alone is not sufficient evidence of a functional ring —
exactly the atlas cautionary lesson, now demonstrated on a new system.

---

## 5. Wet-lab-ready generative winner

`generative_winners_wetlab.json` carries the full single-chain ring sequences. Recommended lead:

**`rhosal_L30_40_0_d0`** — *E. coli* Rho (hexameric single-chain RNA translocase ring)
- Architecture: 6× native Rho motor (res 175-414) joined by 5× identical MPNN-designed connector.
- Designed connector (32 aa): `GLPKEEFYEWFEKHIDEVKFERIESKVVFDLI`
- Ring construct: 1600 aa, single chain.
- Predicted: M1 compact/sequential; M2 6/6 arginine fingers engaged on **both** Boltz-2 and OpenFold3
  (designed handedness); M3 channel open (pore ≥ 6.8 Å radius, admits ssRNA), K184/R366 native-like,
  connector outside the lumen.

Second candidate **`rhosal_L40_52_2_d0`** (42-aa connector `GLPKEEFYKWFEENIDTVFDGFKKEIGAKAIGFYPLVINRVL`,
1650-aa ring) — same cross-predictor + M3 status; a good independent backbone for parallel expression.

Suggested wet-lab checks: express the single-chain ring; SEC/negative-stain EM for hexameric
closure; ATPase activity ± poly-C RNA (Rho is RNA-dependent); RNA-stimulated ATP hydrolysis as the
functional coupler readout. The connector is the only non-native element — a scrambled-connector
control isolates its contribution.

---

## Files
- `panel_descriptors.csv/.json` — 16-protein descriptor table + two-branch/tree predictions
- `screen_3method_results.csv/.json` — master per-protein × three-method M1/M2/M3 closure table
- `rho_3method_scores.csv/.json` — Rho direct/CP/generative full scores incl. OF3 cross-predictor
- `t7gp4_scores.csv` — T7 gp4 direct/CP scores
- `generative_summary.json`, `generative_winners_wetlab.json` — generative outcome + winner sequences
- `screen_comparison.png` — 3-panel comparison figure
- `score_m3_rho.py`, `extend_descriptors.py` — scorers
- structures: `rho_folds/`, `rho_gen_folds/`, `of3_*.cif`, `t7gp4_folds/`
