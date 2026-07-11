# AF3-server jobs — C1 ligand-conditioned P↔H test (parallel to the Boltz-2 run)

_Run these on **alphafoldserver.com** to cross-validate the Boltz-2 NIM results (Claude is running the
Boltz side). The question: **does ATP occupancy shift the gp16 ring from planar (P) toward helical (H)?**
Compare each job's ring geometry (planar disc vs helical staircase) — especially the **axial spread** of
the 5 subunits (planar ≈ small; 7JQQ-helical ≈ up to ~35 Å)._

## The native gp16 protein sequence (327 aa, from 7JQQ) — use this for every "native" job, 5 copies
```
SLFYNPQKMLSYDRILNFVIGARGIGKSYAMKVYPINRFIKYGEQFIYVRRYKPELAKVSNYFNDVAQEFPDHELVVKGRRFYIDGKLAGWAIPLSVWQSEKSNAYPNVSTIVFDEFIREKDNSNYIPNEVSALLNLMDTVFRNRERVRCICLSNAVSVVNPYFLFFNLVPDVNKRFNVYDDALIEIPDSLDFSSERRKTRFGRLIDGTEYGEMSLDNQFIGDSQVFIEKRSKDSKFVFSIVYNGFTLGVWVDVNQGLMYIDTAHDPSTKNVYTLTTDDLNENMMLITNYKNNYHLRKLASAFMNGYLRFDNQVIRNIAYELFRKMR
```

## Jobs (in the AF3 UI: add "Protein" with copies, then "Ligand"/"Ion", then "DNA")

| # | job | protein | ligands / ions | DNA | tests |
|---|-----|---------|----------------|-----|-------|
| 1 | **native_apo** | gp16 ×5 | none | none | baseline: does apo predict planar? |
| 2 | **native_ATP** | gp16 ×5 | ATP ×5 + MG ×5 | none | **does ATP open → helical?** (your Q1) |
| 3 | **native_AGS** | gp16 ×5 | AGS ×5 + MG ×5 | none | non-hydrolyzable (matches 7JQQ) → helical? |
| 4 | **native_ATP_DNA** | gp16 ×5 | ATP ×5 + MG ×5 | dsDNA (below) | ± DNA contrast (your Q3) |
| 5 | **native_apo_DNA** | gp16 ×5 | none | dsDNA (below) | ± DNA contrast, apo |

- **Ligand codes in AF3:** ATP = `ATP`, ATP-γ-S = `AGS`, ADP = `ADP`; **Ion:** `MG` (magnesium). Add 5 of
  each ligand/ion (AF3 lets you set a copy count or add repeatedly).
- **DNA** (jobs 4–5): add two DNA strands (AF3 pairs complementary strands into dsDNA):
  - strand 1: `CGCGAATTCGCGATCGATCG`
  - strand 2: `CGATCGATCGCGAATTCGCG`

## What to compare / send back
For each job, the readout is the **arrangement of the 5 subunits**:
- **Planar (P):** the 5 subunit centres lie in ~one plane (small axial spread), like native apo.
- **Helical (H):** the 5 subunits form a staircase (large axial spread, up to ~35 Å), like 7JQQ.
The **discriminating number** is the axial spread of the 5 subunit centroids. Eyeball first (disc vs
staircase); if you download the CIFs, drop them in `outputs/php_cycle/AF3/` and Claude will score them
with the same planarity/axial-spread metric used on the Boltz outputs (so the two predictors are
directly comparable).

## Honest caveat (same for both predictors)
AF3/Boltz may be biased toward one basin and **under-respond** to the ligand for a large allosteric
motion. A clear ATP→helical shift is strong evidence; a null result is inconclusive and escalates to the
enhanced-sampling MD (Aim C3 in `../../COMPUTATIONAL_AIMS_PHP_CYCLE.md`). Cross-checking AF3 vs Boltz is
exactly to guard against one predictor's bias.
