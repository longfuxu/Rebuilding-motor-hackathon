# OVERNIGHT — systematic circular-permutation closure screen (autonomous, free NIM)

Run unattended overnight. No decisions needed — build, fold, score, save incrementally, produce a landscape map.
Context: cp233_int15_inter10 already closes the ring under both Boltz (3/3 seeds 5/5) and OpenFold3 (5/5).
Goal now: **map the CP closure landscape** — which cut site × linker length closes cross-predictor — to confirm
cp233 is optimal (or find better) and produce a systematic figure. All free NVIDIA NIM; sequential M2, not global pTM.

## Constructs to build (same recipe as the earlier CP task; manifest positions apply)
CP subunit for cut site P = `[gp16 res P+1..330] + internal-linker(GGGGS→15) + [gp16 res 4..P]`. Pentamer = 5 CP
subunits joined by inter-subunit linker (GGGGS→L). R146 within a copy = len(part1)+15+143; Walker-A = +len(part1)+15+21..28.
- **Cut sites (grid across the hinge/C-domain boundary):** 220, 225, 228, 230, 232, 233, 234, 236, 238, 240.
- **Inter-subunit linker:** 10, 15, 20.

## Stage 1 — monomer gate (fold 10 CP monomers, one per site, internal-linker 15)
Fold each CP monomer (Boltz-2, permuted ColabFold MSA on the CP subunit). Superpose onto native gp16; report
N-domain (4–200) + C-domain (229–330) RMSD, R146/Walker intact, pLDDT. Advance a site only if it folds native-like
(N ≤ ~3.5 Å, C not fragmented), as cp233 did.

## Stage 2 — pentamer closure (for advancing sites × inter-linker 10/15/20)
Fold Boltz-2 (≥3 seeds) + OpenFold3, tiled CP MSA. Score **sequential M2** (each copy's R146 → its DESIGNED next
copy, <8 Å; use the CP within-copy positions above) + **M1** (ring closure). Pass = full ring ≥4/5 under BOTH
predictors (like cp233). Save every fold incrementally.

## Deliverable
A **site × inter-linker closure map** (pass/fail + Boltz seeds-close + OF3), the best site(s), and a comparison to
cp233. If a site beats cp233 (more Boltz seeds robust at more linker lengths), flag it. Write to RUN_SUMMARY +
folding_results + a figure. If you finish early: deepen the top 1–2 constructs with more Boltz seeds.

## Do NOT
- Don't run E119Q yet (waits for the user's AF3 result on cp233). Don't touch the ligand-state runs (done, all closed).
- Guardrails: cross-checked not validated; sequential M2 not global pTM; report predictor divergence; apo.
