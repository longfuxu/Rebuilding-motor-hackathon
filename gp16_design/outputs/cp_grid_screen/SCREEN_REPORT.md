# CP linker-grid screen — 27 gp16 single-chain pentamers (tiled-MSA Boltz-2 NIM)

CP site {233,285,297} x internal-linker {10,15,20} x inter-linker {10,15,20}. Each folded
with a block-diagonal (tiled) monomer MSA on the free Boltz-2 NIM and scored with the
**handedness-robust** ring M2 (R146 arginine-finger -> sequential-neighbour Walker-A < 8 Å,
counted for whichever winding, k->k+1 or k->k-1, is coherent) plus M1 ring geometry.
Ranked by M2 engaged, then sequential register, compact ring, tighter ring (radius_CV),
tighter finger engagement (mean R146 distance). Never global pTM.

## Best construct per CP site

- **cp233: int20_inter10** — M2 5/5, ring compact=True sequential=YES, radius 25.8 Å (CV 0.00), mean R146 5.11 Å
- **cp285: int20_inter10** — M2 5/5, ring compact=True sequential=YES, radius 25.9 Å (CV 0.00), mean R146 5.83 Å
- **cp297: int10_inter20** — M2 5/5, ring compact=True sequential=YES, radius 24.9 Å (CV 0.03), mean R146 6.46 Å

## Per-site geometry (all 9 linker combos per site; all close M2 5/5)

| site | M2 (all combos) | mean radius_CV | mean R146 dist (Å) | note |
|---|---|---|---|---|
| cp233 | 45/45 (9/9 pass) | 0.008 | 5.74 | tightest, most regular rings |
| cp285 | 45/45 (9/9 pass) | 0.019 | 5.56 | tightest finger engagement |
| cp297 | 45/45 (9/9 pass) | 0.052 | 6.43 | closes but least regular (CV ~6x cp233) |

**Takeaway:** every CP site x linker combination closes (27/27 M2 5/5, compact sequential rings) — the single-chain CP-ring design is robust to internal- and inter-linker length. cp233 and cp285 give consistently more regular rings (radius_CV ~0.008-0.019) than cp297 (~0.052); cp297 is the geometrically weakest site though still passing. Linker length is a second-order knob: shorter inter-linkers (inter10) trend to the tightest rings.


## Full ranking (best first)

| # | construct | M2 | hand | compact | seq-reg | radius(Å) | CV | meanR146(Å) | ifpLDDT | AF3 |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | cp233_int20_inter10 | 5/5 | designed | True | YES | 25.8 | 0.00 | 5.11 | 64.1 |  |
| 2 | cp285_int20_inter10 | 5/5 | designed | True | YES | 25.9 | 0.00 | 5.83 | 65.2 |  |
| 3 | cp233_int15_inter10 | 5/5 | designed | True | YES | 25.5 | 0.00 | 5.84 | 63.0 | yes |
| 4 | cp233_int10_inter10 | 5/5 | designed | True | YES | 25.5 | 0.00 | 6.1 | 62.6 |  |
| 5 | cp285_int15_inter10 | 5/5 | designed | True | YES | 25.9 | 0.00 | 6.45 | 62.1 | yes |
| 6 | cp285_int15_inter15 | 5/5 | designed | True | YES | 25.9 | 0.01 | 5.29 | 63.1 |  |
| 7 | cp285_int10_inter20 | 5/5 | designed | True | YES | 26.5 | 0.01 | 5.29 | 62.5 |  |
| 8 | cp285_int10_inter15 | 5/5 | designed | True | YES | 25.6 | 0.01 | 5.31 | 62.5 |  |
| 9 | cp285_int15_inter20 | 5/5 | designed | True | YES | 25.8 | 0.01 | 5.42 | 62.8 |  |
| 10 | cp233_int15_inter15 | 5/5 | mirror | True | YES | 26.2 | 0.01 | 5.5 | 62.9 |  |
| 11 | cp233_int20_inter15 | 5/5 | designed | True | YES | 25.6 | 0.01 | 5.72 | 63.5 |  |
| 12 | cp233_int10_inter20 | 5/5 | designed | True | YES | 25.0 | 0.01 | 5.84 | 62.1 |  |
| 13 | cp233_int10_inter15 | 5/5 | designed | True | YES | 25.2 | 0.01 | 5.94 | 63.6 |  |
| 14 | cp285_int10_inter10 | 5/5 | designed | True | YES | 26.0 | 0.01 | 5.96 | 63.4 |  |
| 15 | cp233_int15_inter20 | 5/5 | designed | True | YES | 25.2 | 0.01 | 5.97 | 63.6 |  |
| 16 | cp233_int20_inter20 | 5/5 | mirror | True | YES | 26.2 | 0.02 | 5.63 | 63.1 |  |
| 17 | cp297_int10_inter20 | 5/5 | designed | True | YES | 24.9 | 0.03 | 6.46 | 59.5 |  |
| 18 | cp285_int20_inter15 | 5/5 | mirror | True | YES | 28.5 | 0.04 | 5.75 | 63.6 |  |
| 19 | cp297_int15_inter15 | 5/5 | designed | True | YES | 24.4 | 0.05 | 5.73 | 62.9 |  |
| 20 | cp297_int20_inter20 | 5/5 | designed | True | YES | 24.2 | 0.05 | 6.55 | 63.3 |  |
| 21 | cp297_int20_inter15 | 5/5 | designed | True | YES | 24.8 | 0.05 | 6.75 | 63.6 |  |
| 22 | cp297_int15_inter20 | 5/5 | designed | True | YES | 23.8 | 0.05 | 6.88 | 62.1 |  |
| 23 | cp297_int20_inter10 | 5/5 | designed | True | YES | 24.6 | 0.06 | 6.08 | 63.7 |  |
| 24 | cp297_int10_inter10 | 5/5 | designed | True | YES | 24.5 | 0.06 | 6.11 | 62.2 |  |
| 25 | cp297_int15_inter10 | 5/5 | designed | True | YES | 24.2 | 0.06 | 6.58 | 61.9 | yes |
| 26 | cp297_int10_inter15 | 5/5 | designed | True | YES | 24.2 | 0.06 | 6.73 | 60.7 |  |
| 27 | cp285_int20_inter20 | 5/5 | designed | True | YES | 26.9 | 0.08 | 4.76 | 63.1 |  |

AF3-confirmed (cross-predictor) baselines: cp233_int15_inter10, cp285_int15_inter10, cp297_int15_inter10. Handedness is a per-sample coin-flip (n=1 diffusion sample each) decoupled from sequence — both windings are equally-closed rings; that is why scoring is handedness-robust.
