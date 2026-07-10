# tiled_msa_fold results — gp16 single-chain rings (Boltz-2 NIM, tiled block-diagonal MSA)

**M2(ring)** = trans-R146 arginine-finger -> sequential-neighbour Walker-A < 8 Å; pass >=4/5.
Scored handedness-robust: a near-C5 ring can wind either way in one diffusion sample, so
engagement is counted toward the DESIGNED neighbour (k->k+1) or its MIRROR (k->k-1),
whichever is coherent. **M1** = ring geometry (compact + planar + sequential register).
Never global pTM. Depth-1167 monomer MSA, tiled to 5836 rows, ~55-70 s/fold.

## Key findings
- **All 7 constructs close (M2 5/5, compact planar ring).**
- **Dead seat does NOT open the ring:** cp233 WT, E119Q_1seat and E119Q_5seat are all
  5/5 with all R146 fingers engaged at ~5.8-6.2 Å and identical geometry (radius ~26 Å,
  CV ~0.01). The naive direction-specific M2 shows WT 5/5-forward but E119Q 0/5-forward;
  that is a *ring-handedness* artifact — the E119Q rings are equally closed (5/5 reverse).
  Confirmed by 3x replicates each (handedness_replicates.json): WT wound
  designed/designed/mirror, E119Q wound mirror/designed/designed — i.e. BOTH sequences
  wind both ways, so handedness is a per-sample coin-flip decoupled from the mutation.
  Across all 13 folds (7 constructs + 6 replicates) M2(ring) = 5/5 every time.
- **Tiled MSA is required:** the RFdiffusion native-order rings score 0/5 single-seq but
  **5/5 tiled** (d0, d1) — single-seq cannot evaluate these single-chain rings.
- Secondary CP sites **cp285 and cp297 also close 5/5** under tiled MSA.

| construct | len | M2(ring) | pass | handedness | radius(Å) | CV | compact | ifacepLDDT | wall |
|---|---|---|---|---|---|---|---|---|---|
| cp233_WT | 1750 | 5/5 | PASS | designed(k->k+1) | 25.5 | 0.00 | True | 63.0 | 60.7s |
| cp233_E119Q_1seat | 1750 | 5/5 | PASS | mirror(k->k-1) | 26.5 | 0.01 | True | 64.4 | 54.6s |
| cp233_E119Q_5seat | 1750 | 5/5 | PASS | mirror(k->k-1) | 26.5 | 0.01 | True | 64.3 | 55.5s |
| cp285_int15_inter10 | 1750 | 5/5 | PASS | designed(k->k+1) | 25.9 | 0.00 | True | 62.1 | 54.4s |
| cp297_int15_inter10 | 1750 | 5/5 | PASS | designed(k->k+1) | 24.2 | 0.06 | True | 61.9 | 55.2s |
| rfdiff_ring_L50_d0 | 1835 | 5/5 | PASS | designed(k->k+1) | 28.0 | 0.01 | True | 63.1 | 72.6s |
| rfdiff_ring_L50_d1 | 1835 | 5/5 | PASS | designed(k->k+1) | 27.5 | 0.02 | True | 65.5 | 66.0s |
