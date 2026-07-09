# TASK — circular-permutation single-chain ring (staged, one-shot)

## Why (context)
Direct N-to-C linker single chain is dead (AF3 scrambles at every linker length). Storyline escalation:
N-C linker (failed) → **circular permutation (this task)** → diffusion (if CP fails). CP moves each subunit's
termini to a folding-unit boundary (the hinge), so the inter-subunit junction is short (~39 Å vs native 54–66 Å)
and can't thread wrong. Test in RING context (CP@228 fragmented as a 2-mer in cycle 2, but 2-mers under-
determine the interface — like B1's dimer did; must judge in the full ring).

## Constructs (in `handoff/circular_permutation/`, manifest.json has all positions)
Each CP subunit = `[res P+1..330] + INTERNAL-linker(joins native 330→4) + [res 4..P]`. R146 + Walker-A sit in
the res4..P half (positions in manifest). Sites P = **228, 233, 217** (folding-unit boundaries, rule R1).
- **monomers/** (9): cp{228,233,217}_int{10,15,20}_monomer — 3 sites × internal-linker 10/15/20.
- **pentamers/** (9): cp{228,233,217}_int15_inter{10,15,20}_pentamer — 3 sites × inter-subunit-linker 10/15/20 (internal fixed 15).

## STAGE 1 — fold the 9 monomers first (gate)
Fold each CP monomer (Boltz-2; +OF3 optional). **Question: does the circularly-permuted subunit recover a
native-like gp16 fold?** Metric: superpose the predicted CP monomer onto native gp16 (res 4–330); report
N-domain (4–200) and C-domain (229–330) RMSD, whether R146/Walker-A are intact, and pLDDT. **Pass** = domains
fold in native arrangement (not fragmented, like the native monomer). Per site, note the internal-linker length
that folds best. MSA: run ColabFold on the CP subunit sequence directly (342 aa) — one MSA per site.

## STAGE 2 — for sites whose monomer passed, fold the pentamers (sweep inter-subunit linker)
Build/fold cp{site}_int15_inter{10,15,20}_pentamer (Boltz-2 ≥3 seeds + OpenFold3). MSA: tile the CP-subunit
MSA across the 5 copies. **Score by M2 (designed-sequential) + M1 (closure/sequential) in ring context**, using
the CP-aware scorer:
```
python reproduce/score_m2.py <cif> --copies <copy_starts from manifest, each span = start..start+341> \
   --copy_start_res 1 --r146_incopy <R146_in_copy> --walker_incopy <lo-hi from manifest>
# cp228 example: --copies A:1-342,A:358-699,A:715-1056,A:1072-1413,A:1429-1770 --r146_incopy 260 --walker_incopy 138-145
```
**Pass** = full 5-subunit ring ≥ 4/5 R146 to designed neighbours under BOTH predictors (like B1, unlike the
direct-linker single chain).

## REPORT
Which (site, internal-linker, inter-linker) gives (a) a native-like monomer AND (b) a cross-predictor-closed
pentamer. Compare to B1. Best → hand back for AF3 (user runs). If NO CP construct closes → the storyline's
next rung is diffusion (RFdiffusion). Guardrails: free NIM only; sequential M2 not global pTM; report
divergence; save incrementally.
