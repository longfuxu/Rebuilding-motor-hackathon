# PIVOT — direct-linker single chain is dead; go to 2–3 group cut-splice

## What the AF3 sweep settled (details in handoff/af3_3predictor_sweep.md)
AlphaFold3 (5 models each, **sequences verified correct**) scores the single-chain (GGGGS)-linker pentamer
by designed-sequential M2:

| linker | Boltz-2 (tiled MSA) | OpenFold3 (tiled MSA) | AlphaFold3 (own MSA) |
|---|---|---|---|
| L30 | 5/5 | 5/5 | **0/5 scrambled** |
| L34 | 5/5 | 5/5 | **2/5 scrambled** |
| L20/L25/L36/L50 | — | — | all scrambled (≤2/5) |

**AlphaFold3 scrambles the direct-linker single chain at EVERY length.** The L30/L34 "fix" was 2-predictor
only. So the direct-linker single chain is NOT 3-predictor robust → **abandoned** (user decision).

Ligand-state win: native ring + 5 ATP + 5 Mg → all 5 R146 **symmetric tight (3.14–3.23 Å, pLDDT ~81)** vs apo
~6.6 Å. ATP tightens every finger; fully-loaded ring is symmetric (resolves the 7JQQ partial-occupancy caveat).

## NEW DIRECTION (user): 2–3 covalent groups (cut-splice), all positions defined
Instead of one 5-copy chain (scrambles) or B1's dimer+3 free WT, partition the pentamer into **2–3 covalent
pieces that co-assemble**. Each piece is short (dimer ~684, trimer ~1041 aa) → far less tandem-repeat
scramble, AND every subunit position is covalently defined (better addressability than B1). Precedent:
**single-chain ClpX** (covalently-linked subunits with defined-position W/M mutants; Sosa et al. 2026 / Martin-
Baker-Sauer) is proven to assemble and function — so covalent multi-group linking is legitimate, and part of
the single-chain scramble may be a *predictor artifact* on tandem repeats that short pieces mitigate.

**Constructs staged in `handoff/cut_splice/` (native-order tether, (GGGGS)×6 = the best sweep linker):**
- **Scheme A (2 groups, 3+2):** `trimer_s123` (seats 1-2-3) + `dimer_s45` (seats 4-5). 2 free seams (3→4, 5→1).
- **Scheme B (3 groups, 2+2+1):** `dimer_s12` + `dimer_s34` + `monomer_s5`. 3 free seams.

## TASK (free NIM only; sequential M2; save incrementally)
1. **Fold each scheme with ALL its chains together** (Boltz-2 ≥3 seeds + OpenFold3), tiled core MSA per chain.
   Score every subunit's R146 vs its DESIGNED neighbour (covalent within a piece; co-assembly at the seams).
   Pass = the full 5-subunit ring reforms ≥4/5 under BOTH predictors (like B1 did, unlike the single chain).
2. **Adjust the cut sites if a scheme doesn't close:** (a) sweep the internal linker length; (b) try the other
   grouping (e.g. tetramer+monomer, or move which boundary is the free seam); (c) as a last resort, a
   circular-permutation cut at the hinge (res 201–228, folding-unit boundary per rule R1) to shorten a junction.
3. Compare to B1 (the validated dimer-in-ring) and report which grouping is most robust.
4. Best scheme → hand back for AF3 third-predictor confirmation (user runs AF3).
Guardrails: cross-checked not validated; apo; sequential M2 not global pTM; report divergence.
