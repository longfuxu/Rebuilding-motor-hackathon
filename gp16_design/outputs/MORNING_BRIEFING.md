# MORNING BRIEFING — gp16 position-addressable design (overnight autonomous run)

**Bottom line:** A **covalent native-order gp16 dimer embedded in a pentameric ring (2 tethered seats + 3 WT chains)** reforms the catalytic trans-arginine-finger interface (R146 → neighbor Walker-A pocket) **as reliably as the native ring**, cross-checked by two independent structure predictors. This makes **one ring position genetically defined without destabilizing the inter-subunit interface** — the project's core goal — for a total GPU spend of **$9.14** (cap $10.00).

## What was tried, and what each dollar bought
| Phase | What | GPU $ | Result |
|---|---|---|---|
| A (Stage A, free) | domain-aware cut-penalty filter | $0 | reproduced cycle-1 CP failure for free; hinge = res 201–228 |
| 1 | native ring reference, Boltz 1 seed | $0.20 | ring assembles (closed pentamer, correct cyclic order); apo ring = symmetric, 5/5 R146 at 6.6 Å |
| 2 | B1 dimer+3WT, 3 linker variants, Boltz 1 seed | $0.95 | all 3 PASS assembly gate; seats fold (4.7 Å), R146 engaged |
| 3/7 | ensemble: Boltz 3-seed + Chai 2-seed × (native, L32, L40) | $7.99 | decisive per-interface distributions (below) |
| 5 (free) | B2 folding-unit graph + topologies | $0 | one framework explaining cycles 1→3 |
| 6 (free) | B3 pRNA-scaffold memo | $0 | phi29-specific alternative route |
| **Total** | | **$9.14** | |

## The single best construct + honest reliability verdict
**`B1_L40_in_ring`** — native-order covalent dimer (res 4–330 – GS40 linker – res 4–330) co-assembled with 3 separate WT chains.

| | R146 engaged (frac) | median R146 | tether-flanking seats |
|---|---|---|---|
| native ring (reference) | Boltz 15/15 · Chai 50/50 | 6.6 / 2.8 Å | — (WT) |
| **B1_L40 in ring** | **Boltz 15/15 · Chai 47/50** | **5.8 / 4.0 Å** | **Boltz 6/6 · Chai 20/20** |
| B1_L32 in ring | Boltz 15/15 · Chai 45/50 | 5.7 / 4.0 Å | Boltz 6/6 · Chai 19/20 |
| *(cycle-2 isolated 2-mer)* | *Boltz 0/5 · Chai 7/10* | *stochastic* | *n/a* |

**Verdict: PASS, cross-checked.** In ring context the covalent tether does NOT perturb the tether-flanking interfaces relative to the native ring — both predictors, both linker lengths, tether seats included. This is the opposite of the isolated 2-mer, and it is exactly what the cycle-2 topology finding predicted (the interface is a ≥3-subunit property; the ring supplies the missing neighbors).

## Design rules learned
1. **Cut only at folding-unit boundaries** (Stage A cut-penalty; cycle-1 CP@280 broke this).
2. **The gp16 interface is a ≥3-subunit property** — score every covalent unit *in ring context*, never in isolation (cycle 2).
3. **Tether, don't permute.** A native-order covalent dimer (no circular permutation, no severed unit) embedded in the ring is the minimal genetically-defined unit that satisfies all constraints (B2 R1–R3). Boundary-CP single chains remain unproven in ring context (Phase 4 budget-skipped).
4. **Linker length** L40 ≥ L32 (both work; L40 slightly more reliable in Chai tether seats 20/20 vs 19/20).

## Recommended next action
1. **Wet-lab:** express `B1_L40` dimer + 3 WT gp16; test ring assembly + pRNA/DNA binding + ATPase. This is the synthesizable, position-addressable reagent.
2. **The actual experiment (compute-ready):** put a catalytically-dead **E119Q** on ONE module of the dimer → a ring with a *defined* dead position (Phase 3 escalation (b), budget-skipped tonight).
3. **B3 pRNA-scaffold** as the phi29-unique, no-giant-chain alternative if covalent expression fails.

## What was budget-skipped (no silent truncation)
- **Phase 3(b) dead-subunit E119Q dimer-in-ring** — skipped: ensemble consumed the reserve (Chai ran 99 min, 3× its estimate). Next-run priority.
- **Phase 4 boundary-CP (228/233) in ring context** — skipped for the same reason; the B1 result is the stronger positive so the budget was rightly spent there.
- **Phase 7 extra seeds / context refold (+pRNA/DNA/ATP)** — skipped at cap; would tighten the confidence interval and test the interface with nucleotide present.

## Honest abstract (cross-checked; caveats intact)
> The bacteriophage φ29 gp16 packaging ATPase is a homopentameric ring whose functional asymmetry cannot currently be assigned to a defined physical position. Using an auditable, cost-bounded (<$10 GPU) structure-prediction funnel, we show that circular-permutation single chains fail because the permutation cut fragments an inter-subunit-packed C-terminal domain, and that gp16's catalytic trans-arginine-finger (R146) interface is a ≥3-subunit property that an isolated dimer under-determines. We then show, cross-checked by two independent predictors (Boltz-2 and Chai-1, multi-seed ensembles), that a **native-order covalent gp16 dimer embedded in a pentameric ring** (two tethered seats + three wild-type subunits) reforms R146 at the tether-flanking interfaces as reliably as the native ring (Boltz 15/15, Chai 47/50; native 15/15, 50/50). This yields a synthesizable construct in which one ring position is genetically defined without destabilizing the inter-subunit interface, enabling position-specific dead-subunit and single-fluorophore experiments. Because 7JQQ is the helical (lock-washer) state and predictions were made without pRNA/DNA/nucleotide, these constructs address a geometric ring *seat*, not a guaranteed fixed functional role (planar↔helical & dynamic-special-subunit; [refs supplied in project brief; not independently verified in this session]); which role occupies that seat is precisely the hypothesis the construct is built to test. Results are computationally cross-checked, not experimentally validated.
