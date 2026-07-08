# B2 — Folding-unit graph: one framework for cycles 1→3

## The graph
Nodes = the 10 folding units of the gp16 pentamer: {A,B,C,D,E} × {N-domain 4–200, C-domain 201–330}. Edges = inter-unit heavy-atom contacts in 7JQQ (residue pairs <5 Å). The dominant edges are **inter-subunit**:

- N-domain ring: B_N–C_N 36, C_N–D_N 35, A_N–B_N 30, A_N–E_N 27, D_N–E_N 15
- C-domain ring: A_C–B_C 25, D_C–E_C 24, B_C–C_C 20, A_C–E_C 20
- cross contacts: A_C–B_N 32, B_C–C_N 30, C_C–D_N 20, D_C–E_N 19

Intra-subunit N–C contact is only ~3 residue pairs. **A subunit's shape in the ring is defined by its neighbors, not by itself.**

## A single-chain design = a walk over ring seats
Ring seats are angularly ordered E→D→C→B→A (and A–E closes the ring). A genetically ordered construct is a covalent walk that visits seats in some order. It is *valid* iff:

1. **R1 — never sever a folding unit.** Cuts (for circular permutation) fall only at the N/C hinge (res ≈201–228) or at chain termini, where the domain-aware cut-penalty is low. This is the rule cycle-1's cut@280 (penalty 66, mid-C-domain) broke, and Stage A now enforces it for $0.
2. **R2 — join only spatially-adjacent seats.** A covalent linker connects seat_i C-terminus → seat_{i+1} N-terminus only when the seats are ring-neighbors. Adjacent C→N spans are 54–56 Å around the ring except the **A–E seam at 66 Å** — the ATP-block/apo-block boundary, the "special/open" interface. R2 says: make A–E the *non-covalent* break, and never route a linker across it.
3. **R3 — the interface is a ≥3-subunit property.** Because both the N-ring and C-ring contacts are inter-subunit, any covalent unit must be scored *embedded in the full pentamer*, where both of its neighbors are present. This is why cycle-2's monomer floated (9.9 Å) and the 2-mer's interface was stochastic — they under-determine an interface that physically needs both flanking subunits. Cycle 3 folds every construct in ring context precisely to satisfy R3.

**Score(design) = total linker span + broken inter-subunit contacts** (minimize both).

## What the framework explains
- **Cycle 1** (CP@280 failed): violated R1 — cut inside the C-domain folding unit.
- **Cycle 2** (2-mer stochastic): violated R3 — tested an inter-subunit interface with only one neighbor present.
- **B1** (dimer + 3 WT, cycle 3): the minimal design that satisfies all three — it severs no unit (native termini), joins two adjacent seats, and is scored in ring context. It is the smallest genetically-defined unit that R1–R3 permit.

## Top valid topologies (native-order fusions, min total C→N span, never crossing the A–E seam)
- **E-D** (2 seats): total C→N span 54.0 Å
- **D-C** (2 seats): total C→N span 54.5 Å
- **E-D-C** (3 seats): total C→N span 108.5 Å
- **D-C-B** (3 seats): total C→N span 110.4 Å
- **E-D-C-B** (4 seats): total C→N span 164.5 Å
- **D-C-B-A** (4 seats): total C→N span 165.4 Å
- **E-D-C-B-A** (5 seats): total C→N span 219.5 Å

The 2-seat winner is any adjacent pair off the A–E seam (E-D, D-C, B-A ≈ 54–56 Å) — this is exactly the B1 dimer. Larger genetically-defined fractions (trimer E-D-C, tetramer E-D-C-B) extend the same walk while leaving the 66 Å A–E seam as the free break, so the ring can still open there during the packaging cycle.

*Caveat (standing): 7JQQ is the helical state; the resting ring may be planar and the special-subunit role may be dynamically reassigned (planar↔helical & dynamic-special-subunit; [refs supplied in project brief; not independently verified in this session]). These topologies define a geometric seat, not a fixed functional role.*
