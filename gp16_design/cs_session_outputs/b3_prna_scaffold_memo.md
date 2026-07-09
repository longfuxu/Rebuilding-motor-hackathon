# B3 — pRNA-scaffold route to position-addressability (feasibility memo, no fold)

## The idea
Instead of a giant single chain, use the **pRNA scaffold that phi29 already provides** to template a defined subunit order. If each gp16 subunit is recruited by a *distinct* pRNA molecule, and the pRNAs are forced into a defined cyclic order by engineered loop–loop complementarity, then the subunit at each ring position is set by which pRNA occupies that position — no covalent protein fusion, no folding-unit severing, no 1,700-aa chain.

## Structural basis (measured from 7JQQ, heavy-atom <4.5 Å)
Each gp16 subunit contacts exactly ONE pRNA chain — a clean 1:1 position code already present in the native motor:

| gp16 seat | pRNA chain | gp16 contact residues (n) |
|---|---|---|
| A | K | 19 (e.g. 10, 13, 14, 15, 16, 17…) |
| B | L | 21 (e.g. 10, 11, 13, 14, 15, 16…) |
| C | M | 15 (e.g. 10, 11, 13, 14, 15, 16…) |
| D | N | 20 (e.g. 10, 13, 14, 15, 16, 17…) |
| E | O | 19 (e.g. 10, 13, 14, 15, 16, 17…) |

The gp16 pRNA-contacting residues cluster at the N-terminal/pRNA-binding surface (the 10–17, 37–41, 148–152, 237–239, 268–273 patches annotated in the topology map), i.e. away from the Walker/arginine-finger catalytic machinery — so pRNA-based positioning does **not** compete with the trans-R146 interface that B1 must preserve.

## Engineering route (concept; cites the phi29 pRNA literature)
1. **Interlocking-loop complementarity.** phi29 pRNA oligomerizes via right/left interlocking loops (pRNA interlocking loops; [refs supplied in project brief; not independently verified in this session]). Native loops are self-complementary, giving a symmetric ring with no positional information. Replace the wild-type loop pairs with an **orthogonal set of engineered complementary sequences** (e.g. a 5-membered set where loopᵢ-right pairs only with loopᵢ₊₁-left), so the pRNA ring can assemble in exactly ONE cyclic order.
2. **Couple a subunit identity to each pRNA.** Because each pRNA already binds one gp16 seat, appending a subunit-specific recruitment element (or simply using the WT gp16–pRNA interface) makes the ordered pRNA ring impose an ordered gp16 ring.
3. **Address one position.** Put the mutation/label on the single gp16 whose pRNA carries the unique loop pair for the target seat. This addresses one ring position without touching the other four subunits' sequence.

## Why this is a distinct design rule, not just a fallback
- It **sidesteps R1/R3 entirely**: no folding unit is severed and no inter-subunit protein interface is put under covalent tension — the interface remains fully native, assembled from five separate WT-fold subunits.
- The addressing information lives in the **RNA scaffold**, which phi29 uniquely offers among ring ATPases. This is the phi29-specific route a generic single-chain method (B2) cannot claim.
- Cost/feasibility: pRNA loop mutagenesis is standard in the Guo-lab toolkit; the readout (defined-position label) is the same optical-tweezers/fluorescence assay the project already targets.

## Standing caveat
7JQQ is the helical state; whether an engineered-pRNA ring locks the *resting* (planar) subunit order or is itself remodeled during the planar↔helical cycle is exactly what a position-labeled construct would test (planar↔helical & dynamic-special-subunit; [refs supplied in project brief; not independently verified in this session]). B3 addresses a geometric seat, not a guaranteed fixed functional role.
