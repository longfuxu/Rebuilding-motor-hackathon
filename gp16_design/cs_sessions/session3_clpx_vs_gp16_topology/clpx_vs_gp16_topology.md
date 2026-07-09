# Same framework, different solution: optimal single-chain topology for ClpX vs gp16

## The question
CS-B validated that the gp16 *fold-in-context + M1/M2* framework transfers to a second
ring ATPase, hexameric ClpX (good construct 6/6 engaged; R307A/R307E coupler lesions
0/6; calibrated to 5/6 on the native cryo-EM hexamer 6PP5). This memo answers the
follow-on: **what is ClpX's optimal single-chain topology — direct fusion or circular
permutation (gp16's solution)?** — decided on measured geometry of the functional
hexamer, not assertion.

## What was measured (6PP5, substrate-bound ClpXΔN hexamer, cryo-EM 3.98 Å)
The pore axis is defined by the substrate peptide (chain S) threaded through the ring;
ring order around it is C-D-E-F-A-B, with subunit F (ADP) the spiral seam.

| geometric test | ClpX (6PP5) | gp16 (7JQQ, prior) |
|---|---|---|
| native C-terminus → functional channel | **62 Å (off the pore)** | **6.2 Å — a DNA-channel contact** |
| native N-terminus → functional channel | 47 Å (off the pore) | 36 Å (DNA) / 16 Å (pRNA) |
| direct C(i)→N(i+1) head-to-tail gap | **~20 Å (19–23)** | 53–66 Å |
| circular-permutation alternative | internal join 54.7 Å; short cuts fall *inside* the pore | cut297/280 → 18–24 Å, off the DNA channel |

## The determination
**ClpX's optimal single-chain topology is a direct C→N (head-to-tail) fusion.** Three
independent lines of measured/known evidence agree:
1. **Termini are peripheral.** ClpX's native C-terminus (res414) sits 62–65 Å off the
   substrate translocation pore and the N-terminus (res62) 47–49 Å off it — a direct
   genetic fusion never approaches the functional channel.
2. **The gap is short.** The head-to-tail C→N distance around the ring is only ~20 Å,
   comfortably bridged by a genetically encoded flexible linker.
3. **Circular permutation is strictly worse here.** A CP would have to join the native
   termini across 54.7 Å internally, and the only CP cut sites that give a short
   new-terminus gap (~11–14 Å) sit 5–11 Å from the pore axis — i.e. inside the
   substrate channel. CP buys nothing because ClpX's termini are already off-channel.
4. **Precedent.** Covalently tethered single-chain ClpXΔN pseudohexamers (Martin/Baker/
   Sauer lineage), built by exactly this direct C→N genetic linkage, degrade substrate
   at near-WT rates.

## Same framework, different solution
Both motors are validated by the *same* fold-in-context M1/M2 coupler metric — trans
arginine finger → neighbor Walker-A, engaged < 8 Å, cyclic, with an oligomer-appropriate
closure floor (ClpX ≥5/6, gp16 ≥4/5). But the **optimal topology diverges because the
geometry of the native termini diverges**:

- **gp16 must circular-permute.** Its native C-terminus (res330) is a dsDNA-channel
  contact; a naive C→N fusion would foul the translocation channel, and the head-to-tail
  gap is 53–66 Å. Relocating the join to a surface loop (cut ~297) drops the gap to
  ~18 Å and moves it off the DNA channel.
- **ClpX need not.** Its native termini already sit on the outer rim, off the pore, and
  only ~20 Å apart — direct head-to-tail fusion is both sufficient and optimal.

**The framework is one; the topological solution is dictated by where each motor's
termini sit relative to its functional channel.** That is the "same framework, different
solution" result: a general design principle (size the inter-subunit linker to the
functional geometry) instantiated two different ways by two different rings.

## Provenance
- ClpX geometry measured this session from RCSB 6PP5 (all distances min heavy-atom or
  CA, pore axis = PCA of substrate chain S).
- ClpX catalytic residues (R307 trans arginine finger, K125, E185, R370) grounded in
  CS-B from eLife 2020;52774 + UniProt P0A6H1; not guessed.
- gp16 values from the project's topology_map.json (7JQQ, UniProt P11014).
- Geometric/structural analysis only; no wet-lab validation.

## Artifacts
- clpx_topology_map.json — grounded ClpX topology map (schema matches gp16's topology_map.json)
- clpx_vs_gp16_topology.csv — 19-row side-by-side comparison
- clpx_vs_gp16_topology.png — 2-panel comparison figure
