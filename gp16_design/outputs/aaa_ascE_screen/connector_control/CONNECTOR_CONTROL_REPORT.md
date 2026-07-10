# Connector-necessity control for the Rho generative ring winners

**Question.** The two cross-predictor-confirmed generative Rho rings (`rhosal_L30_40_0_d0`,
`rhosal_L40_52_2_d0`) close with an RFdiffusion-designed rigid connector whose sequence was set by
ProteinMPNN. Was that generative design *necessary*, or would any linker of the same length work?

**Design.** For each winner, the designed connector was replaced by two controls of identical length,
keeping the 6× native Rho motor ring otherwise unchanged:
- **scramble** — the designed connector's residues, composition-matched but order-shuffled (seed 37).
- **gslinker** — a plain (GS)ₙ flexible linker of the same length.

Each was folded with the same Rho tiled MSA on **Boltz-2** and **OpenFold3**, and scored with the
same handedness-robust M1/M2 and the M3 channel gate. A control "closes" only if it clears the same
bar as the winners: M2 6/6 designed-handedness on **both** predictors + M3 clean.

## Results (`connector_control_table.csv`)

| winner | connector | Boltz M2 / M3 | OF3 M2 | verdict |
|---|---|---|---|---|
| rhosal_L30_40_0 | designed (RFdiffusion+MPNN) | 6/6 / clean | 6/6 | CLOSES (reference) |
| rhosal_L30_40_0 | scramble (L32) | 0/6 / flag | 6/6 | **fail** (Boltz 0/6) |
| rhosal_L30_40_0 | **(GS)₁₆ linker (L32)** | 6/6 / clean | 6/6 | **CLOSES (Boltz+OF3, M3 clean)** |
| rhosal_L40_52_2 | designed (RFdiffusion+MPNN) | 6/6 / clean | 6/6 | CLOSES (reference) |
| rhosal_L40_52_2 | scramble (L42) | 6/6 / clean | 6/6 | CLOSES (Boltz+OF3, M3 clean) |
| rhosal_L40_52_2 | **(GS)₂₁ linker (L42)** | 6/6 / clean | 6/6 | **CLOSES (Boltz+OF3, M3 clean)** |

## Conclusion — honest negative on generative necessity

**The generative (RFdiffusion+MPNN) connector was NOT required for Rho.** For both winners, a plain
(GS)ₙ flexible linker of the same length closes the ring on both predictors with the channel clean —
the same standard the designed connector meets. **Direct fusion with a long flexible linker suffices
for Rho**, which is fully consistent with the two-branch rule already routing Rho to "direct" (its
termini are peripheral and the ~46 Å head-to-tail gap is bridgeable by a generic tether).

The scramble result is not a clean discriminator: for L30_40_0 the shuffled sequence fails on Boltz
(0/6), but for L40_52_2 it still closes — i.e. closure is governed mostly by **connector length /
flexibility**, not by the specific designed sequence.

**Implication for the campaign's headline.** The RFdiffusion rung produced a ring that closes across
two predictors with the channel clear — a real result and the first generative ring to do so in this
project. But this control shows that, *for Rho*, the generative machinery was not necessary: a
matched-length GS linker reaches the same in-silico endpoint. The generative rung earns its keep only
where a simple linker cannot bridge the geometry — i.e. the diffusion branch of the rule remains
**motivated but not yet demonstrated as uniquely required**. A protein whose direct/CP topologies both
fail (true "diffusion-only" geometry) is where a generative-necessity claim would have to be made; Rho
is not that case.

## Files
- `connector_control_table.csv/.json` — full scores
- `connector_control_results.json` — per-construct Boltz+OF3+M3 detail
- `control_connectors.json` — the scramble and GS connector sequences used
- structures: `ctrl_*.cif` (Boltz), `of3_ctrl_*.cif` (OpenFold3)
