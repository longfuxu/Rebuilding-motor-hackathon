# Strongest current evidence that the cp233 single-chain gp16 ring is designable/viable

_The weight comes from **convergence of independent methods**, not any single number. All evidence below
is computational; there is no wet-lab validation yet (the honest top-line caveat)._

Construct under test: **cp233_int15_inter10** — a 1750-aa single chain = 5 covalent gp16 copies,
circular-permuted (CP cut at 233), int/inter linkers 15/10.

| # | signal | result | where |
|---|---|---|---|
| 1 | **Cross-predictor fold-in-context** | closes an in-order ring, **5/5 on AF3 + Boltz + OF3** (three independent predictors), using the tiled block-diagonal MSA — robust, unlike the marginal RFdiffusion rings | `outputs/af3_2026_07_09/AF3_INDEX.md`, `AF3_SCORING_REPORT.md`; `pipelines/tiled_msa_fold/` |
| 2 | **It is the native gp16 fold** | US-align to native **1.8 Å / TM 0.94** — the "off" look at first glance was the linker, not the fold | `outputs/` structural-credibility overlay; `reproduce/score_m2.py` |
| 3 | **Catalytic coupling geometry (M2)** | trans arginine-finger R146 → neighbour Walker-A engaged, native-like, **5/5, handedness-robust** (scored both k→k±1) | `reproduce/score_m2.py`; `reproduce/functional_contacts_7jqq.json` |
| 4 | **Predictor-independent MD (OpenMM)** | design stays **closed & symmetric like native apo**, unlike the 7JQQ helical state — a dynamical check that knows nothing about the folding predictors | `md/openmm_validation/` (RESULTS.md, METHODS.md); memory [[md-openmm-validation-result]] |
| 5 | **Interface energetics (MM-GBSA)** | trans-interface **−192 vs native −174 kcal/mol** — as/more favourable than native | `outputs/independent_validation/interface_energy/REPORT.md` |
| 6 | **DNA-channel competent (M3)** | pore **20.8 Å admits dsDNA ≈ native apo 20.6 Å**; distinguishes function beyond mere ring closure (cp285/297 close but their channels are too narrow) | `outputs/m3_m4_m5/` |
| 7 | **Biochemistry-consistent** | force-transmission residue **Y129 validated against the φ29 mutation dataset**; ESM naturalness supports the designed sequence | `outputs/excel_mutation_analysis/`; `reproduce/esm_naturalness.py` |

## What does NOT change this
- **BioEmu equilibrium sampling of the full 1750-aa ring** produced expanded conformations that differ
  from 1–7 (`outputs/bioemu/RESULTS.md`). Recorded as an open observation; BioEmu is a sequence-only
  de-novo generator (it does not start from the fold), so it is a different kind of test than 1–7
  (which start from or converge to the fold). Not resolved; see that run log + its open questions.

## The honest gaps (state these plainly)
- **All evidence is in silico.** No expression, no cryo-EM, no single-molecule assay yet.
- The MD was implicit-solvent, short (2–3 ns); explicit-solvent + longer + steered/targeted MD are the
  documented escalations (`md/STEERED_TARGETED_MD_PLAN.md`).
- "Designable" ≠ "functions as a motor": the coupling/power-stroke case is proxy-based (PRS,
  force-network, soft-mode) and awaits the driven-cycle simulation and wet-lab test.
