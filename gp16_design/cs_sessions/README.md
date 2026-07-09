# Claude Science session outputs — mirrored into the repo (2026-07-09)

The three CS session workspaces for the gp16 / ClpX single-chain-ring project, mirrored
here so the repo is the single source of truth (CS workspaces are per-session local and
not auto-synced). Excluded during mirror: `.venv/`, `.cache/`, `.tmp/`, `__pycache__/`,
`hpc/` (scratch), `handoff/tiled_msa/` (95M regenerable MSAs). Source base:
`~/.claude-science/orgs/34666d1b-.../workspaces/<id>/`.

## session1_main_gp16/  (id 9c944ce5…, 366M, 327 structures)
Main gp16 workspace. Cheap→expensive fold funnel; every construct judged by M1/M2/M4
under ≥2 predictors. **Leads: cp233_int15_inter10 (3-predictor-robust single chain) and
B1_L40_in_ring (covalent-dimer addressable).** Dead-seat E119Q tolerated at all patterns
(Boltz+OF3). Read first: `outputs/RUN_SUMMARY.md`, `outputs/folding_results.md`,
`handoff/cp233_3predictor_result.md`, `handoff/af3_3predictor_sweep.md`,
`handoff/af3_ligand_states.md`. Structures under `outputs/structures/**`.

## session2_denovo_clpxvalidation/  (id 95d796ed…, 8.5M)
"De novo ring motors": ProteinMPNN redesign of cp233 → **cp233_NOVEL** (10/12 pass both
predictors at ~53% identity, ~47% redesigned; D5 flagship). Plus **ClpX retrospective
validation** in `clpx_out/` (good 6/6, R307A/R307E 0/6; native 6PP5 = 5/6 calibration).
Read: `denovo_gp16_report.md`, `clpx_out/clpx_validation_report.md`. Sequences:
`af3_candidates.fasta`, `clpx_out/clpx_constructs.fasta`. Structures: `cp233_NOVEL_D5__*.cif`,
`clpx_out/clpx_*__of3.cif`.

## session3_clpx_vs_gp16_topology/  (id e1736e32…, 3.4M)
ClpX-vs-gp16 topology decision: **ClpX → direct C→N head-to-tail fusion** (termini off-pore,
~20Å gap); **gp16 → circular permutation** (native C-term fouls the DNA channel). Read:
`clpx_vs_gp16_topology.md`, `clpx_topology_map.json`, `clpx_vs_gp16_topology.png`;
source `6PP5.cif`.

## Note on versioning (2026-07-09)
`session1_main_gp16/outputs/structures/` (~341M of raw Boltz/OF3/AF3 `.cif` fold models
from the CP screen / ladder / ligand states) is **git-ignored** — preserved on local disk
at this path but not versioned, to keep the repo lean. Git holds all reports, score tables
(`*.csv`), figures, and the **definitive cp233 lead structures** in `session1_main_gp16/handoff/cp233_af3/`.
The scores in the CSVs are sufficient to reproduce/interpret; re-fold from `fold_inputs/` if a raw model is needed.
