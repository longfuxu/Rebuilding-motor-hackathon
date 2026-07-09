# Rho de-novo connector backbones (RFdiffusion, 2026-07-09)

Generative method (3rd of the topology ladder) applied to Rho helicase: fix two adjacent Rho
motor subunits (res 175-414 each), diffuse a de-novo connector across the ~46 Å inter-subunit gap.
Colab (rhorfd session), sokrypton RFdiffusion, T=50. `run_inference rc=0`, 10 backbones:
- `rhosal_L30_40_{0..4}.pdb` — 30-40 residue connector
- `rhosal_L40_52_{0..4}.pdb` — 40-52 residue connector
`.trb` = RFdiffusion metadata (mapping/contig). `rhosal.progress.json` = per-job contigs.

## Next (feeds the tiled_msa_fold pipeline)
1. ProteinMPNN sequence design on each backbone (motif fixed, connector designed).
2. Tile into the Rho ring, fold with **Rho-specific tiled MSA** (single-seq is confounded — see PROJECT_STATUS §4).
3. score_m2 (Rho M2: R366 -> neighbour Walker-A K184). Compare vs Rho direct/CP (the 3-method table).
