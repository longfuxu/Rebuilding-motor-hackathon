# RFdiffusion de-novo connector — run on Colab Pro (robust route)

Goal (GENERATIVE_DESIGN_PLAN §2): design a rigid, directional de-novo CONNECTOR that links two adjacent gp16
subunits into one chain and LOCKS the sequential register — keeping the catalytic/interface machinery native.
This is the "de novo 连接件" step. It needs a GPU; run it in the standard RFdiffusion Colab (Colab Pro), then feed
the outputs into the FREE ProteinMPNN→fold→score loop.

## Inputs (provided here)
- Motif PDB: `motif_AB_adjacent_subunits.pdb` — two ADJACENT gp16 subunits (chains A, B, res 1-332 each) carved from
  the closed native ring, so their relative geometry IS the trans interface. RFdiffusion holds these fixed and grows
  a connector; the resulting chain must keep R146(A)→Walker-A(B) engaged (the M2 interface).

## Notebook
Open the ColabDesign RFdiffusion notebook (GPU): 
`https://colab.research.google.com/github/sokrypton/ColabDesign/blob/main/rf/examples/diffusion.ipynb`
Upload `motif_AB_adjacent_subunits.pdb`.

## Two design modes (do Mode A first — simpler, lower risk)
- **Mode A — partial diffusion to RIGIDIFY (easier, lower risk).** Instead of the 2-subunit motif, upload the folded
  **cp233** single-chain model and run PARTIAL diffusion (`partial_T` ≈ 8–20 of 50) to idealize/stiffen only the
  hinge (201–228) and the permutation seam, holding Walker-A 24-31, Walker-B 115-119, R146, ATP pocket, and the
  C-domain trans strand FIXED. Contig = the full chain fixed except the noised seam/hinge window. 100–200 designs.
- **Mode B — de-novo CONNECTOR between the two subunits (what you asked).** Use `motif_AB_adjacent_subunits.pdb`;
  fix chain A and chain B as motif; diffuse a NEW connector (~10–30 residues) bridging A's C-terminal region to B's
  N-terminal region, routed AWAY from the DNA channel. In the notebook's `contigs`, express as: [chain A fixed]
  + [free/diffused segment 10-30] + [chain B fixed] (set the exact contig using the notebook's own examples — I'll
  help finalize it once you're in the notebook and can see its contig UI). 100–200 designs, a few connector lengths.

## Parameters (both modes)
- `iterations` default; `num_designs` 100–200; keep several seeds. Save backbones (PDB).
- File-size note: RFdiffusion weights are multi-GB — the notebook downloads them into the runtime automatically; don't
  transfer them from local.

## Downstream (FREE, no GPU — the same loop as everything else)
1. ProteinMPNN (proteinmpnn-nim, free) on each RFdiffusion backbone; **freeze Walker-A/B, R146, ATP pocket**;
   LigandMPNN for positions near the ATP/Mg pocket. 8 seqs/backbone.
2. Fold each with Boltz-2 + OpenFold3 (free NIM, tiled MSA).
3. Score with `reproduce/score_m2.py` (M1 ring closure + M2 designed-sequential R146; NOT global pTM).
4. Best designs → AlphaFold3 (manual, 3rd predictor). A design that closes ≥4/5 under both predictors with a rigid
   register and native catalytic residues = a de-novo-connector addressable ring (a generative comparator to cp233).

## Honest note
RFdiffusion is NOT on NVIDIA NIM and cannot run unattended in Claude Science; it needs your GPU (Colab Pro suffices).
The `google-colab-cli` route (driving Colab from Claude Code) is not set up here and I can't verify its reliability
for a multi-GB GPU tool — the browser-notebook route above is the dependable path. Everything downstream is free.
