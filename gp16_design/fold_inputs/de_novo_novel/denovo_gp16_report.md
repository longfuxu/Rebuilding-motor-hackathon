# De novo gp16 ring motors — same fold, different sequence (overnight campaign)

**Goal.** Design gp16 single-chain ring motors that keep the fold and catalytic
machinery but differ maximally in sequence from native gp16, and confirm ring
closure with two independent structure predictors.

**Method (all free NVIDIA NIM, no Modal).**
- Two backbones: **cp233** (circular-permutation single chain, closes under all 3
  predictors) and the **native pentamer ring** (non-covalent, control).
- ProteinMPNN-NIM inverse folding, temperatures T={0.1,0.2,0.3,0.5}, soluble model, Cys omitted.
- Two freeze regimes: **CONS** = freeze catalytic + trans-interface (135 native res);
  **NOVEL** = freeze catalytic only (31 res), interface open (novelty ceiling).
- Catalytic freeze = Walker-A 24–31, Walker-B 115–119, R146, and the full ATP/Mg
  first shell (28 residues, ≤5 Å of AGS/Mg in 7JQQ) — frozen in **both** regimes.
- 12 designs per pool folded by **Boltz-2** and **OpenFold3** with tiled block-diagonal
  MSA (1167-homolog core MSA placed into each of 5 ring copies), scored with
  reproduce/score_m2.py (M1 ring geometry + sequential order; M2 trans-R146↔Walker-A
  engagement; sequential M2, global pTM not used).
- PASS gate = **both** predictors M2 ≥ 4/5 **and** M1 sequential = YES.

**Result (48 designs, 96 folds).**

| pool | median id% | Boltz closes | OF3 closes | PASS both |
|---|---|---|---|---|
| cp233 de novo (NOVEL) | 54 | 10/12 | **12/12** | **10/12** |
| cp233 conservative (CONS) | 73 | 9/12 | 5/12 | 4/12 |
| native-ring de novo | 58 | 4/12 | 3/12 | 1/12 |
| native-ring conservative | 73 | 3/12 | 3/12 | 0/12 |

**15 designs pass the two-predictor de novo gate; 10 are cp233_NOVEL at ~53–54%
identity (46–47% of the sequence redesigned), closing 5/5 on both predictors with
all catalytic residues intact.**

**Conclusion.** On the covalent circular-permutation (cp233) scaffold, ProteinMPNN
alone produces gp16 motors that fold into the correct sequential ring under two
independent predictors while differing from native gp16 in ~half of all core
residues — same fold/function, different sequence. The non-covalent native ring does
not tolerate de novo redesign (0–1/12 pass), so the covalent chain topology, not the
native sequence, is what makes the ring predictor-robust.

**Top 5 candidates for AF3 (single-chain 1750-aa cp233 construct):**
- cp233_NOVEL_D5: identity 53.2% (novelty 47%), Boltz M2 5/5, OF3 M2 5/5
- cp233_NOVEL_D2: identity 53.4% (novelty 47%), Boltz M2 5/5, OF3 M2 5/5
- cp233_NOVEL_D3: identity 53.5% (novelty 46%), Boltz M2 5/5, OF3 M2 5/5
- cp233_NOVEL_D8: identity 53.8% (novelty 46%), Boltz M2 5/5, OF3 M2 5/5
- cp233_NOVEL_D7: identity 53.9% (novelty 46%), Boltz M2 5/5, OF3 M2 5/5

*Cross-checked by two predictors, not wet-lab validated; apo (no ATP/pRNA/DNA);
geometric closure only. Catalytic residues verified native in all candidates.*
