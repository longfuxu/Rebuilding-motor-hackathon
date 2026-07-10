# AF3 folds — consolidated index (all AlphaFold3 runs, 2026-07-09)
3rd independent predictor. M2 handedness-robust, tiled where noted. Detail: `AF3_SCORING_REPORT.md`, `../e119q_functional/REPORT.md`.

| construct | type | AF3 M2 | verdict |
|---|---|---|---|
| cp233_WT (apo) | cp233 baseline | 5/5 | closed |
| cp233_WT + ATP·Mg | cp233 baseline +ligand | 5/5 | closed; 5 ATP+5 Mg; catalytic baseline |
| cp233_E119Q_1seat (apo / +ATP·Mg) | addressable dead seat | 5/5 / 5/5 | assembly-tolerant |
| cp233_E119Q_5seat (apo / +ATP·Mg) | addressable dead seat | 5/5 / 5/5 | assembly-tolerant; Mg intact (silent knockout) |
| cp233_novel_d5 | de-novo (~53% id) | 5/5 unanimous | cleanest de-novo |
| cp233_novel_d2, d3 | de-novo | 5/5 | pass |
| cp233_novel_d7, d8 | de-novo | 2/5 | genuinely scrambled (handedness-checked) |
| cp285_int15_inter10 | 2nd CP site | 5/5 | closes (but M3 channel too narrow) |
| cp297_int15_inter10 | 3rd CP site | 5/5 | closes (but M3 narrow + cuts DNA-contact) |
| rfdiff_ring d0 / d2 / d1 | generative (Mode B) | 0/5 / 2/5 / 2/5 | does NOT close on AF3 (predictor-split vs Boltz 5/5) |

**Takeaways:** cp233 family robust across 3 predictors incl. ligand; E119Q assembly-tolerant + catalytically silent; cp285/297 close but fail M3; RFdiffusion ring predictor-split (Boltz yes, AF3 no).
