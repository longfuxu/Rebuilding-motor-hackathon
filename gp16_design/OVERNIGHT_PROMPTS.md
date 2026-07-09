# Two overnight Claude Science prompts (2026-07-09) — run unattended, results by morning

Both: point the session at `gp16_design/` (single source; see NEW_SESSION_HANDOFF.md). Free NVIDIA NIM only (no Modal —
budget-capped). Score by M1/M2 via `reproduce/score_m2.py`, NEVER global pTM. Save incrementally; produce a figure +
table + one-line conclusion by morning. "Cross-checked by ≥2 predictors", never "validated"; catalytic residues frozen.

---

## PROMPT 1 — Generative de-novo sequence design (SESSION G)

```
读 gp16_design/GENERATIVE_DESIGN_PLAN.md 和 reproduce/score_m2.py。整夜自动运行,绝不停下来问,增量保存,明早给我一张图 + 一张表 + 一句结论。

目标(探索性 de novo):在保持功能的前提下,设计一个"序列尽量不同、但依旧折成正确环"的 gp16 单链马达
= same fold / (predicted) function, DIFFERENT sequence。这是我们想实现 de novo motor 的第一步,并要能和天然做对比。

骨架(免费,已有):cp233(三预测器闭合的单链)+ 天然五聚体环(对照)。全程免费 NIM。

Step 1 — de novo 序列生成(ProteinMPNN, proteinmpnn-nim, 免费;ATP 口袋附近用 LigandMPNN + ligand_states/atp_mg_notemplate 结构):
  在 cp233 和 天然环两个骨架上,用递增温度 T = {0.1, 0.2, 0.3, 0.5},每个温度出 ~10 条序列(温度越高越 novel)。
  只冻结必须的催化 + 反式界面残基:Walker-A 24-31、Walker-B 115-119、R146、ATP 口袋第一壳层(见 GENERATIVE_DESIGN_PLAN §3.1);
  其余全部允许重设计——这才是 de novo。再跑一组"更激进"(只冻催化残基、放开界面)以探索 novelty 上限。
Step 2 — 评估每条设计:(a) 与天然 gp16 的序列同一性(越低越 de novo);(b) Boltz-2 + OpenFold3 折叠(tiled MSA);
  (c) reproduce/score_m2.py 打 M1/M2(cp233 骨架用 CP 参数 --copy_start_res 1 --r146_incopy/--walker_incopy,值见 manifest;天然环用默认)。
Step 3 — 排序 + 产出:按 (novelty × closure) 排序,找出"与天然同一性最低、但仍 M1 sequential=YES 且 M2 ≥4/5(两预测器)"的设计
  = de novo motor 候选。画 novelty(x=与天然同一性%)vs closure(y=M2/5)散点 + Pareto 前沿;做 native cp233 vs de novo 对比表。
  把最好的 3-5 条(最 novel 且闭合)单独列出发我跑 AF3(第三预测器)。
守则:sequential M2 不看 global pTM;催化残基永远冻结;免费 NIM 不用 Modal;整夜增量保存。
```

---

## PROMPT 2 — Systematic circular-permutation site screen (SESSION CP)

```
读 gp16_design/cs_session_outputs/candidate_cp_sites.csv 和 reproduce/score_m2.py。整夜自动运行,绝不停下来问,增量保存,明早给我一张地形图 + 一张表 + 一句结论。

目标:系统性筛选环状置换切点(不只是 228-234),向审稿人证明我们做了系统 screen 而不是碰巧命中一个位点,并理解不同切点的功能特性。

切点集合 = candidate_cp_sites.csv 里 func_clear=True 的候选,跨整个蛋白,排除已跑的 220-240 细网格:
  {60, 75, 89, 124, 164, 178, 183, 205, 217, 248, 258, 280, 285, 297}
  这些从低 cut_penalty 的折叠单元边界(205/217)到高 cut_penalty 的结构域内部(280/285/297,预期会碎)——
  目的正是建立并检验设计规则:低 cut_penalty / 折叠单元边界 → 闭合;高 penalty / 域内切 → 失败(碎折叠单元)。
内部 linker 固定 15;亚基间 linker 先 15,过关的再补 10 / 20。
CP 亚基 = [res P+1..330] + (GGGGS→15) + [res 4..P];五聚体 = 5 亚基 + 亚基间 linker。
打分位置(从切点 P 算):R146 = 488−P;Walker-A = 366−P .. 373−P;--copy_start_res 1。

Stage 1 — 每个切点折 CP 单体(Boltz-2 + 置换 ColabFold MSA),叠到天然 gp16,报 N-domain(4-200)/C-domain(229-330) RMSD;
  N≤3.5Å 且 C 不碎(不 fragment)才进 Stage 2。
Stage 2 — 过关切点折 CP 五聚体(Boltz 3 seeds + OpenFold3,tiled MSA),按 sequential M2 + M1 打分;两预测器 ≥4/5 = PASS。
产出:切点 × 是否闭合 的地形表 + 一张 cut_penalty(或"是否折叠单元边界")vs 闭合 的相关图(这就是"系统性 + 规则"的证据),
  和已知的 cp232-234 甜点区放一起对比;标出任何新的闭合位点及其功能特性。
守则:免费 NIM;sequential M2 不看 global pTM;整夜增量保存。
```
