# gp16 单链设计：三方法对比 + 下一步任务清单（2026-07-09）

作者视角：Claude Code。所有判定用 sequential M2（reproduce/score_m2.py）+ M1，**不看 global pTM**；apo；交叉预测器 = filter 不是 validation。

---

## Task 1 —— 三大类单链构建方法对比：谁最好？哪些还没 AF3？

三种“把五个 gp16 亚基做成一条链”的拓扑生成方法，各自最优结构 + 三预测器判定：

| 方法 | 最优结构 | Boltz-2 | OpenFold3 | AlphaFold3 | 判定 |
|---|---|---|---|---|---|
| **① N-C 首尾直连**（passive linker） | pentamer **L30**（(GGGGS)×6，1755 aa） | 5/5 | 5/5 | **0/5 · SCRAMBLED** | ❌ AF3 否决 |
| **② Circular Permutation**（切点重连） | **cp233_int15_inter10**（1750 aa） | 5/5（3 seeds） | 5/5 | **5/5 · sequential YES** | ✅ **胜者** |
| **③ RFdiffusion de-novo 连接件** | **sal_L50** 铺环（1835 aa） | 0/5，但环紧凑+register对（single-seq） | — | **未跑** | ⚠ 拓扑对、催化界面未 engage；待 AF3 |

### 结论：**最好的单链是 Circular Permutation（cp233_int15_inter10）**，三预测器唯一全 5/5。

**为什么 CP 赢、N-C 输、RFdiffusion 卡在中间——一条清晰的机制线：**
- **① N-C 直连输在 AF3**：L30 在 Boltz+OF3 都 5/5，看起来成了，但 AF3（独立 MSA）把链穿成了**紧凑但错序**的环（radius 28.1、平面、radius_CV 0.00，但 designed 邻居相距 34–55 Å，R146 埋向**错误**伙伴 ~7 Å）。两预测器的“闭合”是假象，AF3 暴露了 scramble。**这条路死了**——也正说明为什么必须三预测器 + interface-resolved M2。
- **② CP 赢**：cp233 把每个亚基的末端搬到 C-domain 折叠单元边界 → 亚基间接头很短 → 链**没法穿错序**。三预测器全 5/5、register 全对。**definitive 单链 lead**。第二个独立 C-domain 位点 cp285 也过关（系统性、非碰巧）。
- **③ RFdiffusion 卡中间**：de-novo 连接件铺 5 份 → 折成**紧凑环 + register 全对（M1 sequential YES）**，这一点**比 N-C 强**（N-C 会 scramble，它不会）；**但**环半径 27–30 Å 比 cp233 的 25.4 Å 松 2–5 Å，R146 够不到邻居 WalkerA → **M2 0/5**（single-seq Boltz）。拓扑锁住了、催化界面没锁住。

**旁支（同属生成式，但不是新拓扑）：**
- **Mode A（partial diffusion 加固 cp233）= 负结果**：pT 8/10/12/20 全部把 R146→WalkerA 的 Cα 从 8.1 松到 8.4 Å，从不收紧。cp233 已接近最优,加固没用。
- **ProteinMPNN/LigandMPNN de-novo 序列**（cp233_NOVEL，~53% identity，两预测器 5/5）是在 **CP 骨架**上换序列的“新颖性”结果，不是新拓扑；归到 CP 名下。

### 要你去跑 AF3 的结构（只有一个方法没跑过 AF3）
- **RFdiffusion 铺环**：`gp16_design/outputs/rfdiffusion_modeB/tiled_ring/rfdiff_tiled_ring_for_AF3.fasta`（4 条 1835-aa 序列）。← **发给你了，去跑 AF3**。
- N-C L30、CP cp233 都已 AF3-complete，无需再跑。
- 说明：铺环在 Boltz single-seq 是 M2 0/5（长 shot），但 AF3 用自己的 MSA 可能收紧，值得花 1 个 AF3 job 验证——若 AF3 也 0/5，就正式确认“RFdiffusion 拓扑对但催化不达标”，cp233 是无争议的赢家。

---

## Task 2 —— 除了湿实验，第二个独立的“干实验”验证方法

**核心洞察**：AF3 / Boltz / OpenFold3 **不是真正独立**——它们都训练在 PDB + coevolution 上，共享同一套“学到的先验”。真正独立的干实验必须**不依赖这套学习先验**。只有一类满足：**基于物理力场的方法**。

推荐（按独立性/价值排序）：
1. **★ 分子动力学 MD（物理力场，Amber/CHARMM/OpenMM）** —— 最独立。折叠预测器只给一个静态快照;MD 问的是完全不同的问题:**这个设计的环在物理力场下稳不稳、随时间开不开**。这就是 Task 4。**这是 Task 2 的主答案。**
2. **Rosetta / physics 能量**（FastRelax + interface ΔΔG + packing/clash 评分）——半独立(knowledge-based 力场,和深度学习先验不同源)。便宜,可对每个界面算 trans 界面能。
3. **MM-GBSA / MM-PBSA 界面结合自由能**——对 R146→WalkerA trans 界面算结合能,物理独立。
4. **蛋白语言模型 naturalness（ESM pseudo-perplexity）**——对**序列**打分,与**结构**预测器正交:设计序列像不像天然蛋白?便宜的 sanity check(主要对 de-novo 序列如 cp233_NOVEL 有意义)。

**建议组合**:MD(Task 4,主)+ Rosetta interface ΔΔG(便宜的第二正交信号)。两个都基于物理、都不共享 AF3 的先验 = 真正独立的第二方法。

---

## Task 3 —— 死座点突变实验：**同意暂时不做**

你说得对,单个氨基酸(E119Q,近等排 Glu→Gln)对结构的破坏可忽略——所以“环仍闭合”几乎是必然,信息量低(这和我对 CS §11 的分析一致:E119Q 的“容忍”只是**装配容忍**,不是功能失活;真正有信息的是 R146A 定位对照)。**跳过,不占算力。**

---

## Task 4 —— OpenMM MD:天然 apo / 7JQQ / 最优设计 三体系对比

见下方 Prompt。**这同时就是 Task 2 的执行。**

---

## Task 5 —— ClpX vs gp16:同一个方法找最优,还是各自独立?

**我的判断(要写进论文)**:**框架(方法)是同一个,最优解(拓扑)是每个蛋白各自的。**
- **共享的是框架**:枚举多种单链拓扑(直连 / CP / 生成式)→ 用 interface-resolved M-metrics 在多预测器下打分 → 按 escalation ladder 升级。**这套“fold-in-context”流程对两个蛋白完全一样。**
- **不共享的是赢家**:gp16 需要 **CP**(直连被 AF3 否决);而 single-chain ClpX(Baker & Sauer 2005,PMID 16237435)是用 **直接串联融合 + 接头**做成的——**直连对 ClpX 就够了**。所以两个蛋白的最优拓扑不同。
- **论文点(强)**:贡献不是“某个构建”,而是**可迁移的方法**——同一套框架,对不同环 ATPase 各自找到其最优连接策略(gp16→CP,ClpX→直连)。这正是“泛化”的正确含义:**方法迁移,解不迁移**。要证实 ClpX 这半边,需要真的把框架跑在 ClpX 上(= 我上轮给 CS-B 写的 ClpX 回顾性验证任务)。

---

## 元问题:每个任务用什么工具、在哪个 session 做、Prompt 是什么

| 任务 | 工具 | Session | 谁执行 |
|---|---|---|---|
| **T1 对比 + 分享结构** | Claude Code(已完成) | 本 session | 我已做;**AF3 你跑**(RFdiffusion 铺环 fasta) |
| **T2 设计独立干实验** | Claude Code(已设计)→ 落地即 T4 | = T4 | 我(设计)+ 见 T4 |
| **T3 死座** | —— | —— | **跳过** |
| **T4 MD(apo/7JQQ/设计)** | Claude Code 驱动 **Colab Pro + OpenMM** | **新开一个 session**(算力密集,别和 CS 折叠混) | 我(Claude Code)驱动,或你跑 notebook |
| **T5 ClpX vs gp16** | Claude Code(分析已给)+ **Claude Science**(ClpX 数据) | ClpX = **新开一个 CS session**(新蛋白、干净上下文) | 分析我已给;数据见 CS-B prompt(NEXT_CS_TASKS.md) |

**“Claude 还是 Claude Artifacts”**:全部是 **Claude(Code/Science)** 的计算/分析任务,**不是 Artifacts**。Artifacts(可视化网页)只在你想把“三方法对比 + MD 曲线”做成一页交互对比图时才用,可选、非必需。

---

## 可直接复制的 Prompt

### 【T4 / T2】MD 独立验证 —— 给一个新的 Claude Code(驱动 Colab Pro + OpenMM)session
```
用 OpenMM 在 Colab Pro(A100,别用 Modal)跑 MD,做“独立于结构预测器的物理验证”,比较三体系:
(A) 天然 gp16 环 apo 态(用 7JQQ 去掉 ATP/Mg/pRNA/DNA 的蛋白环,或 apo 预测环);
(B) 7JQQ 天然 ATP-bound 螺旋态(带 ATPγS+Mg);
(C) 我们最优单链设计 cp233_int15_inter10(gp16_design/outputs/structures/af3_sweep/cp233/inter10/model_0)。
每体系:加氢、溶剂化(TIP3P)+ 中和离子、能量最小化、NVT+NPT 平衡、production(算力允许的最长;1750-aa 全原子很贵——
先隐式溶剂或短跑几 ns 看稳定性,再决定是否上显式长跑,或用 CG/Go-model 跑全长)。
读数(全部随时间):整体 Cα-RMSD、per-residue RMSF、环闭合度(每个 sequential 界面 R146→邻居 WalkerA 距离)、
径向对称性(radius_CV)、界面接触数。
判据:设计环(C)的 MD 稳定性/闭合度应接近天然 apo(A);若快速开环/塌陷、界面解离,而 A 稳定 → 设计不合理。
7JQQ(B)作为“上了 ATP 的螺旋/非对称态”参照,看设计更像 apo 对称态还是能响应。
守则:OpenMM;Colab Pro;增量保存轨迹+能量;报 RMSD/RMSF 曲线 + 三体系闭合度对比图 + 一句结论。
先明确算力预算:1750 aa 显式溶剂 MD 极贵,务必先估 ns/day,不够就降到隐式溶剂或粗粒化,别空跑。
```

### 【T5】ClpX 回顾性验证(证明“方法迁移、解不迁移”)—— 给一个新的 Claude Science session
见 `gp16_design/NEXT_CS_TASKS.md` 的 **CS-B**(single-chain ClpX 回顾性验证:已知 good 构建应闭合/coupler 咬合、known-bad 应失败;六聚体 → M2 判据改 ≥5/6;找 ClpX arginine-finger/Walker-A/sensor 残基编号别猜;拿不到 ClpX 就退 T4 gp17)。跑完把 ClpX 的最优拓扑(直连?CP?)和 gp16(CP)并排,写“同框架、不同解”的对比表。

### 【T1】AF3 —— 你自己在网页跑
把 `rfdiff_tiled_ring_for_AF3.fasta` 里的序列(挑 1–2 条,如 d0)提交 AF3(5 models),导出后我用 score_m2 打分,补齐 RFdiffusion 那一格。
