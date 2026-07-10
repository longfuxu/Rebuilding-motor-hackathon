# CS UPDATE prompt — switch directed evolution to GREEDY EPISTATIC ACCUMULATION
_粘贴给正在跑定向进化的 CS session。这条更新的核心:不要独立测单点,要在"当前最优背景"上逐轮叠加。_

---

**方法更新(重要,请从下一轮起改用):directed evolution 必须用 GREEDY EPISTATIC ACCUMULATION,而不是独立单点 + 几个拍脑袋的双突变。**

**为什么**:我们要设计更强的马达,但"更强"往往来自**两个及以上突变的协同(epistasis)**,单点测不出来。独立测单点、或只测几个随机双突变(我们之前只测了 4 个双突变,采样严重不足),既会漏掉能叠加的组合,也会得出"doubles 不 stack"这种被采样不足误导的结论。真正的定向进化(以及 MLDE/ALDE)是**在当前最优变体上逐个累积**。

**算法(每一轮)**:
1. **background = 上一轮的最佳变体**(第 1 轮 background = WT;之后 = 累积到现在的最优,例如 M307W → M307W+X → M307W+X+Y …)。
2. 在这个 background 上,**枚举所有允许的"再加一个单点突变" background+m**(m 遍历允许突变集;不要重复已固定的位点)。
3. 每个 background+m:tiled-MSA 折(免费 Boltz-2 NIM)→ **门控打分** GATES(M1 环 ∧ M2 handedness-robust 5/5 ∧ M5 ATP口袋 ∧ soft-mode 保留[防死砖] ∧ biochem[排除致死]) × power(coordination-PRS + force-network Y129 + M3 grip)。ESM-2 plausibility 过滤(避免 "holes")。
4. **选在该 background 上最优的 background+m 作为下一轮的新 background**。**保 2–3 条平行谱系**(从不同的好起点,如 M307W 力线 和 S99W 力线,防止陷入单一局部最优;必要时用 active-learning/UQ 选下一批)。
5. 累积 **10–20 轮**,或**连续 2 轮 power 增益 < 噪声地板(WT replicate 定的 floor)就停**,报告"到累积几个突变时增益饱和"。
6. 每 3–4 轮把当前最优跑 **MM-GBSA + 交叉预测器(AF3/OF3)** 核实。

**关键点**:每轮的候选是"在当前最优上再加一个",所以到第 k 轮你手里是一个**累积了 k 个突变的变体**,而且每个新突变都是在前面突变**存在的前提下**选出来的——这才真正测"两个/多个突变合在一起是否更强"。

**我们系统(gp16 cp233_int15_inter10)的具体设置**:
- **起点 background**:M307W(round-1 力位点)+ 并行一条 S99W(独立的第二力位点)。
- **允许突变集(可动)**:biochem-TOLERANT — S99、S127、N128、K56、Q222、K233、L289、V100、I130、Phe6*、R53*(*=仅保守替换);rigidify 用 →Phe/Trp/Tyr/Pro。
- **绝对禁动(biochem 致死/近致死)**:Y129、E58、K105、K328、R330、I28、Y32、S106、N126、N158、R234、K294;催化残基 Walker-A/B + R146 也不动。
- **判据脚本**:复用 `gp16_design/outputs/directed_evolution/score_variant.py`(门控 + power proxy);打分永远用 M1–M5,**绝不看 pTM/pLDDT**;M2 必须 handedness-robust。
- **产出**:每轮 ranked CSV + 累积曲线(power vs 累积突变数)+ 最终 **5–10 个"实验可测的更强变体"清单**(累积突变、预测 coord/force 增益、过噪声地板否、soft-mode 保留、跨预测器一致、MM-GBSA、机制一句话)。诚实报告:哪条轴(force/coordination)能稳过噪声,饱和在几点突变。
