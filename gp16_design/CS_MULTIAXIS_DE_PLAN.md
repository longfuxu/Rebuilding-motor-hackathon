# CS 多轴定向进化计划 —— 先 stiffness,再 grip / power-stroke / coupling geometry
_一整晚跑。粘贴给正在跑定向进化的 CS session(或新开一个)。**核心:顺序执行,一相收敛+写报告再下一相,别混别急。**_

---

## 大原则(所有相通用)
- **顺序,不并行,不着急**:先跑完当前的 stiffness 优化(Phase 0),两条线都收敛并写完报告,**再**进 Phase 1;之后一相一相做,每相独立产出。**四个轴优化的是不同的目标,别弄混。**
- **每一相都用 GREEDY EPISTATIC ACCUMULATION**:固定当前最优变体为背景 → 枚举所有允许的"+1 突变" background+m → tiled-MSA 折 + 打分 → 选最优 background+m 作下一轮背景 → 累积;**连续 2 轮 power 增益 < 噪声地板(WT replicate 定)就停**。保 2–3 条平行谱系防局部最优。
- **共享硬门槛(GATES,每相都必须过,不随相变)**:M1 闭环 ∧ M2 handedness-robust(反式手指没脱)∧ M5 ATP 口袋完整 ∧ biochem(排除 SEVERE/LETHAL)∧ 构建能折。**只有"优化目标 objective"随相不同。**
- **绝对禁动(biochem 致死/近致死)**:Y129、E58、K105、K328、R330、I28、Y32、S106、N126、N158、R234、K294;催化 Walker-A/B + R146 也不动。可动优先用 biochem-TOLERANT(S99、S127、N128、K56、Q222、K233、L289、V100、I130、F6*、R53*)+ 各相相关区域的耐受位点。
- **每相赢家跨预测器核**(Boltz-2 + OpenFold3;无 AF3 权重就 OF3;重计算 MM-GBSA/显式 MD 可上 GCP:Spot + 用完删 VM + credit 计费)。
- **判据脚本复用** `gp16_design/outputs/directed_evolution/score_variant.py`(门控 + 各 proxy);打分**永远 M1–M5,绝不看 pTM/pLDDT**;M2 必须 handedness-robust。功能残基/接触集见 `reproduce/functional_contacts_7jqq.json`。

## 这个计划要回答的关键科学问题
同一轴(rigidification/刚度)在 ~1 个突变就饱和了(M307W;N128W 只 within-noise,不 stack)。**不同轴(grip / power-stroke / coupling geometry)会不会各自独立提升、而且彼此能叠加?** 能叠 → 多轴组合出真正更强的马达;每轴也饱和 → proxy 引导的点突变接近上限(诚实结论)。**这才是"设计更强马达"该问的下一步。**

---

## Phase 0 —— 先跑完 stiffness(力传导刚度)优化【正在跑,收尾】
把正在跑的 greedy 累积(M307W 线 + S99W 线)**跑到两线都收敛**。写报告:每线累积曲线 + 最优 + "刚度轴饱和在几点突变"的结论。**这一相彻底结束、报告写完,才进 Phase 1。** 记下 **STIFF_BEST**(跨预测器确认后的刚度轴最优变体及其关键突变)。
- 目标(已在用):force-network(Y129 betweenness)+ PRS coordination。
- 预期(诚实):很可能确认"刚度不 stack,~1 突变封顶"。

## Phase 1 —— DNA/substrate-grip 轴（"抓得更牢",不是"更硬")
- **目标(最大化)**:抓 DNA 的能力 —— 把 7JQQ 的 dsDNA 叠进设计通道,量 **gp16–DNA 界面结合能(MM-GBSA)更有利** + M3 通道抓手残基更好地朝孔内咬 DNA 骨架。
- **允许突变**:通道内衬 / DNA 接触面**附近**的 biochem-TOLERANT 位点;可在朝孔内的中性/耐受位点**加抓手**(→K/R 咬磷酸骨架)。DNA 接触集:55–60,82–83,98–100,125–130,292–293,330(**但这些里的 lethal 不动**)。
- **门槛**:M1∧M2∧M5∧biochem ∧ **通道仍通(M3 孔径够 dsDNA ~20Å,别把通道堵死/抓死)**。
- 赢家跨预测器 + MD 稳定性核。产出 **GRIP_BEST**。

## Phase 2 —— power-stroke amplitude 轴【⚠ 和刚度相反:这一轴是"更能动"】
- **目标(最大化)**:power-stroke 幅度 —— ENM/ANM 最低频**功能模**(环开合 / planar↔helical 易位坐标)的**幅度 / 与易位坐标的 overlap**。马达靠大构象变化 translocate DNA,这一轴要**增强功能软模**。
- **⚠ 别和 Phase 0 弄混**:Phase 0 里 soft-mode 是"别低于 native"的**门槛(gate)**;这里 soft-mode(功能模幅度)是**优化目标(objective)** —— 要尽量**高于** native,同时保持耦合。
- **允许突变**:hinge / lever 铰链区(域间连接、power-stroke 支点)的 biochem-TOLERANT 位点 —— 让铰链更柔顺、行程更大(如在过刚铰链引 Gly、去掉限制行程的接触)。**不是全局变柔**,是定点松开功能铰链。
- **门槛(关键张力:更能动但仍耦合)**:必须仍 M2 耦合 ∧ M1 闭环 ∧ M5 ∧ biochem。太柔会脱耦,门槛会挡住。
- 赢家跨预测器 + MD(看是否真出现更大幅度的功能运动)。产出 **STROKE_BEST**。

## Phase 3 —— coupling geometry 轴（催化偶联几何预组织)
- **目标(最优化)**:反式 R146 → 邻居 Walker-A 的**几何预组织** —— 让手指更好地咬进邻居 ATP 口袋(M2 距离更短更一致、sequential register 更干净)。
- **允许突变**:界面上 R146 / 邻居口袋**周边**的 biochem-TOLERANT 位点(优化手指 reach/register 的);**不动 R146 / Walker 本身**。
- **门槛**:M1∧M5∧biochem ∧ soft-mode 保留(别为了咬紧而僵死)。
- 赢家跨预测器核。产出 **GEOM_BEST**。

## Phase 4 —— 多轴组合【核心检验:不同轴能不能叠加】
取 STIFF_BEST + GRIP_BEST + STROKE_BEST + GEOM_BEST 各自的关键突变(它们动的是**不同残基**),**组合**(两两、三三、全组合)→ tiled-MSA 折 + 全门槛 + **四个轴的 proxy 都量**。问:**同轴不能 stack,不同轴能不能?** 过所有门槛、且 ≥2 个轴 proxy 同时稳过噪声的组合 = 多轴更强马达候选。
- 产出**最终 5–10 个"实验可测的更强马达"清单**,每个标:累积突变、各轴 proxy 增益(过噪声否)、跨预测器一致、MM-GBSA、机制一句话。

## 产出与预算
- 目录 `gp16_design/outputs/directed_evolution/multiaxis/`:每相 ranked CSV + 收敛曲线 + `PHASE{0..4}_REPORT.md`,外加 **`FINAL_REPORT.md`**(四轴各自结果 + 组合 + 最终清单 + **诚实结论**:哪些轴可动、哪些饱和、不同轴能否叠、最终能否稳定超 native)。
- 免费 Boltz-2/OF3 NIM;每变体 ~60s。一整晚:跑完 Phase 0 + Phase 1–3 各几轮 + Phase 4 组合足够。重计算上 GCP(Spot+删+credit)。
- **再次强调:一相收敛 + 写完该相 REPORT,再进下一相。不要并行,不要跳步,不要把 Phase 2 的"越软越好"和别相的"soft-mode 门槛"搞混。**
