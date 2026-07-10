# 两个 Claude Science 通宵 Session prompts（2026-07-09）

粘贴到两个独立的 CS session。守则(两者通用):免费 NIM 折叠(Boltz-2/OF3);**tiled MSA**(大单链环必须,single-seq 假阴性);**handedness-robust M2**(镜像绕向别误判);判据用 M1–M5(`gp16_design/reproduce/score_m2.py` + `pipelines/tiled_msa_fold/`),**绝不看 global pTM/pLDDT**;不写 "validated";增量保存+commit。复用已建管线,别重造。

---
## SESSION A —— 定点进化"更强马达"(10–20 轮,产实验候选)
你在延续一个已跑到 round-2 的 in-silico 定向进化引擎(`gp16_design/outputs/directed_evolution/`,读 REPORT.md + score_variant.py + rank_variants.py 复用)。目标:**跑 10–20 轮,找出最可能提升"马达功率"的突变位点,产出 5–10 个可直接进湿实验的、更 powerful 的 cp233 变体。**

打分(已建,复用):`score = GATES(M1环 ∧ M2 handedness-robust ∧ M5 ATP口袋[含 E58/K105/S106/N158] ∧ soft-mode保留[防死砖] ∧ biochem[排除 SEVERE/LETHAL,见 outputs/excel_mutation_analysis/REPORT.md]) × power(coord-PRS + force-network Y129 + M3 grip + MM-GBSA界面能)`。

每轮:(1) 从上一轮最优 + biochem-TOLERANT 位点(Phe6*/R53*/K56/S99/**S127**/**N128**/Q222/K233,*=保守)+ round-1/2 赢家(**M307W/M307F**、N128F、L289P)提出突变(单点+双点+三点);(2) 用 ESM-2 plausibility 过滤(避免 "holes",见 lit/insilico_directed_evolution_methods.md);(3) tiled-MSA 折 + gated 打分;(4) 选 top-k(贪心 + 保多样性)进下一轮;(5) 每 3–4 轮把当前 top-3 跑 MM-GBSA + 交叉预测器(AF3/OF3)核。
**硬约束**:绝不动 Y129/E58/K105/K328/R330/I28/Y32/S106/N126/N158/R234/K294(biochem 致死);productive 位点是 Y129↔M307 clasp 的 **M307** + 旁边的 **S127/N128**(不是 Y129 本身)。
**产出** `outputs/directed_evolution/overnight/`:每轮 ranked CSV + 收敛曲线(power vs round)+ **FINAL: 5–10 个变体的实验清单**(突变、预测的 coord/force 增益、过噪声地板否、soft-mode 保留、跨预测器一致否、MM-GBSA)+ 每个变体"为什么更强"的机制一句话。诚实报告增益幅度(round-1 单突变增益温和;看多点/多轮能否稳过噪声)。

---
## SESSION B —— 扩 AAA+/ASCE 家族筛选(三方法 × 每蛋白;**务必找到一个能闭合的 RFdiffusion 生成式设计**)
我们只有 4 个蛋白(gp16→CP、gp17→直连、ClpX→直连、Rho→测试中),且**还没有一个 RFdiffusion 生成式设计能闭合成环**。目标:**系统筛更多 AAA+/ASCE 环 ATP酶,对每个跑三方法(N-C 直连 / CP / RFdiffusion),用 M1/M2 找出各自能闭合的方法;重点找到 ≥1 个 RFdiffusion 生成式能闭合的环,交湿实验。**

起点:复用 buildability atlas(`outputs/buildability_atlas/`,读 ATLAS_REPORT.md + descriptors.csv;2-分支规则 + 拓扑描述子管线)。**扩 panel**(取有环结构的):FtsK/SpoIIIE、RuvB、katanin/spastin、Vps4、p97/VCP、proteasomal Rpt、DnaB、T7 gp4、HslU、Lon、ClpA/ClpB、dynein AAA、NSF、Pex1/6、MCM、ORC……(从 PDB 生物学装配体取,残基/末端从结构查证,别猜)。
每个蛋白:(1) 算拓扑描述子(末端-通道距离、直连 C→N gap、寡聚态、末端堵不堵通道)→ 2-分支规则预测方法;(2) **三方法都建**:直连(linker 按 gap)、CP(切在低破坏 loop、避开功能残基——用 M3 式接触检查)、**RFdiffusion 生成式连接件**(固定相邻两亚基 motif、扩散刚性 connector、ProteinMPNN 设序、铺环);(3) tiled-MSA 折 + handedness-robust M1/M2 + M3(通道)打分;(4) 记录每蛋白"哪个方法闭合"。
**⚠ 关键教训**(见 atlas):**apo 单折看不见堵孔——直连的已知错构建也会假闭合**;所以生成式/直连的闭合必须叠 **M3(通道通畅 + 功能残基不被扰)** 才算数,并尽量跨预测器(AF3)核。
**重点交付**:**≥1 个 RFdiffusion 生成式设计的、M1/M2 闭合 + M3 通道 OK + 跨预测器一致的环** → 给它一段 wet-lab-ready 说明(序列、方法、预测判据)。
**产出** `outputs/aaa_ascE_screen/`:panel 描述子表 + 每蛋白×三方法的 M1/M2/M3 结果表 + SCREEN_REPORT.md(哪些蛋白哪种方法成、生成式成功案例、拓扑规则是否扩展/被证伪)+ 生成式赢家的实验清单。
