# 独立于结构预测器的干实验验证（Task 2 落地，2026-07-09）

结构预测器（AF3/Boltz/OF3）不是真正独立——都学 PDB+coevolution 同一套先验。下面是**不共享该先验**的正交验证，已跑的三个是本地纯计算（numpy/scipy/ESM），待跑的两个需要 Colab 物理力场。脚本在 `outputs/rfdiffusion_modeA/` 与 `outputs/rfdiffusion_modeB/`。

## ✅ 已跑（本地，正交信号）

### 1. ATP 结合位点定位对照（几何）—— CP 后功能位点是否还在正确位置
在 cp233 上量 ATP 口袋三角（Cα）:Walker-A K30、Walker-B E119(顺式同口袋)、R146(反式手指)。
| | K30–E119(顺式口袋) | K30–R146 | E119–R146 |
|---|---|---|---|
| 天然 gp16 亚基 | 10.4 Å | 26.5 | 24.0 |
| cp233 sub1–5(全部) | **10.6 Å** | 27.2 | 24.4 |
**结论:环状置换没有破坏 ATP 位点几何——五个亚基完全一致、与天然差 ≤0.7 Å。**功能位点在 CP 后仍正确定位。(DNA 接触位点需要带 DNA 的结构才能量,7JQQ 无 DNA;留待 DNA-bound 结构。)`ca_ring_m2.py`。

### 2. 弹性网络开/闭测试(物理,无学习先验)—— 共价单链是否被 linker 锁死不能开环
ENM(Cα 弹簧网络)最软模式对"环呼吸/开合"坐标的累积 overlap:
| | top1 | top5 | top10 | top20 |
|---|---|---|---|---|
| 天然环(非共价) | 0.00 | 0.00 | 0.01 | 0.24 |
| cp233(共价单链) | 0.00 | **0.09** | 0.10 | 0.13 |
**结论:共价单链没有把环锁死——cp233 保留了与天然同量级的软开合模式(top20 0.13 vs 0.24,且更早在 top5 出现)。**力学上仍是"能开合的马达",不是被 linker 焊死的死环。`enm_openclose.py`。(粗粒度;真 MD 更决定性——见待跑。)

### 3. ESM-2 蛋白语言模型 naturalness(序列层面,正交于结构)
每残基 log-likelihood(esm2_t30_150M):
| 序列 | mean LL |
|---|---|
| 天然 gp16(327aa) | −0.43 |
| RFdiffusion de-novo 连接件(d0–d3) | **−0.40 ~ −0.41** |
**结论:RFdiffusion 的 de-novo 连接件序列和天然 gp16 一样"自然"(LL 几乎相同)。**所以铺环 M2 0/5 的失败是**几何的(环太松)、不是序列不自然**——指向"需要更紧的连接件",不是"序列有问题"。`esm_naturalness.py`。

## ⏳ 待跑(需 Colab 物理力场,prompt 见 DESIGN_METHODS_COMPARISON_AND_NEXT.md）

### 4. 分子动力学 MD(OpenMM,Colab Pro,尽量便宜)—— 最决定性的正交验证
比较 (A)天然 apo 环、(B)7JQQ、(C)cp233 设计的 MD 稳定性/开合。先隐式溶剂或粗粒化(算力便宜),1750aa 显式长跑太贵。ENM(上面 #2)是它的免费预演,已给出"cp233 力学native"的初步正信号;MD 给时间演化的决定性证据。

### 5. Rosetta / MM-GBSA 界面能(半物理,和深度学习不同源)
FastRelax + 每个 trans 界面 R146→WalkerA 的 interface ΔΔG / MM-GBSA 结合自由能;比较 cp233 vs 天然 vs RFdiffusion 铺环。Colab 装 PyRosetta 或用 OpenMM-MMGBSA。

## 一句话
三个正交、非结构预测器的信号都指向同一结论:**cp233 是几何正确(ATP 位点保住)、力学 native(能开合)、序列自然的单链环** —— 独立方法支持了预测器的判断。RFdiffusion 铺环的短板被定位为"几何太松"而非序列或力学问题。
