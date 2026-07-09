# Rho 三方法 fold 测试 + 一个关键方法学发现（2026-07-09）

目的：把三种单链方法(N-C 直连 / CP / RFdiffusion 生成式)在 Rho(第二个"gap 长、末端在孔外"的六聚体环马达)上真的 fold 出来对比,理解各自 limitation/advantage(用户要求:即使 CP 能折也要测 diffusion)。
Rho = E. coli 转录终止因子(P0AG30, 419aa, 六聚体)。M2 = 反式 arginine finger **R366** → 邻居 Walker-A(**K184**,P-loop 179-186),≥5/6。residue 均从结构/UniProt 确认,未猜。

## 已 fold(Boltz single-seq NIM)—— 但结果被 confound
| 构建 | 长度 | M2 | M1 |
|---|---|---|---|
| Rho 直连(6×motor res175-414, (GGGGS)×8) | 1640 | **0/6**(手指 ~24Å 偏) | 紧凑环、register YES |
| Rho CP330(join gap 16Å off-pore) | 1580 | **0/6** | 紧凑环、register YES |
| Rho CP390(join gap 15Å off-pore) | 1580 | **0/6** | 紧凑环、register YES |
| **对照: gp16 cp233(已知好, tiled-MSA 下 5/5)** | 1750 | **0/5** ← | 紧凑环 |

## ★ 关键方法学发现(good finding)
**单序列(single-seq)结构预测无法评估这类大单链环——连已验证的 cp233 在 single-seq 下都掉到 0/5。**
所以上面 Rho 三个 0/6 **不能区分方法优劣**,是被 single-seq 的弱先验压平的,不是构建本身失败。
→ **公平对比必须用 tiled MSA**(就是让 cp233 拿到 5/5 的那套 cofactor.fold + block-diagonal tiled MSA)。single-seq 只能当"能不能折成环"的粗筛,不能当 coupler-engagement 的判据。这条直接写进论文的"方法 limitation"。

## 还没做 / 下一步(要 tiled MSA)
1. **RFdiffusion 生成式 Rho 连接件**:motif 已建(`rho_motif_AB.pdb`,2 个相邻 motor domain,C→N gap 46Å)。跑 RFdiffusion de-novo connector → ProteinMPNN → 铺 6× → 得生成式 Rho 环。(需 GPU;MD 占用中,排队。)
2. **三个 Rho 构建(直连/CP/生成式)全部用 tiled MSA 折 + score** —— 这才是公平对比。需要 Rho 专属 tiled MSA(对 Rho motor domain 跑 MMseqs2 homolog → tile 6×),用 cofactor.fold 或交给 CS 的 tiled-MSA 管线。
3. 判据同 gp16:M1 sequential=YES + M2 ≥5/6(≥2 预测器)。看直连/CP/生成式哪个在 Rho 上闭合 = 三方法在"gap 长、末端外周"体系上的优劣。

数据:`rho_direct.fasta` `rho_cp.fasta` `rho_motif_AB.pdb` `rho_scores.json` `cp_scores/`。
