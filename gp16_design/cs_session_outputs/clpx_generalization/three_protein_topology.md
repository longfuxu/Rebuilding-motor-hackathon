# 三个环 ATPase × 三种单链方法 —— 方法迁移、解不迁移(2026-07-09)

同一框架(fold-in-context + interface-resolved M2 + escalation ladder 直连→CP→生成式),
按序试到能用为止。三个蛋白落在三个不同的最优——正好让三种方法都有合理的"最优蛋白"。

| 轴 | ClpX | gp16 | Rho(本轮新测,3ICE) |
|---|---|---|---|
| 家族 | AAA+ HslU/ClpX unfoldase | ASCE 打包马达 | RecA-type RNA 易位酶 |
| 寡聚 | 六聚体 | 五聚体 | 六聚体 |
| 末端离通道轴 | 47–65 Å(孔外周) | C端 res330 **扎在 DNA 通道**(6.2 Å) | 57–60/40–45 Å(**孔外周**) |
| 直连 C→N gap | **~20 Å** | 53–66 Å | 65 Å(全长)/ **41 Å(motor domain)** |
| 直连是否堵通道 | 否 | **是** | 否 |
| 最优拓扑 | **直接融合** | **环状置换(CP)** | **生成式 de-novo 刚性连接件(预测)** |
| 为什么 | gap 短 + 末端在孔外 + 有近-WT 单链先例 | C端堵 DNA 通道、gap 又长 → CP 把接头挪出通道并缩短 | 末端在孔外(CP 的"挪出通道"优势用不上)**但 gap 41–65 Å 太长**,软 linker 会 scramble(如 gp16 L40+)→ 需要**刚性 de-novo 连接件**跨 41 Å |

**核心结论(论文点)**:三种连接策略各自是某个蛋白的最优,判据是**两条几何量**——(1)末端相对通道的位置(堵不堵),(2)直连 gap 长度。ClpX(短、外周)→直连;gp16(长、堵通道)→CP;Rho(长、外周)→生成式刚性连接件。**这就补齐了"生成式必需"的那一格,三种方法才都站得住。**

**待确认(fold 证实)**:把 Rho motor-domain(~res150–414 ×6 ≈ 1400 aa)分别做成 (a) 直连+GS linker、(b) CP、(c) RFdiffusion de-novo 连接件,折叠打 Rho 版 M2(trans arginine finger→邻居 Walker-A),证明只有 (c) 闭合。Rho 催化残基需查(RecA Walker-A P-loop ~179–186、arginine finger)。这是 CS-C 的具体化。

数据:`rho_3ice.pdb`(RCSB 3ICE),环序 D-C-B-A-F-E,孔半径 ~7 Å。
