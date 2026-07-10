# T4 gp17 — 单链环框架测试:直连(direct fusion)能闭合 ✓（2026-07-09）

> 框架第 3(4)个蛋白测试。问题:gp17 能否用 N-C 直连 / CP / diffusion 任意一种闭合成单链环?
> **答案:能——DIRECT FUSION 直接闭合**(R245 反式手指 5/5,native-like)。是框架继 ClpX 之后**第二个"直连"蛋白**。

## 结论
- **VERDICT:直连 C→N fusion 闭合一个单链 gp17 ATPase 环。**
- gp17 是 T4 large terminase(病毒 dsDNA 打包 ATPase,和 phi29 gp16 同科)。UniProt P17312。

## 拓扑(native 五聚体 3EZK,34Å cryo-EM)
- 五聚体(和 gp16 同)。**ATPase C 端(res360)离通道 32.8Å——不像 gp16 的 C 端 6.2Å 扎进 DNA 通道**;N 端(res10)13.3Å 靠近通道壁。
- **直连 C→N gap = 34Å**(介于 ClpX ~20Å 和 Rho ~41Å 之间)。
- 功能残基(全部从结构/文献 ground,未猜):Walker-A 161-167 / K166,Walker-B 251-256 / D255,催化羧基 E256,ATP 结合 138/143/202(UniProt P17312 + 2O0H ATP 接触 + Rao lab 文献)。

## 方法(建 + 折)
- **直连构建**:ATPase 域 res1-360 ×5,20-aa GS linker,共 1880aa。
- **折叠**:tiled block-diagonal MSA 经**免费 Boltz-2 NIM**(NIM 通过 polymers[].msa 接受预算 a3m)。另加 3 个对照。

## M2 结果——是否闭合?**YES**
| 构建 | M1 | M2(R245 反式手指→邻居催化 E256) | 备注 |
|---|---|---|---|
| **直连(tiled MSA)** | sequential YES,compact planar(radius 32.7Å,CV 0.01,planarity 0.2Å) | **5/5**(4.5-5.2Å) | 闭合 |
| native 5 链参考(独立折) | 同几何,radius 32.0Å | 5/5(3.7-4.0Å) | 设计 native-like,非 linker 假象 |
| single-seq 对照(同构建无 MSA) | planarity 6.8Å | 0/5(~52Å) | confounded(验证 tiled MSA + 判据) |
| 单体 vs 晶体 2O0H | — | RMSD 1.09Å | pipeline 验证 |

## 诚实的关键点(honest subtlety)
- **gp17 的反式手指是 R245,不是 R162**——R162 是 Walker-A 近端/顺式(→自己 K166 5.9Å)。R245 是 Rao lab 4 个候选(R162/R245/R321/R406)之一;fold + native 参考显示 R245 是伸到邻居活性位点的那个。
- 判据锚定在 **biology-agnostic 的 M1 register + native 环复现**,不是假设某个反式手指。

## 拓扑表新增行
**gp17(五聚体):ATPase C 端 32.8Å 孔外,直连 C→N gap 34Å → DIRECT FUSION(测试:闭合;R245 反式手指 5/5,native-like)。** 框架第 2 个直连蛋白(继 ClpX),把直连扩展到 34Å gap;与 gp16 形成干净对比(同病毒打包 ATPase 科、同五聚体,但 gp16 因 C 端扎 DNA 通道必须 CP)。

## Deferred / 诚实边界
- 不需要 GPU/RFdiffusion(直连已闭合)。便宜的下一步:linker 长度 sweep(10/15/25aa)+ 跨预测器(AF3/OF3)。
- **不写 "validated"**——这是 prediction + native 参考 + 阴性对照,无 MD/湿实验。
- 机读 verdict:`gp17_topology_map.json`;脚本:`topology_analyze.py`, `build_seqs.py`, `get_msa.py`, `run_fold.py`, `score_gp17.py`;结果:`fold_results/`, `structures/`。
