# 判据体系(M1–M5)+ 论文 Limitations —— 2026-07-09

回应:(a) 除 M1/M2 外需要 M3/M4(与 pRNA、DNA 的接触位点等更重要的功能指标);(b) 把所有 limitation 记进论文。

## A. 判据体系:从"装配"到"功能"的五层
> 判据都用几何/物理,**永远不看 global pTM/pLDDT**(见 Limitation #2)。功能接触残基从**实验结构 7JQQ**(gp16 五聚体 + dsDNA + pRNA 五聚体全复合物)ground,存 `reproduce/functional_contacts_7jqq.json`。

| 指标 | 问的问题 | 判据 | 状态 |
|---|---|---|---|
| **M1 装配** | 是不是一个闭合、按顺序的环? | 五亚基 compact planar 环 + sequential_consistent | ✅ 有(score_m2.py) |
| **M2 催化偶联(灵魂)** | 每个催化手指搭对下一个亚基的油门没? | trans R146→邻居 Walker-A <8Å,**handedness-robust**,沿环 k→k±1 | ✅ 有(需 handedness-robust,见 Limitation #3) |
| **M3 DNA 通道能力(打包马达最关键)** | 中央孔能过 dsDNA 吗?抓 DNA 的残基还在位吗?linker/末端堵不堵孔? | (a) 20 个 DNA 接触残基**在位且朝向孔内**;(b) 中央孔径够 dsDNA(~20–24Å);(c) CP 切点/linker/末端**不落在也不堵** DNA 接触面 | ⏳ 提议;残基已 ground |
| **M4 pRNA 装配面** | 装配到 pRNA 上的表面还完整、还在外侧吗? | 29 个 pRNA 接触残基 surface-exposed + native-like 位置,没被 CP 搬家/linker 埋掉 | ⏳ 提议;残基已 ground |
| **M5 ATP 口袋完整** | ATP/Mg 口袋几何对吗? | Walker-A(24–31)/Walker-B(D118/E119)+ Mg/攻击水配位 native-like | ✅ 部分有(ATP-geom 独立检查) |

**M3/M4 的 ground 数据(7JQQ 4.5Å 接触,native 编号)**:
- **DNA 接触(20,中央通道抓手)**:55–60, 82–83, 98–100, 125–130, 292–293, **297**, 330。
- **pRNA 接触(29,装配面)**:10–17, 37–45, 148–152, 182, 237–240, 254, 267–273。

**M3 已经给出一个即时区分**:**cp233 切点 233 不碰 DNA 也不碰 pRNA(干净);但 cp297 切在 297 = 一个 DNA 接触残基上** → cp297 有破坏 DNA 抓手的风险。这正是 M3 的价值:**cp233 > cp297 不只是因为 M2,还因为切点不动功能面。** cp233 的 trim(去掉 native 1–3/331–332)也不碰任何接触残基。

**优先级**:打包马达的两个"功能灵魂"是 **M2(催化)+ M3(过 DNA)**;M4(pRNA)管装配;M5 与 M2 部分重叠。→ 建议筛选/评估时:M1 门槛 → M2 handedness-robust → **M3 DNA 通道(下一个要建的)** → M4/M5。

## B. 论文 Limitations(全部记录,诚实边界)
1. **单序列预测失效(single-seq confound)**:预测器对这类大单链环用单序列会假阴性——连已验证的 cp233 都掉到 M2 0/5(tiled MSA 下 5/5)。**必须用 tiled MSA**;single-seq 只能当"能不能成环"的粗筛。
2. **pLDDT/pTM 陷阱(实锤)**:RFdiffusion d0 铺环的界面 pLDDT 全场最高(81)却 M2 最差(0/5,arginine finger 在 4/5 拷贝被 MPNN 改没)。→ **高置信度 ≠ 界面对**;只能用 M1/M2 几何判据。
3. **Handedness 简并**:近 C5 环在采样里会随机绕两个方向(镜像),方向特异的 M2 会把"镜像绕向但闭合"的环误判成开环 0/5(E119Q 就栽在这)。→ **共价环 M2 必须 handedness-robust(同时读 k→k+1 和 k→k-1)**。
4. **生成式这一级是框架硬边界**:RFdiffusion 能生成环形骨架,但**同时满足 sequential register + 保住催化残基**很难。gp16 铺环 = **predictor-split**(Boltz+tiled 5/5,AF3 2/5 scrambled,同序列),即"一个预测器认可、另一个不认可"的边际设计 —— 对比 cp233 三预测器一致。这是升级阶梯的真实上限。
5. **预测器非真正独立**:AF3/Boltz/OF3 都学 PDB + coevolution 同一套先验 → 跨预测器一致是 **filter,不是 validation**。真正独立的信号靠物理(MD、MM-GBSA、ENM)、序列(ESM)、几何(ATP-pocket)。
6. **MD 范围**:隐式溶剂 + 短单轨迹(1–3 ns);显式 TIP3P + replicas 是 Phase-2(需 CUDA OpenMM,~10× 成本)。判据在共享窗内已解析,但非长时间尺度。
7. **全计算,无湿实验**:所有结论是 **prediction + native 参考 + 阴性对照**,**不写 "validated"**;ATPase 活性、打包、装配需湿实验。
8. **预测器不能建 ATP 驱动的开环动态**:apo/ADP/ATP 预测都闭合;环的开合机制要靠 MD 或实验,不是 predictor。
9. **M1/M2 是几何 proxy**:抓的是装配 + 催化偶联,**不直接抓 DNA 易位、pRNA 装配、ATPase 动力学** → 需 M3/M4/M5 + 湿实验补。

## C. 来源 & 复现
- 功能接触残基:`reproduce/functional_contacts_7jqq.json`(7JQQ A–E vs F,G/K–O,4.5Å)。
- M1/M2 脚本:`reproduce/score_m2.py`;handedness-robust 版:`pipelines/tiled_msa_fold/`。
- 本文并入 `METHODOLOGY.md`(主文档)的 Limitations 与判据章节。
