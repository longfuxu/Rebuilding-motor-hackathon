# CS 三 session 进度盘点 + 下一步最有科学意义的方向 + 三个 prompt（2026-07-09）

## 三个 session 现状（都已完成各自任务）
1. **`9c944ce5` gp16 主线**:cp233 三预测器 lead + 全蛋白 CP 系统 screen(C-domain 是闭合区,cp285 第二位点)+ E119Q 死座可寻址(§11:各座都容忍、R146A 定位对照)。
2. **`95d796ed` de-novo + ClpX**:cp233_NOVEL(~53% identity,10 条两预测器过关)+ ClpX 回顾性验证(good 6/6、R307A/E lesion 0/6)。
3. **`e1736e32` ClpX 拓扑**:ClpX vs gp16 拓扑对比("同框架、不同解":gp16→CP、ClpX→直连)。

## ★ 接下来最有科学意义的一件事:从"能装配"迈到"有功能"
**判断:结构那条线已经很扎实了**(闭合、可寻址、系统性、能泛化)。**最薄弱、也最值钱的缺口是"功能"——目前所有证据都是 apo + 装配层面**。论文的标题主张是"一个基因可寻址的**功能性**环马达",但我们还没碰"功能"。所以最有科学意义的下一步 = **给它加一层功能/机制证据**:
- 死座的"可寻址"现在只是**装配容忍**(E119Q 近等排,环仍闭合几乎必然)。真正有信息的是:**上了 ATP/Mg 之后,E119Q 是否在它那个座局部破坏催化口袋**(Mg/攻击水配位),而邻座完好?→ 这把"装配容忍"升级成"位置特异的催化口袋破坏"= 真正的可寻址失活证据。
- 加上一个**并行的方法学收尾**:single-seq 折叠会把已验证的 cp233 都压到 0/5(confound),所以生成式三方法(直连/CP/RFdiffusion)的公平对比**必须用 tiled MSA**重折。

下面三个 prompt = 三条并行的最高价值任务,分给三个 session。守则同前:免费 NIM、sequential M2、apo 之外用 ligand-loaded、不写 "validated"、增量保存。

---

### Prompt ① 给 `9c944ce5`(gp16 主线)——死座的**功能定位**(ligand-loaded,最高价值)
```
在 cp233_int15_inter10 上做 E119Q 死座的"功能定位"实验(不是装配层面):
把每个亚基的 ATP·Mg(用 7JQQ 的 ATPγS+Mg 或 ADP·Mg,template OFF,和之前三态对照同法)放进催化口袋,
折 WT、E119Q{seat3}、E119Q{all5} 三个构建(Boltz-2 + OpenFold3,tiled MSA),per-seat 量:
(a) ATP/Mg 口袋是否成形:Walker-A K30→ATP β/γ-磷酸距离、Walker-B(D118/E119)→Mg/攻击水配位;
(b) E119Q 是否**在它那个座**把 Mg/攻击水配位打乱(距离显著变大/配位数掉),而**邻座保持完好**;
(c) R146 反式手指 vs 邻座 ATP 是否仍咬合(M2)。
判据:E119Q 座的催化口袋(Mg/water)局部塌、邻座不塌 = 位置特异的催化失活(不只是装配容忍)。
这是把"可寻址"从装配层面升到功能层面的关键证据。cp233 within-copy 位点:R146=255、Walker-A=133-140、E119=228。
最好的 1-2 个发用户跑 AF3。守则:sequential M2、ligand-loaded、增量保存、不写 validated。
```

### Prompt ② 给 `95d796ed`(de-novo/泛化)——**T4 gp17** 第三个蛋白(扩泛化)
```
把 fold-in-context + 拓扑升级阶梯框架应用到第三个环 ATPase:噬菌体 T4 large terminase gp17
(病毒 dsDNA 打包马达,和 phi29 gp16 同科——同科泛化,比 ClpX 跨科更直接)。
步骤(和 ClpX/gp16 一样,残基编号从结构/UniProt 查证,别猜):
1. 取 T4 gp17 功能结构(cryo-EM 打包复合物)+ UniProt 序列;定位 Walker-A、Walker-B、trans arginine finger、寡聚态。
2. 拓扑分析:末端离 DNA 通道多远(堵不堵)+ 直连 C→N gap 多长 → 预测最优拓扑(直连/CP/生成式),
   和 gp16(CP)、ClpX(直连)、Rho(测试中)并排进拓扑表。
3. 建对应最优拓扑的单链构建,tiled MSA 折(Boltz+OF3),score(gp17 版 M2:arginine finger→邻居 Walker-A,
   floor 按寡聚态定 ≥n-1/n);good 闭合 = 框架第三例通过。
产出:gp17 拓扑 + fold 结果 + 更新"四蛋白(gp16/ClpX/Rho/gp17)× 三方法"拓扑规则表。守则同上。
```

### Prompt ③ 给 `e1736e32`(拓扑/折叠)——**tiled-MSA 公平重折**生成式构建(修 confound)
```
关键方法学问题:single-seq Boltz 把已验证的 cp233 都折成 M2 0/5(tiled MSA 下是 5/5),
所以任何生成式大单链环用 single-seq 打分都不可信。用 tiled block-diagonal MSA(给 cp233 5/5 的那套 cofactor.fold)
重折下面这些,做**公平的三方法对比**:
1. gp16 RFdiffusion 铺环(sal_L50,1835aa,gp16_design/outputs/rfdiffusion_modeB/tiled_ring/):single-seq 是 0/5,
   tiled MSA 是否翻?
2. Rho 三方法(直连 1640aa / CP330 / CP390,gp16_design/outputs/rho_generative/):都要 Rho 专属 tiled MSA
   (对 Rho motor domain res175-414 跑 MMseqs2 homolog → block-diagonal tile 6×);M2 用 R366→邻居 WalkerA(K184)。
3. cp285/cp297(gp16 备选 CP 位点,AF3_SEQUENCES_TODO.txt)。
每个用 Boltz-2 + OpenFold3,score_m2(M1 sequential + M2)。产出:统一的"蛋白 × 方法 × tiled-MSA M2"对比表——
这是论文三方法对比的核心数据(把之前 single-seq 的 confounded 结果替换掉)。守则:sequential M2、增量保存。
```

## 我(Claude Code)并行在做 / 已完成
- **✅ MD 物理验证已完成**(CS session 跑的,GBSA-OBC2):cp233 像天然 apo 一样保持闭合+对称,不同于 7JQQ ATP 螺旋态 → 通过独立于预测器的物理稳定性检验(`md/openmm_validation/`)。**这就是"功能/机制"那条线的第一条腿**——所以 Prompt① 的 ligand-functional-locality 是自然的下一步(第二条腿)。GCP 重跑 MD 因此取消(已有结论)。DCD 全轨迹没存,真·轨迹动画待补(见 PROJECT_STATUS §6)。
- **Rho 生成式 RFdiffusion**:正在 Colab 上跑(motif res175-414,连接件 30-52,跨 46Å gap);出 backbone 后交 Prompt③ 的 tiled-MSA 折叠。
