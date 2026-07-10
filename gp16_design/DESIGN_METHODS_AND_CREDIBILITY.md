# 设计方法学 + 可信度（linker/CP 系统筛选 + 独立验证）— 2026-07-09

回答"两个 linker 长度、CP 切点是否系统筛选过?还需不需要再优化?"以及"cp233 到底可不可信"。

## 1. 两个 linker 是什么
cp233 单链环 = 5 份"环状置换(CP)后的 gp16"用 linker 串成一条链。每份需要**两个** linker:
- **int(NC linker / CP linker)**:环状置换在天然序列内部切一刀,产生新的 N/C 端,再用一个 linker 把**天然的旧 C 端(332)接回旧 N 端(1)区**,让单个亚基仍然折成天然折叠。构建名里的 `int15` = 这个 linker 15 aa(GS 富集:GGGGSGGGGSGGGGS)。
- **inter(subunit linker)**:把相邻的两份亚基共价连起来(前一份的尾 → 后一份的头)。`inter10` = 10 aa。

## 2. 都系统筛选过了 ✓（答案:是,不用现在再大规模优化)
### CP 切点:20 个位点系统筛选(`cs_sessions/session1_main_gp16/outputs/cp_screen_site_summary.csv`)
- 全蛋白扫了 **20 个候选切点**,每个打两个分:**penalty(切点处结构破坏度,越低越好——优先切在 loop、不打断二级结构)** + **是否闭合(跨预测器)**。
- 结果:**233 胜出(penalty 13,最低;闭合 3/3 预测器)**;297 次之(penalty 60,闭合);217/228 也折了当对照。→ 切点不是拍脑袋,是"最小结构破坏 + 能闭合"选出来的,且 cp233 是三预测器 5/5 的 definitive lead。
### linker 长度:两个 linker 各 3 档都建过
- int(NC linker):**10 / 15 / 20** aa 都建了(monomers/)。
- inter(subunit linker):**10 / 15 / 20** aa 都建了(pentamers/)。
- cp233 上 **int15_inter10 胜出**(三预测器 5/5;MD 稳定;结构叠合 1.8 Å)。

### 还需不需要再优化?
- **cp233:不需要再大规模筛——切点(20 选 1)+ 两个 linker(各 3 档)都扫过,已是三预测器 lead,现在够了。** 若要精修可在 int/inter ±2-3 aa 做细网格,但不是瓶颈。
- **cp285 / cp297:验证较浅**——主筛只把 217/228/233 折成了 pentamer;285/297 是备选切点,linker 组合没系统扫。→ 这正是 `tiled_msa_fold` 要补的(见 §4):用 tiled MSA 公平折 cp285/297,决定要不要升级为第二/第三位点。

## 3. 可信度证据栈(结构在最前,详见 `outputs/structural_validation/STRUCTURAL_CREDIBILITY.md`)
1. **结构叠合**:亚基 vs 天然 apo **RMSD 1.80 Å / TM 0.94**;vs 实验 7JQQ **3.72 Å / TM 0.79**;整环 44 vs 42 Å、5/5 归位。"看着 off" = 单链化 linker,不是坏折叠。
2. **三预测器**:Boltz+OF3+AF3 全 5/5 M2(tiled MSA)。
3. **MD 物理稳定**:cp233 像 apo 一样保持闭合(`md/openmm_validation/`)。
4. **独立干实验(不共享预测器先验)**:
   - ✅ ATP 口袋几何(CP 后五亚基一致、≤0.7 Å)
   - ✅ ENM 开/闭(没被 linker 焊死)
   - ✅ ESM-2 naturalness(linker 序列自然)
   - ✅ MD(见上)
   - ❌ **Rosetta / MM-GBSA 界面能** —— 唯一还没跑的正交信号(最后一条腿)。需 PyRosetta 或 OpenMM-MMGBSA;跑完 5 条正交证据齐活。

## 4. `tiled_msa_fold` 是什么(为什么需要它)
**问题**:单序列(single-seq)结构预测**无法评估这类大单链环**——连已验证的 cp233 在 single-seq 下都掉到 M2 0/5(tiled MSA 下是 5/5)。因为环里 5 份是同一个域的拷贝,single-seq 给不出足够的共进化信号。
**tiled_msa_fold = 一个可复用管线**:(1) 对**单个亚基域**跑 MSA(homolog 搜索),(2) 把它**沿对角线平铺 N 份**成 block-diagonal MSA(每份拷贝各自有信号、互不串扰),(3) 送 **免费 Boltz-2 NIM** 折叠,(4) `reproduce/score_m2.py` 打分(M1 sequential + M2 trans-finger)。
**用途**:公平折所有生成式大环——cp285/297、gp16 RFdiffusion 铺环、Rho 三方法、E119Q 死座系列——得到论文核心的"蛋白 × 方法 × M2"对比表(把 single-seq 的 confounded 结果替换掉)。建好后可无限复用。
