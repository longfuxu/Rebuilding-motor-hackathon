# PROJECT STATUS — 单链环 ATPase 设计（READ THIS FIRST）

> **如果只读一个文件,读这个。** 这是全项目的当前状态 + 完整任务图 + 进度 + 下一步。
> 细节见文末"详情文档"。完整时间顺序日志见 `../HANDOFF.md`(很长,是流水账)。
> **工作守则:每完成一步就写进度 + commit**,任何 agent 都能接着做。判据永远用 sequential M2(`reproduce/score_m2.py`)+ M1,**绝不看 global pTM**;apo;交叉预测器=filter 不是 validation;不写"validated"。

_最后更新: 2026-07-09_

## 1. 项目一句话
造一条**可搭建、单条链、带一个基因可寻址位点**的环形 ATPase 分子马达(像 single-chain ClpX),并证明一套**可迁移的框架**:fold-in-context + interface-resolved M2 coupler 指标 + **升级阶梯(N-C 直连 → 环状置换 CP → RFdiffusion 生成式)**——按序试到能闭合为止。

## 2. 当前 leads(gp16)
- **cp233_int15_inter10(单链,definitive lead)**:环状置换,**三预测器(Boltz+OF3+AF3)全 5/5**(tiled MSA)。ATP 位点几何 CP 后保住(定位对照 ✓)。
- **B1_L40_E119Q**:多链(共价二聚体+3 WT),可寻址死座 lead。

## 3. 核心科学:三方法 × 三蛋白矩阵("同框架、不同解")
判据两条几何量:**末端离功能通道多远(堵不堵)** + **直连 C→N gap 多长**。
| 蛋白 | 末端 vs 通道 | 直连 gap | 最优拓扑 | fold 证据 |
|---|---|---|---|---|
| **ClpX**(六聚体) | 62Å 孔外 | ~20Å | **直连** | 回顾性验证过(good 6/6, R307A/E lesion 0/6);CP 更差 |
| **gp16**(五聚体) | C端 6.2Å **扎 DNA 通道** | 53–66Å | **CP** | cp233 三预测器 5/5;N-C 直连 AF3 0/5 scrambled;RFdiffusion 铺环 register 对但 M2 待 AF3 |
| **Rho**(六聚体) | 40–60Å 孔外 | 41Å(motor) | **测试中** | 直连+CP 都测了(见 §5 confound) |

## 4. ★ 关键方法学发现(写进论文的 limitation)
**单序列(single-seq)结构预测无法评估这类大单链环**——连已验证的 gp16 cp233 在 single-seq Boltz 下都掉到 **M2 0/5**(tiled MSA 下是 5/5)。所以任何用 single-seq 折的大环 M2 都不可信,**必须用 tiled MSA**。single-seq 只能当"能不能折成环"的粗筛。

## 5. 正在进行 / 刚完成
- **✅ MD(物理独立验证,Task 4)——已完成(2026-07-09)**。OpenMM 隐式溶剂 (GBSA-OBC2, OpenCL 平台),对比 apo / 7JQQ-helical / cp233 设计。**结论:cp233 设计像天然 apo 一样保持闭合+对称(radius_CV 0.009、界面接触 584→657 反而收紧),不开环不塌陷——明显不同于 7JQQ 上 ATP 的螺旋态(卡 3/5、seam 持续)→ 设计通过了独立于结构预测器的物理稳定性检验**。产出:`md/openmm_validation/`(RESULTS.md/METHODS.md/figs/CSVs/code)。Caveats:隐式溶剂 + 短单轨迹(1–3 ns,A100 会 ~20–30min 自终止;判据在共享 0–1 ns 窗内已解析)。显式 TIP3P 需 CUDA OpenMM,列为 Phase-2。
- **Rho 生成式 RFdiffusion——Claude Code 正在 Colab 跑**(motif res175-414,连接件 30-52,跨 46Å gap);出 backbone 后交 tiled-MSA 折叠(§7.2)。
- **Rho 三方法 fold(single-seq,confounded)**:直连 0/6、CP330 0/6、CP390 0/6 —— 全被 single-seq 压平(见 §4),**不能区分方法**。设计已建好(`outputs/rho_generative/`)。
- **CS 三 session 已镜像进 repo**:`cs_sessions/`(reports/leads/scores;341M 原始 fold 模型 git-ignore,存本地)。下一步 CS 任务见 `NEXT_CS_TASKS_v2.md`。
- **独立干实验(非预测器)已做**:ATP 位点几何、ENM 开闭、ESM naturalness——三个正交信号都支持 cp233(`outputs/INDEPENDENT_VALIDATION.md`)。

## 6. ⚠ 待办:MD 动画(用户要求,给未来 agent)
**要求:把三条轨迹(A_apo / B_7jqq / C_cp233-design)各做成动画/movie,用于 tracing + 报告。**
- **⚠ 现状(2026-07-09):MD 已跑完并得出物理结论,但 DCD 全轨迹没保存**——`md/openmm_validation/results/<X>/` 只有 `<X>_start.pdb` + `<X>_final.pdb` + timeseries CSV,没有 `.dcd`。所以**真正的轨迹动画还产不出来**。
- 两条路:(a) **start→final 形变 morph**(用现有 start/final PDB,PyMOL `morph`,便宜、无 GPU,能直观展示"闭合→仍闭合"稳定);(b) **重跑一遍 MD 并保存 DCD**(在 `md/openmm_validation/code/md_driver.py` 里加 `DCDReporter`,每 ~10ps 存一帧),再 nglview/VMD/PyMOL 渲染 GIF/MP4。物理判据已定,动画属"报告可视化",优先级次于 §7。

## 7. 下一步(open)
1. **AF3**(用户跑):`AF3_SEQUENCES_TODO.txt` 里 11 条(5 条 cp233_NOVEL de-novo + 4 条 RFdiffusion 铺环 + cp285/cp297)。跑完交回 score_m2。
2. **Rho 公平对比(需 tiled MSA)**:把 Rho 直连/CP/**生成式**三构建用 cofactor.fold + Rho 专属 tiled MSA 折 + score(single-seq 不行)。生成式:motif 已建 `rho_motif_AB.pdb`,跑 RFdiffusion de-novo connector(需 GPU)。
3. **RFdiffusion 铺环 gp16**:tiled-MSA 重折(single-seq 是 M2 0/5,不公平)。
4. MD 动画(§6)。

## 8. 详情文档(按需深入)
- `DESIGN_METHODS_COMPARISON_AND_NEXT.md` — 三方法对比 + 各任务 prompt/工具/session。
- `outputs/RFDIFFUSION_OVERNIGHT_RESULTS.md` — RFdiffusion 全过程(Mode B 连接件→MPNN→铺环)。
- `outputs/INDEPENDENT_VALIDATION.md` — 三个正交干实验。
- `outputs/rho_generative/RHO_METHODS_COMPARISON.md` — Rho 三方法 + single-seq confound。
- `cs_session_outputs/clpx_generalization/` — ClpX 回顾性验证 + clpx_vs_gp16 拓扑。
- `NEXT_CS_TASKS.md` — 给 CS 的任务。`AF3_SEQUENCES_TODO.txt` — AF3 队列。
- **CS 输出在别的本地工作区**(不同步 repo),见 memory `cs-workspace-locations`;主动 grep `/Users/longfu/.claude-science/orgs/*/workspaces/*/`。
- 运维坑:memory `colab-cli-overnight-gotchas`(exec/download hang→本地 python 硬超时;STOP 会话省额度;MD 用 OpenCL 平台非 CUDA)。

## 9. 运维现状

- **GCP GPU(可靠替代 Colab,已配好 2026-07-09):** 项目 `longfu-protein-gpu`,billing=`0161E2-E9BA4A-8F5A6A`(Google Gemeni 信用额度账户,~$45/月,**必须用它,不产生自费**),$40 预算告警。GPU 配额已批:GPUS_ALL_REGIONS=1、L4/T4 us-central1=1。gcloud 在 `~/google-cloud-sdk/bin/gcloud`(全路径),已登录。**守则:Spot GPU + 用完立刻 delete VM**(idle 烧额度)。细节见 system memory `gcp-gpu-account`。下一步 GPU 活(Rho RFdiffusion / tiled-MSA / MD)可搬到这里,比 Colab 稳。
- Colab 会话:MD 在 `mdrun2`(L4)。**跑完必 stop**(idle A100 ~11–13 units/h)。A100 目前缺货,MD/RFdiffusion 用 L4/T4 回退。
- 无 git remote(本地 commit)。commit 作者=longfuxu,无 AI 署名(硬性)。
