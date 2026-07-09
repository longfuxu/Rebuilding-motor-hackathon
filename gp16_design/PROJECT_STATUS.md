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
- **MD(物理独立验证,Task 4)——正在 L4 上跑**(A100 缺货)。OpenMM 隐式溶剂 (OBC2, OpenCL 平台),对比 apo / 7JQQ / cp233 设计的 Cα-RMSD + 环呼吸半径随时间。脚本 `outputs/rho_generative/../` 见 §7。**⚠ 要产出动画见 §6。**
- **Rho 三方法 fold(single-seq,confounded)**:直连 0/6、CP330 0/6、CP390 0/6 —— 全被 single-seq 压平(见 §4),**不能区分方法**。设计已建好(`outputs/rho_generative/`)。
- **独立干实验(非预测器)已做**:ATP 位点几何、ENM 开闭、ESM naturalness——三个正交信号都支持 cp233(`outputs/INDEPENDENT_VALIDATION.md`)。

## 6. ⚠ 待办:MD 动画(用户要求,给未来 agent)
**要求:MD 跑完后,把三条轨迹(A_apo / B_7jqq / C_cp233-design)各做成一个动画/movie,用于 tracing + 报告。**
- 需要的文件:每个体系的轨迹 `mdout/<X>/<X>.dcd` + 起始拓扑 `mdout/<X>/<X>_start.pdb`。
- ⚠ 当前 MD orchestrator(`md_orchestrator2.py`)默认只抓分析 JSON,**没抓 DCD**——要产动画必须在会话 stop 前把 DCD 也下下来(改 capture 那行的 tar 加 `mdout/*/*.dcd`,或单独 `colab download`)。
- 做法:下载 DCD+start.pdb 后,用 nglview / VMD / PyMOL / mdtraj 渲染;或转 GIF/MP4。对比动画:设计环是否像天然一样稳定,还是塌陷/开环。

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
