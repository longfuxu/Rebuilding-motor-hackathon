# CS prompt — existing AF3 structures + GCP compute for AF3/analysis
_粘贴给 CS。两件事:(1) 现有 AF3 结构在哪、直接读别重折;(2) 需要跑 AF3 或重计算时用 GCP。_

---

**你们不需要重折已经有的结构 —— 现有 AF3 结果已分析完,直接读取:**
- **主目录**:`gp16_design/outputs/af3_2026_07_09/`
  - `AF3_INDEX.md` —— 所有 AF3 fold 的一张总表(M2 判据 + verdict)。
  - `AF3_SCORING_REPORT.md` —— 9 个生成式/CP fold 的详细打分(cp233_novel d2–d8、cp285、cp297、rfdiff d0/d1/d2),对抗性验证过。
  - 每个 `fold_2026_07_09_*/` 里是 5 个 AF3 model 的 `.cif` + `*summary_confidences*.json`(full_data 大文件没存)。
  - 已覆盖:**cp233_WT(apo + ATP·Mg)、E119Q(1seat/5seat,apo + ATP·Mg)、cp233_novel d2/d3/d5/d7/d8、cp285、cp297、rfdiff_ring d0/d1/d2**。
- **E119Q 功能分析**:`gp16_design/outputs/e119q_functional/REPORT.md`(死座 3 预测器 + ligand 都 5/5 装配容忍;WT E→Mg 4.2Å = E119Q Q→Mg 4.2Å → 结构上沉默的催化敲除)。
- **cp233 主 lead 的 AF3(三预测器 5/5)**:`gp16_design/cs_sessions/session1_main_gp16/handoff/cp233_af3/`。
- **打分复用**:`gp16_design/reproduce/score_m2.py`(M1/M2)或 `pipelines/tiled_msa_fold/`(handedness-robust)。cp233 landmarks:copies `A:1-342,A:353-694,...,A:1409-1750`,r146_incopy=255,walker=133-140,E119=228+352k。rfdiff:copies of 377,r146=143,walker=21-28。**先查 AF3_INDEX,已有的别重折。**

**需要跑新的 AF3、或做重计算(MD/MM-GBSA/大批折叠)时,可以用 GCP GPU(有 credit,不产生自费):**
- **工具**:`gp16_design/gcp_pipeline/gcp_gpu.py`(一条命令:开 Spot GPU VM → 装 → 跑 → 下结果 → **删 VM**)。gcloud 全路径 `/Users/longfu/google-cloud-sdk/bin/gcloud`,project `longfu-protein-gpu`。
- **★ 计费必须用 credit 账户** `0161E2-E9BA4A-8F5A6A`(Google 开发者额度,不产生自费)。**成本纪律(硬)**:只用 **Spot** GPU + **跑完立刻 delete VM**(idle GPU 烧额度);跑完 `gcloud compute instances list` 确认 0 实例。细节见 system memory `gcp-gpu-account`。
- **跑 AF3 本身的诚实说明**:AF3(AlphaFold3)推理需要 **DeepMind 授权的模型权重**(github.com/google-deepmind/alphafold3,学术用途申请)+ GPU;拿到权重后可在 GCP GPU 上跑,给 `gcp_gpu.py` 加一个 af3 recipe(clone alphafold3 + 权重 + `run_alphafold.py`)。**若暂无 AF3 权重**:用**免费 Boltz-2 / OpenFold3 NIM**(`health.api.nvidia.com/v1/biology/mit/boltz2/predict`,key 在主 checkout `.env`)当第 2、3 个预测器——无需 GPU、无需权重,已是我们跨预测器的主力;AF3 网页版(alphafoldserver.com)留给少量关键构建。
- **重计算适合上 GCP 的**:显式溶剂 MD、MM-GBSA 全套、RFdiffusion 生成式、大批 tiled-MSA 折叠。轻量的(几个折叠、ENM/PRS、打分)用免费 NIM + 本地就够,别开 VM。

**一句话**:先读 `AF3_INDEX.md` 复用已有 AF3 结构;要新 AF3 或重算就用 GCP(Spot + 用完删 + credit 计费),没有 AF3 权重就用免费 Boltz-2/OF3 NIM 顶替。
