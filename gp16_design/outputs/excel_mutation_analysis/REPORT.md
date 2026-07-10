# phi29 生化突变数据 × 我们的判据 —— 校准 & 为什么 CP233(2026-07-09)

数据:`phi29_mapped_residues.xlsx`(92 残基,bulk packaging/ATPase/slippage 表型,非单分子)。清洗后 49 条有信息 → `cleaned_residue_phenotype.csv`。脚本 `analyze_mutations.py`。

## Q1 —— proxy 验证(生化关键残基 vs 我们 proxy-flagged)
- **预测重要残基 21 个,52%(11/21)被我们的 proxy 命中**(DNA/pRNA 接触、Y129、R146、Walker)。
- **★ 强验证:Y129(我们算的 4.6× 力传导枢纽)生化上 Y129A = SEVERE(近致死)** —— 计算预测的"力枢纽"被 wet 独立证实为功能必需。
- **Tyr297 flagged**(= 我们 M3 说 cp297 切坏的 DNA 接触)、**R146 flagged**(反式手指)。
- **DNA 接触残基分"承重 vs 冗余"**:承重(N126A LETHAL、Y129A SEVERE、R330A SEVERE);冗余(K56/S99/S127/N128 都 TOLERANT)。→ 不是所有接触残基等价,M3 应加权。

## ★ Proxy 缺口(生化关键但我们 proxy 漏掉)—— 诚实,要补
主要是**催化活性位点(cis/trans)残基**:Tyr32、E58(glu switch)、Lys105/Ser106(trans active site/first catalysis)、Asn158、以及 grip 的 Arg234/Lys294/Lys328。
→ **M5(ATP 口袋)判据应扩到完整催化集**(不止 Walker-A/B + R146,还要 glu-switch E58、trans 催化 K105/S106、N158)。这是论文的一个 honest fix。
(注:functional_contacts JSON 把 DNA 接触写 55-60 但漏了 58;E58A 实测 SEVERE = glu-switch,范围 vs 枚举集不一致,需对齐。)

## Q2 —— 为什么 CP233 更好(生化支撑)
| CP 位点 | 切点残基类别 | 在接触面? | 最近关键残基 |
|---|---|---|---|
| **CP233** | **TOLERANT** | **否** | 234(相邻,|d|=1)—— 唯一 caveat |
| CP285 | UNTESTED | 否 | 294(|d|=9) |
| **CP297** | PREDICTED-ONLY | **是(在接触面)** | 294(|d|=3) |
→ **cp233 切在耐受、非接触面的位点;cp297 切在 DNA 接触面、靠近关键残基** —— 生化独立支持 cp233 > cp297(和 M3"cp297 通道太窄+切 DNA 接触"一致)。唯一小 caveat:233 紧邻关键 Arg234。

## Q3 —— 定向进化的允许突变先验(★ 关键校准)
- **TOLERANT(可安全突变)**:Phe6*、R53*、Lys56、Ser99、Ser127、**Asn128**、Gln222、Lys233(*=仅保守替换)。
- **SEVERE(近致死,别碰)**:E58、Lys105、**Tyr129**、Lys328、Arg330。
- **LETHAL(绝对别碰)**:Ile28、Tyr32、Ser106、Asn126、Asn158、Arg234、Lys294。
- **★★ 对定向进化的直接后果:我们的 rigidify 靶点里 Y129 是近致死(Y129A SEVERE)——不能突变 Y129 本身;应 rigidify 它周围(N128 是 TOLERANT)。** 即:加固力传导路径要动 128、130 一带,不要动 129。

## 一句话
生化数据(a)**独立验证** Y129 力枢纽 + Tyr297/R146;(b)**暴露 proxy 缺口**(催化活性位点 → 补 M5);(c)**支撑"CP233 更好"**;(d)给定向进化一份**tolerant/lethal 先验**(关键:别动 Y129,动它旁边)。
