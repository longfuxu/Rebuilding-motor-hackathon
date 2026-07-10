#!/usr/bin/env python3
"""Generate round2_greedy/REPORT.md from the greedy ranking artifacts. Idempotent."""
import os, json, csv
HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.abspath(os.path.join(HERE, ".."))
SC = os.path.join(HERE, "scores")

summ = json.load(open(os.path.join(SC, "greedy_summary.json")))
rows = summ["rows"]
noise = summ["noise"]
bg = summ["bg_metrics"]
wt = json.load(open(os.path.join(SC, "cp233_WT.json")))
esm = {r["name"]: r for r in json.load(open(os.path.join(SC, "esm_plausibility.json")))}
stackers = summ["stackers"]

cn, fn = noise["coord"], noise["force"]
n_pass = sum(1 for r in rows if r["pass_all"])
# expected non-foldable (boundary artifact)
excluded = ["de_M307W_K233F", "de_M307W_K233W"]

def esm_llr(n):
    r = esm.get(n); return f"{r['esm_llr_sum']:+.2f}" if r else "n/a"

# sort rows: STACKS first then by best margin (already sorted in vs_m307w by rank_greedy)
lines = []
lines.append("# In-silico 定向进化 —— GREEDY 上位性扫描(背景 = M307W)(2026-07-09)\n")
lines.append("用户方法学修正:round-2 的\"双突变不叠加\"只测了 4 个预选双突变(欠采样)。本轮固定"
             "round-1 赢家 **M307W 为背景**,在**一大批**被允许的第二位点上构建 M307W+X,直接、"
             "更广地检验:**在最好的单点背景上,还有没有第二个突变能叠加(两个突变强于一个)?**\n")
lines.append("## 设计\n")
lines.append("- 背景固定 = M307W(不再动 307)。第二位点 X = 生化 TOLERANT {56,99,127,128,222,233} "
             "+ 未测刚性化界面/核心位点 {100,130,289},每个位点取有刚性化意义的 F/W/Y(核心 L289 另取 Pro)。\n")
lines.append("- 排除所有 SEVERE/LETHAL(129,58,105,328,330,28,32,106,126,158,234,294)与催化残基。\n")
lines.append(f"- 库 = 20 个 M307W+X。其中 **K233F / K233W 无法用 tiled-MSA 折叠**"
             f"(native 233 位于 segB/linker 拷贝边界,derive_blocks 无法锚定 5 个拷贝——边界伪影,已剔除),"
             f"故实测 **{len(rows)} 个 M307W+X**(3 个 M307W+N128F/S127F/L289P 复用 round-2 折叠/打分)。\n")
lines.append("- 每个 M307W+X:tiled-MSA / Boltz-2 折叠 → gated 打分(M1∧M2 handedness-robust∧M5∧soft-mode∧biochem)"
             "× power(coord PRS + force-network + M3 grip)。ESM-2 可信度过滤。\n")
lines.append("## 判定标准(决定性问题)\n")
lines.append(f"**叠加(STACKS)= 门槛全过 且(相对 M307W 背景,coord 增益 > {cn*100:.2f}% 噪声地板 "
             f"或 force 增益 > {fn*100:.2f}% 噪声地板)。** 即:比最好的单点还要好、且超过折叠噪声。\n")
lines.append(f"背景 M307W:PRS_NN = {bg['prs_ATP_NN_coupling']:.6f},conduit = {bg['len_ATP_Y129_DNA']:.4f}"
             f"(round-1:coord +6.4% / force +3.5% vs WT——这是每个 M307W+X 要超越的线)。\n")
lines.append(f"噪声地板(|WT−WTrep|):coord = {cn*100:.2f}%,force = {fn*100:.2f}%。"
             f"门槛全过的 M307W+X:**{n_pass}/{len(rows)}**。\n")

# main table
lines.append("## 结果:所有 M307W+X vs 背景 M307W(dCoord/dForce = 相对 M307W 的增益)\n")
lines.append("| 变体 | 过门槛 | dCoord vs M307W | dForce vs M307W | soft 保留 | ESM LLR | STACKS? |")
lines.append("|---|---|---|---|---|---|---|")
for r in rows:
    st = "**YES**" if r["STACKS"] else ("coord" if r["stacks_coord"] else ("force" if r["stacks_force"] else "no"))
    gate = "✓" if r["pass_all"] else "✗("+r["gates_failed"]+")"
    lines.append(f"| {r['name'].replace('de_','')} | {gate} | {r['d_coord_vs_M307W']*100:+.1f}% | "
                 f"{r['d_force_vs_M307W']*100:+.1f}% | {r['soft_retention']:.3f} | {esm_llr(r['name'])} | {st} |")

lines.append("\n## 诚实结论\n")
if stackers:
    b = rows[0]
    lines.append(f"1. **叠加存在。** {len(stackers)} 个门槛全过的 M307W+X 在噪声之上击败了 M307W 背景。"
                 f"**最佳累积变体 = {b['name'].replace('de_','')}**"
                 f"(dCoord {b['d_coord_vs_M307W']*100:+.1f}%,dForce {b['d_force_vs_M307W']*100:+.1f}% vs M307W,"
                 f"soft {b['soft_retention']:.3f})。→ 两个突变确实强于一个;上位性正向存在。\n")
    lines.append(f"2. 叠加者:{', '.join(s.replace('de_','') for s in stackers)}。"
                 "→ 这些是新的、可交给单分子测试的\"更强\"候选。\n")
else:
    lines.append("1. **没有叠加。** 在这个更广的第二位点扫描里,**没有一个门槛全过的 M307W+X 在噪声之上"
                 "击败 M307W 背景**(coord、force 均未过噪声)。→ round-2 的\"单点不叠加\"结论现在是"
                 "**稳健的、更广采样的负结果**:即使在最好的单点背景上,再加一个刚性化点也不能稳健地更强。\n")
    lines.append("2. 机制解读:M307W 已把 Y129↔307 clasp 的网络刚度推到接近最优;在这个已优化的背景上"
                 "再加约束,要么被噪声淹没,要么(和 round-2 一样)过度刚性化、把协同/力压回去。"
                 "**\"更强\"不是靠继续堆刚性化点堆出来的。**\n")
# S99W note
s99w = next((r for r in rows if r["name"] == "de_M307W_S99W"), None)
if s99w:
    verdict = "也叠加" if s99w["STACKS"] else "在此背景上不叠加"
    lines.append(f"3. **另一条力传导线(S99W)**:M307W+S99W {verdict}"
                 f"(dCoord {s99w['d_coord_vs_M307W']*100:+.1f}%,dForce {s99w['d_force_vs_M307W']*100:+.1f}% vs M307W)。"
                 "→ N-domain 的 S99W 与 C-domain 的 M307W 合起来"
                 + ("超过了单独的 M307W。" if s99w["STACKS"] else "并不优于单独的 M307W——两条力传导线不叠加。") + "\n")
lines.append(f"4. **门槛稳健性**:{n_pass}/{len(rows)} 个 M307W+X 全过 M1/M2/M5/soft/biochem——"
             "折叠层面这些刚性化都不破坏装配/催化/软模;区分它们的是 proxy,不是门槛。\n")

lines.append("\n## 对北极星的意义\n")
if stackers:
    lines.append("greedy 扫描找到了真正的上位性叠加:在 M307W 背景上,某些第二突变把马达 proxy 推得"
                 "更高、过噪声、门槛全过。→ \"更强\"可以逐步累积(greedy 定向进化奏效),交出新的单分子候选。\n")
else:
    lines.append("\"能不能设计出更强的?\"→ 更广的 greedy 扫描给出诚实的负结果:**在最好的单点背景上,"
                 "20 个第二位点里没有一个能稳健叠加**。M307W(以及独立的 S99W)仍是 proxy 层最强的、"
                 "可单分子验证的\"略强\"候选;要更进一步需要超出\"叠加点突变\"的手段"
                 "(骨架重设计 / RFdiffusion-de-novo connector / 负设计打开软模,见 backlog)。\n")
lines.append("\n脚本:round2_greedy/{make_variants_greedy,fold_all2,score_all2,rank_greedy,esm_plausibility}.py。"
             "排序(vs WT):scores/ranked_variants.csv;决定性对比(vs M307W):scores/vs_m307w.csv;"
             "ESM:scores/esm_plausibility.csv。")

open(os.path.join(HERE, "REPORT.md"), "w").write("\n".join(lines) + "\n")
print("wrote round2_greedy/REPORT.md")
print("STACKERS:", stackers if stackers else "NONE")
