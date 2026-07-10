#!/usr/bin/env python3
"""Replace the MM-GBSA table placeholder/running-rows in round2/REPORT.md with the finished
numbers from scores/mmgbsa/mmgbsa_summary.json, and append a one-line interpretation."""
import os, json, re
HERE = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.join(HERE, "REPORT.md")
SUMM = os.path.join(HERE, "scores", "mmgbsa", "mmgbsa_summary.json")
LABEL = {"cp233_WT": "cp233_WT", "de_S99W": "de_S99W", "de_M307W_N128F": "de_M307W_N128F",
         "de_M307F_S127F": "de_M307F_S127F", "de_M307W": "de_M307W(R1 基准)"}
ORDER = ["cp233_WT", "de_S99W", "de_M307W_N128F", "de_M307F_S127F", "de_M307W"]

d = json.load(open(SUMM))
rows = {r["name"]: r for r in d["rows"]}
lines = ["", "| 系统 | mean dG_bind (kcal/mol) | Δ vs WT |", "|---|---|---|"]
for n in ORDER:
    r = rows.get(n, {})
    dg = r.get("mean_dG_bind"); dd = r.get("delta_vs_WT")
    dgs = f"{dg:.1f}" if dg is not None else "n/a"
    dds = ("0" if n == "cp233_WT" else (f"{dd:+.1f}" if dd is not None else "n/a"))
    lines.append(f"| {LABEL[n]} | {dgs} | {dds} |")
# interpretation
wt = rows.get("cp233_WT", {}).get("mean_dG_bind")
s99 = rows.get("de_S99W", {}).get("delta_vs_WT")
m307 = rows.get("de_M307W", {}).get("delta_vs_WT")
interp = "\n> 判读:MM-GBSA(独立于 ENM/DL 的第三方界面能读数)"
if s99 is not None and m307 is not None:
    tighter_s99 = "更紧" if s99 < 0 else "更松"
    interp += f"。de_S99W 界面能 vs WT = {s99:+.1f} kcal/mol({tighter_s99});M307W = {m307:+.1f}。"
    if s99 < m307:
        interp += "S99W 的界面比 M307W 还紧,与其 force-proxy 打平/略优一致。"
    elif abs(s99 - m307) < 5:
        interp += "S99W 与 M307W 的界面紧度在 MM-GBSA 上基本相当,佐证二者打平。"
    else:
        interp += "M307W 的界面更紧。"
table = "\n".join(lines) + "\n" + interp + "\n"

s = open(REPORT).read()
# replace everything between the MM-GBSA "running" table and the next "## " heading
pat = re.compile(r"(> more negative dG_bind = 更紧的亚基间耦合。)(.*?)(\n## 对北极星)", re.S)
s2 = pat.sub(lambda m: m.group(1) + table + m.group(3), s)
open(REPORT, "w").write(s2)
print("REPORT.md MM-GBSA table filled.")
print(table)
