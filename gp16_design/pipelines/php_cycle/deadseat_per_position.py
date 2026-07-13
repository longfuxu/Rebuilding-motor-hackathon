#!/usr/bin/env python3
"""PER-POSITION dead-seat scan (follow-on program Aim 1b): fold cp233 single chain with ONE R146A dead trans-finger
at EACH of the 5 ring positions (+WT), on the free Boltz-2 NIM, 3 replicates, handedness-robust M2
scoring. Answers the mixed-ring design question: does a single dead seat's coordination effect depend
on WHICH seat (is a position 'special'?), or is it position-independent (as the count-gradient implied)?

Output: outputs/php_cycle/special_subunit/{manifests/, <name>.result.json, PER_POSITION_RESULTS.md}
Run: NVIDIA_API_KEY read from main-checkout .env. python3 deadseat_per_position.py
"""
import importlib.util, os, json, threading, re

REPO = "/Users/longfu/Developer/claude-science-hackthon/.claude/worktrees/export-gp16-results/gp16_design"
TDIR = f"{REPO}/pipelines/tiled_msa_fold"
spec = importlib.util.spec_from_file_location("t", f"{TDIR}/tiled_msa_fold.py")
t = importlib.util.module_from_spec(spec); spec.loader.exec_module(t)

OUT = f"{REPO}/outputs/php_cycle/special_subunit"; MDIR = f"{OUT}/manifests"
os.makedirs(MDIR, exist_ok=True)
KEY = [l.split("=", 1)[1].strip().strip('"') for l in
       open("/Users/longfu/Developer/claude-science-hackthon/.env") if l.startswith("NVIDIA_API_KEY")][0]
NREP = 3
lock = threading.Lock()


def build_manifests():
    wt = t.load_manifest(f"{TDIR}/manifests/cp233_WT.json")
    _a3m, meta = t.do_tile(wt); seq = wt["sequence"]
    starts = [int(tok.split(":")[1].split("-")[0]) for tok in meta["copies"].split(",")]
    r146 = meta["r146_incopy"]
    positions = [lo + r146 - 1 for lo in starts]      # 1-based chain residue of R146 in each copy
    assert all(seq[p - 1] == "R" for p in positions), "R146 anchor mismatch"
    def mut(s, p):
        l = list(s); l[p - 1] = "A"; return "".join(l)
    names = ["WT"] + [f"R146A_pos{i+1}" for i in range(5)]
    seqs = [seq] + [mut(seq, positions[i]) for i in range(5)]
    for nm, sq in zip(names, seqs):
        json.dump({"name": nm, "sequence": sq, "n_copies": 5, "monomer_a3m": wt["monomer_a3m"]},
                  open(f"{MDIR}/{nm}.json", "w"))
    print(f"built {len(names)} manifests; R146 at chain residues {positions}", flush=True)
    return names


def run_one(name):
    m = t.load_manifest(f"{MDIR}/{name}.json")
    a3m, meta = t.do_tile(m)
    resp, dt = t.fold_nim(m["sequence"], a3m, KEY, samples=NREP)
    if "__err" in resp:
        with lock: print(f"[{name}] FOLD ERR {resp['__err'][:120]}", flush=True)
        json.dump({"name": name, "ok": False, "err": resp["__err"]}, open(f"{OUT}/{name}.result.json", "w"), indent=2)
        return
    reps = []
    for i, st in enumerate(resp.get("structures", [])):
        cif = f"{OUT}/{name}_rep{i}.cif"; open(cif, "w").write(st["structure"])
        v, out = t.score_ring(cif, meta["copies"], meta["r146_incopy"], meta["walker_incopy"])
        m1 = next((l for l in out.splitlines() if l.startswith("# M1:")), "")
        cv = (re.search(r"radius_CV ([\d.]+)", m1) or [0, ""])[1]
        reps.append({"rep": i, "engaged": v["engaged"], "handedness": v["handedness"],
                     "m2_forward": v.get("m2_forward"), "m2_reverse": v.get("m2_reverse"),
                     "radius_CV": cv, "sequential": "sequential_consistent: YES" in m1, "M1": m1})
    eng = [r["engaged"] for r in reps]
    best = max(reps, key=lambda r: (r["engaged"], -float(r["radius_CV"] or 9)))
    if os.path.exists(f"{OUT}/{name}_rep{best['rep']}.cif"):
        os.replace(f"{OUT}/{name}_rep{best['rep']}.cif", f"{OUT}/{name}.cif")
    for r in reps:
        p = f"{OUT}/{name}_rep{r['rep']}.cif"
        if os.path.exists(p): os.remove(p)
    res = {"name": name, "ok": True, "wall_s": dt, "engaged_values": eng, "engaged_max": max(eng),
           "engaged_min": min(eng), "engaged_mean": round(sum(eng) / len(eng), 2),
           "best_rep": best["rep"], "best_m2_forward": best["m2_forward"], "best_m2_reverse": best["m2_reverse"],
           "best_radius_CV": best["radius_CV"], "best_sequential": best["sequential"]}
    json.dump(res, open(f"{OUT}/{name}.result.json", "w"), indent=2)
    with lock: print(f"[{name}] engaged {eng} (max {max(eng)}) CV {best['radius_CV']} seq {best['sequential']} {dt}s", flush=True)


def summarize(names):
    rows = []
    for nm in names:
        p = f"{OUT}/{nm}.result.json"
        if os.path.exists(p): rows.append(json.load(open(p)))
    L = ["# Per-position dead-seat scan — cp233 single chain, R146A at each of the 5 seats",
         "", "One dead trans-finger (R146A) placed at each ring position, folded on Boltz-2 NIM (3 reps),",
         "M2 = trans-R146→neighbour Walker-A <8 Å, handedness-robust (engaged_max of 3 reps). WT baseline = 5.",
         "", "| construct | engaged_max | engaged (3 reps) | radius_CV | sequential | interpretation |",
         "|---|---|---|---|---|---|"]
    wt = next((r for r in rows if r["name"] == "WT" and r.get("ok")), None)
    for r in rows:
        if not r.get("ok"):
            L.append(f"| {r['name']} | FOLD ERR | | | | {r.get('err','')[:40]} |"); continue
        interp = "baseline (all fingers)" if r["name"] == "WT" else \
                 ("one finger lost (local, 5→4 as expected)" if r["engaged_max"] == 4 else
                  ("cooperative extra loss" if r["engaged_max"] < 4 else "no coordination loss (robust seat)"))
        L.append(f"| {r['name']} | **{r['engaged_max']}** | {r['engaged_values']} | {r['best_radius_CV']} | "
                 f"{r['best_sequential']} | {interp} |")
    engs = {r["name"]: r["engaged_max"] for r in rows if r.get("ok") and r["name"] != "WT"}
    if engs:
        vals = list(engs.values())
        spread = max(vals) - min(vals)
        L += ["", f"**Position dependence:** engaged_max across the 5 positions = {engs}. Spread = **{spread}**.",
              "- spread 0 ⇒ position-INDEPENDENT (every seat loses exactly one finger, 5→4) — no 'special' seat by this readout.",
              "- spread ≥1 ⇒ position-DEPENDENT — some seat's dead finger costs more/less coordination (a candidate special seat).",
              "", "_Caveat: Boltz-2 static fold + M2 geometric proxy; predictor basin-bias (all cp233 fold planar);",
              "3 reps only. This is a design-prioritisation prediction for the mixed-ring experiment, not proof._"]
    open(f"{OUT}/PER_POSITION_RESULTS.md", "w").write("\n".join(L))
    print("\n".join(L[4:]))


def main():
    names = build_manifests()
    sem = threading.Semaphore(3); ths = []
    def w(n):
        with sem:
            try: run_one(n)
            except Exception as e:
                with lock: print(f"[{n}] EXC {repr(e)[:140]}", flush=True)
    for n in names:
        th = threading.Thread(target=w, args=(n,)); th.start(); ths.append(th)
    for th in ths: th.join()
    summarize(names)
    print("PER_POSITION_DONE", flush=True)


if __name__ == "__main__":
    main()
