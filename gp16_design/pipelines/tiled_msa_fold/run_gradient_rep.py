#!/usr/bin/env python3
"""Dead-seat coordination gradient with REPLICATES: fold each construct with
diffusion_samples=3 (the NIM returns 3 structures per call) and score every replicate
handedness-robust. Reporting engaged as a distribution defeats the single-sample
ring-handedness/closure noise. Saves the best (max-engaged) structure as <name>.cif."""
import importlib.util, os, glob, json, threading

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("t", os.path.join(HERE, "tiled_msa_fold.py"))
t = importlib.util.module_from_spec(spec); spec.loader.exec_module(t)

REPO = "/Users/longfu/Developer/claude-science-hackthon/.claude/worktrees/export-gp16-results/gp16_design"
MDIR = f"{REPO}/outputs/deadseat_gradient/manifests"
OUT = f"{REPO}/outputs/deadseat_gradient"
NREP = 3
key = os.environ["NVIDIA_API_KEY"]
lock = threading.Lock()


def run_one(name):
    m = t.load_manifest(f"{MDIR}/{name}.json")
    a3m, meta = t.do_tile(m)
    resp, dt = t.fold_nim(m["sequence"], a3m, key, samples=NREP)
    if "__err" in resp:
        with lock:
            print(f"[{name}] FOLD ERR {resp['__err']}", flush=True)
        json.dump({"name": name, "ok": False, "err": resp["__err"], **meta},
                  open(f"{OUT}/{name}.result.json", "w"), indent=2)
        return
    reps = []
    for i, st in enumerate(resp.get("structures", [])):
        cif = f"{OUT}/{name}_rep{i}.cif"; open(cif, "w").write(st["structure"])
        v, out = t.score_ring(cif, meta["copies"], meta["r146_incopy"], meta["walker_incopy"])
        m1 = next((l for l in out.splitlines() if l.startswith("# M1:")), "")
        import re
        cv = (re.search(r"radius_CV ([\d.]+)", m1) or [0, ""])[1]
        rad = (re.search(r"radius ([\d.]+)", m1) or [0, ""])[1]
        seq = "sequential_consistent: YES" in m1
        reps.append({"rep": i, "engaged": v["engaged"], "handedness": v["handedness"],
                     "m2_forward": v["m2_forward"], "m2_reverse": v["m2_reverse"],
                     "radius_CV": cv, "radius": rad, "sequential": seq,
                     "score_full": out, "M1": m1})
    eng = [r["engaged"] for r in reps]
    best = max(reps, key=lambda r: (r["engaged"], -float(r["radius_CV"] or 9)))
    # keep best structure as <name>.cif; drop the others
    if os.path.exists(f"{OUT}/{name}_rep{best['rep']}.cif"):
        os.replace(f"{OUT}/{name}_rep{best['rep']}.cif", f"{OUT}/{name}.cif")
    for r in reps:
        p = f"{OUT}/{name}_rep{r['rep']}.cif"
        if os.path.exists(p):
            os.remove(p)
    res = {"name": name, "ok": True, "wall_s": dt, "nrep": len(reps),
           "engaged_values": eng, "engaged_max": max(eng), "engaged_min": min(eng),
           "engaged_mean": round(sum(eng) / len(eng), 2),
           "best_rep": best["rep"], "best_engaged": best["engaged"],
           "best_handedness": best["handedness"], "best_radius_CV": best["radius_CV"],
           "best_radius": best["radius"], "best_sequential": best["sequential"],
           "best_M1": best["M1"], "best_score_full": best["score_full"],
           "replicates": [{k: r[k] for k in ("rep", "engaged", "handedness", "radius_CV",
                                             "sequential")} for r in reps], **meta}
    json.dump(res, open(f"{OUT}/{name}.result.json", "w"), indent=2)
    with lock:
        print(f"[{name}] engaged {eng} (max {max(eng)}) best_CV {best['radius_CV']} "
              f"seq {best['sequential']} {dt}s", flush=True)


def main():
    names = sorted(os.path.splitext(os.path.basename(p))[0] for p in glob.glob(f"{MDIR}/*.json"))
    print(f"folding {len(names)} constructs x {NREP} replicates (concurrency 3)", flush=True)
    sem = threading.Semaphore(3); threads = []
    def worker(n):
        with sem:
            try:
                run_one(n)
            except Exception as e:
                with lock:
                    print(f"[{n}] EXC {repr(e)[:140]}", flush=True)
    for n in names:
        th = threading.Thread(target=worker, args=(n,)); th.start(); threads.append(th)
    for th in threads:
        th.join()
    print("GRADIENT_REP_DONE", flush=True)


if __name__ == "__main__":
    main()
