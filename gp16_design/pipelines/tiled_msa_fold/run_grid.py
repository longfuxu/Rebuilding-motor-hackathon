#!/usr/bin/env python3
"""Drive the 27-construct CP linker-grid screen through tiled_msa_fold.
Reuses the 3 already-folded int15_inter10 CIFs; folds the other 24 with light
concurrency; handedness-robust scoring; drops the big .tiled.a3m to save space."""
import importlib.util, os, glob, json, shutil, threading

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("t", os.path.join(HERE, "tiled_msa_fold.py"))
t = importlib.util.module_from_spec(spec); spec.loader.exec_module(t)

REPO = "/Users/longfu/Developer/claude-science-hackthon/.claude/worktrees/export-gp16-results/gp16_design"
MDIR = f"{REPO}/outputs/cp_grid_screen/manifests"
OUT = f"{REPO}/outputs/cp_grid_screen"
TILED = f"{REPO}/outputs/tiled_fold"
REUSE = {"cp233_int15_inter10": f"{TILED}/cp233_WT.cif",
         "cp285_int15_inter10": f"{TILED}/cp285_int15_inter10.cif",
         "cp297_int15_inter10": f"{TILED}/cp297_int15_inter10.cif"}
os.makedirs(OUT, exist_ok=True)
key = os.environ["NVIDIA_API_KEY"]
lock = threading.Lock()


def reuse(name, cif_src):
    m = t.load_manifest(f"{MDIR}/{name}.json")
    _, meta = t.do_tile(m)                                    # meta only (copies/landmarks); no fold
    cif = f"{OUT}/{name}.cif"; shutil.copy(cif_src, cif)
    v, out = t.score_ring(cif, meta["copies"], meta["r146_incopy"], meta["walker_incopy"])
    res = {"name": name, "ok": True, "reused_from": cif_src, "ring": v,
           "M2": f"# M2(ring): {v['engaged']}/{v['n']} [{v['handedness']}] fwd{v['m2_forward']} rev{v['m2_reverse']}",
           "M1": next((l for l in out.splitlines() if l.startswith("# M1:")), ""),
           "M4": next((l for l in out.splitlines() if l.startswith("# M4:")), ""),
           "score_full": out, **meta}
    json.dump(res, open(f"{OUT}/{name}.result.json", "w"), indent=2)
    with lock:
        print(f"[reuse {name}] {res['M2']}", flush=True)


def fold(name):
    res = t.do_run(f"{MDIR}/{name}.json", OUT, samples=1)
    a3m = f"{OUT}/{name}.tiled.a3m"
    if os.path.exists(a3m):
        os.remove(a3m)                                        # regenerable via `tile`; save space
    return res


def main():
    names = sorted(os.path.splitext(os.path.basename(p))[0] for p in glob.glob(f"{MDIR}/*.json"))
    todo = [n for n in names if n not in REUSE]
    for n, src in REUSE.items():
        if os.path.exists(src):
            reuse(n, src)
        else:
            todo.append(n)
    print(f"folding {len(todo)} constructs (concurrency 3)", flush=True)
    sem = threading.Semaphore(3); threads = []
    def worker(n):
        with sem:
            try:
                fold(n)
            except Exception as e:
                with lock:
                    print(f"[{n}] EXC {repr(e)[:120]}", flush=True)
    for n in todo:
        th = threading.Thread(target=worker, args=(n,)); th.start(); threads.append(th)
    for th in threads:
        th.join()
    print("GRID_DONE", flush=True)


if __name__ == "__main__":
    main()
