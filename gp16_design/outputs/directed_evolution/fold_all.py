#!/usr/bin/env python3
"""Batch-fold every directed-evolution variant manifest through the tiled-MSA Boltz-2
NIM pipeline (the M1/M2 gate). Serial (avoids NIM rate limits), resumable (skips folds
whose <name>.result.json already reports ok:true). Reuses the existing validated WT fold.
"""
import os, sys, json, glob, shutil, time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
PIPE = os.path.join(REPO, "pipelines", "tiled_msa_fold")
sys.path.insert(0, PIPE)
MANI = os.path.join(HERE, "manifests")
FOLDS = os.path.join(HERE, "folds")
# .env lives in the MAIN checkout (worktrees don't carry it); try a few locations
ENV_CANDIDATES = [
    "/Users/longfu/Developer/claude-science-hackthon/.env",
    os.path.abspath(os.path.join(REPO, "..", ".env")),
]


def load_env_key():
    if os.environ.get("NVIDIA_API_KEY"):
        return
    for env in ENV_CANDIDATES:
        if not os.path.exists(env):
            continue
        for line in open(env):
            line = line.strip()
            if line.startswith("NVIDIA_API_KEY"):
                os.environ["NVIDIA_API_KEY"] = line.split("=", 1)[1].strip().strip('"').strip("'")
                return
    raise RuntimeError("NVIDIA_API_KEY not found in env or .env candidates")


def done(name):
    rj = os.path.join(FOLDS, f"{name}.result.json")
    if not os.path.exists(rj):
        return False
    try:
        return bool(json.load(open(rj)).get("ok"))
    except Exception:
        return False


def main():
    os.makedirs(FOLDS, exist_ok=True)
    load_env_key()
    import tiled_msa_fold as tf

    # reuse the existing validated WT fold for cp233_WT (identical sequence + settings)
    src = os.path.join(REPO, "outputs", "tiled_fold")
    for ext in ("cif", "result.json"):
        s = os.path.join(src, f"cp233_WT.{ext}")
        d = os.path.join(FOLDS, f"cp233_WT.{ext}")
        if os.path.exists(s) and not os.path.exists(d):
            shutil.copy(s, d)
            print(f"reused existing WT fold: {ext}")

    manifests = sorted(glob.glob(os.path.join(MANI, "*.json")))
    for mp in manifests:
        name = json.load(open(mp))["name"]
        if done(name):
            print(f"[skip] {name} already folded", flush=True)
            continue
        print(f"[fold] {name} ...", flush=True)
        t0 = time.time()
        try:
            res = tf.do_run(mp, FOLDS)
            print(f"[fold] {name} -> ok={res.get('ok')} in {round(time.time()-t0)}s", flush=True)
        except Exception as e:
            print(f"[fold] {name} FAILED: {repr(e)[:160]}", flush=True)
        time.sleep(3)   # be polite to the NIM
    print("ALL FOLDS DONE", flush=True)


if __name__ == "__main__":
    main()
