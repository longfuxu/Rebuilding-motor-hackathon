#!/usr/bin/env python3
"""Batch-fold every round-2 directed-evolution manifest through the tiled-MSA Boltz-2 NIM.
Serial (avoids NIM rate limits), resumable (skips folds whose <name>.result.json reports
ok:true). Reuses the round-1 WT + WTrep folds (copied into round2/folds already).

Each fold is written to disk immediately by tiled_msa_fold.do_run (the API has been flaky,
so nothing is held in memory). Rotates between the two NVIDIA keys on failure.
"""
import os, sys, json, glob, time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))   # gp16_design/
PIPE = os.path.join(REPO, "pipelines", "tiled_msa_fold")
sys.path.insert(0, PIPE)
MANI = os.path.join(HERE, "manifests")
FOLDS = os.path.join(HERE, "folds")
ENV_CANDIDATES = [
    "/Users/longfu/Developer/claude-science-hackthon/.env",
    os.path.abspath(os.path.join(REPO, "..", ".env")),
]


def load_keys():
    keys = []
    for env in ENV_CANDIDATES:
        if not os.path.exists(env):
            continue
        for line in open(env):
            line = line.strip()
            for k in ("NVIDIA_API_KEY", "NVIDIA_API_KEY2"):
                if line.startswith(k + "="):
                    v = line.split("=", 1)[1].strip().strip('"').strip("'")
                    if v and v not in keys:
                        keys.append(v)
        if keys:
            break
    if not keys:
        raise RuntimeError("no NVIDIA_API_KEY found in .env candidates")
    return keys


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
    keys = load_keys()
    print(f"loaded {len(keys)} NVIDIA key(s)", flush=True)
    import tiled_msa_fold as tf

    manifests = sorted(glob.glob(os.path.join(MANI, "*.json")))
    todo = [m for m in manifests if not done(json.load(open(m))["name"])]
    print(f"{len(manifests)} manifests, {len(todo)} to fold", flush=True)

    for mp in manifests:
        name = json.load(open(mp))["name"]
        if done(name):
            print(f"[skip] {name} already folded", flush=True)
            continue
        ok = False
        for ki, key in enumerate(keys):
            os.environ["NVIDIA_API_KEY"] = key
            print(f"[fold] {name} (key #{ki+1}) ...", flush=True)
            t0 = time.time()
            try:
                res = tf.do_run(mp, FOLDS)
                if res.get("ok"):
                    print(f"[fold] {name} -> ok in {round(time.time()-t0)}s", flush=True)
                    ok = True
                    break
                else:
                    print(f"[fold] {name} not-ok: {str(res.get('err'))[:160]}", flush=True)
            except Exception as e:
                print(f"[fold] {name} EXC (key #{ki+1}): {repr(e)[:160]}", flush=True)
            time.sleep(5)
        if not ok:
            print(f"[fold] {name} FAILED on all keys -- will retry next run", flush=True)
        time.sleep(3)
    remaining = [json.load(open(m))["name"] for m in manifests
                 if not done(json.load(open(m))["name"])]
    print(f"ALL FOLDS DONE. remaining unfolded: {remaining}", flush=True)


if __name__ == "__main__":
    main()
