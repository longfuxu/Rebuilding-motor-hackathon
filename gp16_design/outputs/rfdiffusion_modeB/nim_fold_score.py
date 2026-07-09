#!/usr/bin/env python3
"""Fold designed single-chain connector constructs with the free NVIDIA Boltz-2 NIM,
then score the A->B trans interface with reproduce/score_m2.py.

Input: a directory of designed PDBs (RFdiffusion+MPNN 'best_design*.pdb', full side chains)
OR a FASTA of full designed sequences. Layout is deterministic: subunit A = out res 1..327
(gp16 4-330), connector = middle, subunit B = last 327 -> copies "A:1-327,A:{N-326}-{N}"
with --copy_start_res 4. M2 = R146(A)->WalkerA(B), engaged < 8A.

Usage: NVIDIA_API_KEY=... python nim_fold_score.py <indir_or_fasta> <outdir> [max_concurrent]
"""
import sys, os, json, time, urllib.request, urllib.error, threading, glob

KEY = os.environ["NVIDIA_API_KEY"]
URL = "https://health.api.nvidia.com/v1/biology/mit/boltz2/predict"
REPO = "/Users/longfu/Developer/claude-science-hackthon/.claude/worktrees/export-gp16-results/gp16_design"
SCORE = f"{REPO}/reproduce/score_m2.py"
AA3TO1 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
          "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
          "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}


def seq_from_pdb(path):
    seq = {}
    for L in open(path):
        if L[:4] == "ATOM" and L[12:16].strip() == "CA":
            seq[int(L[22:26])] = AA3TO1.get(L[17:20].strip(), "X")
    return "".join(seq[k] for k in sorted(seq))


def load_inputs(src):
    items = {}
    if os.path.isdir(src):
        for p in sorted(glob.glob(os.path.join(src, "*.pdb"))):
            items[os.path.splitext(os.path.basename(p))[0]] = seq_from_pdb(p)
    else:
        name = None
        for L in open(src):
            if L.startswith(">"):
                name = L[1:].strip().split()[0]
                items[name] = ""
            elif name:
                items[name] += L.strip()
    return items


def fold_one(name, seq, outdir, results, lock):
    payload = json.dumps({"polymers": [{"id": "A", "molecule_type": "protein", "sequence": seq}]}).encode()
    req = urllib.request.Request(URL, data=payload,
                                 headers={"Authorization": f"Bearer {KEY}", "Content-Type": "application/json"})
    t0 = time.time()
    try:
        resp = json.load(urllib.request.urlopen(req, timeout=1800))
    except urllib.error.HTTPError as e:
        with lock:
            results[name] = {"ok": False, "err": f"HTTP {e.code}: {e.read()[:300].decode(errors='replace')}"}
        return
    except Exception as e:
        with lock:
            results[name] = {"ok": False, "err": repr(e)[:300]}
        return
    dt = round(time.time() - t0, 1)
    # NIM Boltz-2 response format (confirmed via fold_scripts/nim_fold.py)
    try:
        struct = resp["structures"][0]["structure"]
    except (KeyError, IndexError, TypeError):
        open(os.path.join(outdir, f"{name}.raw.json"), "w").write(json.dumps(resp)[:5000])
        with lock:
            results[name] = {"ok": False, "err": "no structures[0].structure", "keys": list(resp.keys()), "dt": dt}
        return
    metrics = {k: resp.get(k) for k in ("confidence_scores", "ptm_scores", "iptm_scores",
                                        "complex_plddt_scores", "complex_pde_scores")}
    cif = os.path.join(outdir, f"{name}.cif")
    open(cif, "w").write(struct)
    N = len(seq)
    copies = f"A:1-327,A:{N-326}-{N}"
    import subprocess
    sc = subprocess.run(["python", SCORE, cif, "--copies", copies, "--copy_start_res", "4"],
                        capture_output=True, text=True)
    with lock:
        results[name] = {"ok": True, "N": N, "clen": N - 654, "dt": dt, "metrics": metrics,
                         "score_m2": sc.stdout, "score_err": sc.stderr[-300:]}
        json.dump(results, open(os.path.join(outdir, "fold_score.json"), "w"), indent=2)
    print(f"[{name}] folded N={N} clen={N-654} in {dt}s")


def main():
    src, outdir = sys.argv[1], sys.argv[2]
    maxc = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    os.makedirs(outdir, exist_ok=True)
    items = load_inputs(src)
    print(f"folding {len(items)} sequences, concurrency {maxc}")
    results, lock, sem = {}, threading.Lock(), threading.Semaphore(maxc)
    threads = []
    def worker(n, s):
        with sem:
            fold_one(n, s, outdir, results, lock)
    for n, s in items.items():
        t = threading.Thread(target=worker, args=(n, s)); t.start(); threads.append(t)
    for t in threads:
        t.join()
    json.dump(results, open(os.path.join(outdir, "fold_score.json"), "w"), indent=2)
    ok = sum(1 for r in results.values() if r.get("ok"))
    print(f"DONE: {ok}/{len(items)} folded")


if __name__ == "__main__":
    main()
