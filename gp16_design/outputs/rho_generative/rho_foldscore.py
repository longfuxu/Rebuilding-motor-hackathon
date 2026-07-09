#!/usr/bin/env python3
"""Fold a Rho single-chain construct with Boltz-2 NIM, score sequential M2
(Rho arginine finger R366 -> neighbor Walker-A ~179-186), 6 copies."""
import sys, os, json, time, urllib.request, urllib.error, subprocess

KEY = os.environ["NVIDIA_API_KEY"]
URL = "https://health.api.nvidia.com/v1/biology/mit/boltz2/predict"
REPO = "/Users/longfu/Developer/claude-science-hackthon/.claude/worktrees/export-gp16-results/gp16_design"
SCORE = f"{REPO}/reproduce/score_m2.py"

fasta, outdir, copies, r_inc, wa = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
os.makedirs(outdir, exist_ok=True)
recs = []; cur = None
for l in open(fasta):
    if l.startswith(">"): cur = l[1:].strip().split()[0]; recs.append([cur, ""])
    elif recs: recs[-1][1] += l.strip()
results = {}
for name, seq in recs:
    print(f"folding {name} ({len(seq)} aa)...", flush=True)
    payload = json.dumps({"polymers": [{"id": "A", "molecule_type": "protein", "sequence": seq}]}).encode()
    resp = None
    for att in range(5):
        try:
            req = urllib.request.Request(URL, data=payload, headers={"Authorization": f"Bearer {KEY}", "Content-Type": "application/json"})
            resp = json.load(urllib.request.urlopen(req, timeout=1800)); break
        except urllib.error.HTTPError as e:
            if e.code == 429: time.sleep(20 * (att + 1)); continue
            results[name] = {"ok": False, "err": f"HTTP {e.code}: {e.read()[:200].decode(errors='replace')}"}; resp = "ERR"; break
        except Exception as e:
            if att == 4: results[name] = {"ok": False, "err": repr(e)[:150]}; resp = "ERR"
            time.sleep(10)
    if resp in (None, "ERR"):
        print(f"  FAIL {results.get(name)}", flush=True); json.dump(results, open(f"{outdir}/rho_scores.json", "w"), indent=1); continue
    cif = f"{outdir}/{name}.cif"; open(cif, "w").write(resp["structures"][0]["structure"])
    plddt = resp.get("complex_plddt_scores")
    sc = subprocess.run(["python", SCORE, cif, "--copies", copies, "--copy_start_res", "1",
                         "--r146_incopy", r_inc, "--walker_incopy", wa], capture_output=True, text=True)
    m2 = [l for l in sc.stdout.split("\n") if "# M2:" in l]
    m1 = [l for l in sc.stdout.split("\n") if "# M1:" in l]
    results[name] = {"ok": True, "N": len(seq), "plddt": plddt, "M2": m2[0] if m2 else "", "M1": m1[0] if m1 else "", "full": sc.stdout}
    print(f"  {m2[0] if m2 else '?'}\n  {m1[0] if m1 else '?'}", flush=True)
    json.dump(results, open(f"{outdir}/rho_scores.json", "w"), indent=1)
print("RHO_FOLDSCORE_DONE", flush=True)
