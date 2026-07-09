#!/usr/bin/env python3
"""Tile each sal_L50 de-novo connector 5x into a single-chain native-order ring,
fold with Boltz-2 NIM, score sequential M2 across the 5 copies.

Ring = 5 x [gp16 res4-330 (327 aa)] joined by 4 x [designed L50 connector (50 aa)]
     = 5*327 + 4*50 = 1835 aa.
Copies (each subunit gp16 4-330, --copy_start_res 4):
  A:1-327, A:378-704, A:755-1081, A:1132-1458, A:1509-1835
"""
import os, sys, json, time, urllib.request, urllib.error, subprocess

KEY = os.environ["NVIDIA_API_KEY"]
URL = "https://health.api.nvidia.com/v1/biology/mit/boltz2/predict"
REPO = "/Users/longfu/Developer/claude-science-hackthon/.claude/worktrees/export-gp16-results/gp16_design"
SCORE = f"{REPO}/reproduce/score_m2.py"
GP16 = open("/Users/longfu/.claude/jobs/48b7a543/tmp/gp16_4_330.seq").read().strip()  # 327 aa
FASTA = "/Users/longfu/.claude/jobs/48b7a543/tmp/mpnn_work/all_designs.fasta"
OUT = "/Users/longfu/.claude/jobs/48b7a543/tmp/tiled_ring_out"
os.makedirs(OUT, exist_ok=True)
COPIES = "A:1-327,A:378-704,A:755-1081,A:1132-1458,A:1509-1835"


def load_designs():
    recs = []; cur = None
    for l in open(FASTA):
        if l.startswith(">"): cur = l[1:].strip().split()[0]; recs.append([cur, ""])
        elif recs: recs[-1][1] += l.strip()
    return recs


def tile(design_seq):
    N = len(design_seq)
    # subunit A = 0..326, connector = 327..(N-327-1), subunit B = N-327..N-1
    connector = design_seq[327:N - 327]
    assert len(connector) == N - 654, (len(connector), N)
    subs = [GP16] * 5
    ring = subs[0]
    for k in range(1, 5):
        ring += connector + subs[k]
    return ring, len(connector)


def fold(name, seq):
    payload = json.dumps({"polymers": [{"id": "A", "molecule_type": "protein", "sequence": seq}]}).encode()
    for attempt in range(5):
        try:
            req = urllib.request.Request(URL, data=payload,
                                         headers={"Authorization": f"Bearer {KEY}", "Content-Type": "application/json"})
            t0 = time.time()
            resp = json.load(urllib.request.urlopen(req, timeout=1800))
            return resp, round(time.time() - t0, 1)
        except urllib.error.HTTPError as e:
            if e.code == 429:
                print(f"  [{name}] 429, backoff", flush=True); time.sleep(20 * (attempt + 1)); continue
            return {"__err": f"HTTP {e.code}: {e.read()[:200].decode(errors='replace')}"}, 0
        except Exception as e:
            print(f"  [{name}] {repr(e)[:80]}, retry", flush=True); time.sleep(10)
    return {"__err": "429 after retries"}, 0


def main():
    only_l50 = [r for r in load_designs() if r[0].startswith("sal_L50")]
    print(f"tiling {len(only_l50)} L50 designs into 1835-aa rings", flush=True)
    results = {}
    for name, ds in only_l50:
        ring, clen = tile(ds)
        print(f"[{name}] ring {len(ring)} aa (connector {clen})", flush=True)
        resp, dt = fold(name, ring)
        if "__err" in resp:
            results[name] = {"ok": False, "err": resp["__err"]}; print(f"  FAIL {resp['__err']}", flush=True)
            json.dump(results, open(f"{OUT}/ring_results.json", "w"), indent=1); continue
        cif = f"{OUT}/{name}_ring.cif"
        open(cif, "w").write(resp["structures"][0]["structure"])
        plddt = resp.get("complex_plddt_scores")
        sc = subprocess.run(["python", SCORE, cif, "--copies", COPIES, "--copy_start_res", "4"],
                            capture_output=True, text=True)
        m2 = [l for l in sc.stdout.split("\n") if "# M2:" in l]
        m1 = [l for l in sc.stdout.split("\n") if "# M1:" in l]
        results[name] = {"ok": True, "ring_aa": len(ring), "dt": dt, "plddt": plddt,
                         "M2": m2[0] if m2 else "", "M1": m1[0] if m1 else "", "full": sc.stdout}
        print(f"  {m2[0] if m2 else '?'}", flush=True)
        print(f"  {m1[0] if m1 else '?'}", flush=True)
        json.dump(results, open(f"{OUT}/ring_results.json", "w"), indent=1)
    print("TILED_RING_DONE", flush=True)


if __name__ == "__main__":
    main()
