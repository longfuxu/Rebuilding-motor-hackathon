#!/usr/bin/env python3
"""Fold one fusion-ladder construct via the free NVIDIA BioNeMo Boltz-2 NIM.
Saves the predicted mmCIF + a metrics JSON into BOTH the repo and the CS workspace,
so Claude Science can read the same result. Uses NVIDIA free credits (no Modal).

Usage:  NVIDIA_API_KEY=... python3 nim_fold.py <name>    e.g. trimer | tetramer | pentamer
"""
import sys, os, json, time, urllib.request

NAME = sys.argv[1]
KEY  = os.environ["NVIDIA_API_KEY"]
URL  = "https://health.api.nvidia.com/v1/biology/mit/boltz2/predict"

REPO = "/Users/longfu/Developer/claude-science-hackthon/.claude/worktrees/export-gp16-results/gp16_design"
CSWS = "/Users/longfu/.claude-science/orgs/34666d1b-db58-4bbc-a186-ff6de1eed1e6/workspaces/9c944ce5-d054-499c-9fbd-2be1bcf2e0cb"

fa = f"{REPO}/fold_inputs/ladder/gp16_ladder_{NAME}.fasta"
seq = "".join(l.strip() for l in open(fa) if not l.startswith(">"))
payload = json.dumps({"polymers": [{"id": "A", "molecule_type": "protein", "sequence": seq}]}).encode()

t0 = time.time()
req = urllib.request.Request(URL, data=payload,
        headers={"Authorization": f"Bearer {KEY}", "Content-Type": "application/json"})
try:
    resp = json.load(urllib.request.urlopen(req, timeout=1800))
except urllib.error.HTTPError as e:
    print(f"HTTP ERROR {e.code}: {e.read()[:400].decode(errors='replace')}")
    sys.exit(1)
dt = time.time() - t0

cif = resp["structures"][0]["structure"]
metrics = {k: resp.get(k) for k in ("confidence_scores", "ptm_scores", "iptm_scores",
                                    "complex_plddt_scores", "complex_pde_scores")}
metrics["name"] = NAME; metrics["len_aa"] = len(seq); metrics["wall_s"] = round(dt, 1)
metrics["predictor"] = "boltz2-nim"; metrics["msa"] = "single-seq (screen)"

for base in (f"{REPO}/outputs/structures/ladder", f"{CSWS}/outputs/structures/ladder"):
    os.makedirs(base, exist_ok=True)
    open(f"{base}/gp16_ladder_{NAME}__boltz2nim.cif", "w").write(cif)
    json.dump(metrics, open(f"{base}/gp16_ladder_{NAME}__boltz2nim_metrics.json", "w"), indent=2)

print(f"OK {NAME}: {len(seq)} aa in {dt:.0f}s | "
      f"conf={metrics['confidence_scores']} ptm={metrics['ptm_scores']} "
      f"plddt={metrics['complex_plddt_scores']}")
print("saved CIF+metrics to repo and CS workspace: outputs/structures/ladder/")
