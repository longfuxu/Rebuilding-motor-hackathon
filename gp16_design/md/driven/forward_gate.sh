#!/usr/bin/env bash
# Forward-pass GATE detector for the both=1 umbrella run. Exits (one notification) the moment the
# FORWARD pass is complete (NWIN forward windows present) so the caller can check forward-PMF quality
# and decide whether to let the reverse pass continue or kill it. Also exits on VM-gone / python-crash.
# Does NOT download (the main watcher owns downloads) and does NOT tear anything down.
set -uo pipefail
GC="$HOME/google-cloud-sdk/bin/gcloud"; PROJECT="longfu-protein-gpu"; GCP_VM="${GCP_VM:-protein-gpu}"
HERE="$(cd "$(dirname "$0")" && pwd)"
ZONE="${GCP_ZONE:-$(cat "$HERE/.umbrella_zone" 2>/dev/null || echo us-central1-a)}"
NWIN="${NWIN:-24}"; INTERVAL="${INTERVAL:-300}"
S(){ $GC compute ssh "$GCP_VM" --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap --command "$1" 2>/dev/null; }
ts(){ date -u +%H:%M:%S; }
for cycle in $(seq 1 120); do
  alive=$($GC compute instances list --project "$PROJECT" --filter="name=$GCP_VM" --format='value(name)' 2>/dev/null || true)
  if [ -z "$alive" ]; then echo "[$(ts)] GATE ABORT: VM gone (preemption/deleted)"; exit 0; fi
  nfwd=$(S "~/be/bin/python -c \"import json,os;w=json.load(open(os.path.expanduser('~/out/umb/window_data.json')))['windows'];print(sum(1 for x in w if x.get('pass_')=='fwd'))\" 2>/dev/null || echo 0")
  echo "[$(ts)] cycle $cycle | forward windows = ${nfwd:-0}/$NWIN"
  if [ "${nfwd:-0}" -ge "$NWIN" ] 2>/dev/null; then echo "[$(ts)] FORWARD_COMPLETE ($nfwd/$NWIN) — gate: check forward PMF now"; exit 0; fi
  if S 'grep -q "UMBRELLA DONE" ~/umb.log 2>/dev/null'; then echo "[$(ts)] run already DONE before gate fired"; exit 0; fi
  running=$(S 'pgrep -c -f "[t]md_umbrella.py" 2>/dev/null || echo 0')
  if [ "${running:-0}" = "0" ]; then echo "[$(ts)] GATE ABORT: python not running (crash?)"; exit 0; fi
  sleep "$INTERVAL"
done
echo "[$(ts)] GATE timed out"
