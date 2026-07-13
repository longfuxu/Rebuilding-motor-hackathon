#!/usr/bin/env bash
# Long-running watcher for the detached A100 umbrella run. Every INTERVAL: verify the VM is alive,
# tail umb.log, count windows, and INCREMENTALLY download ~/out/umb -> outputs/php_cycle/C3_umbrella
# (Spot-preemption salvage). Exits (one completion notification) when the run finishes ("UMBRELLA
# DONE"), the python dies, or the VM is gone. Interim progress: Read this script's task output file.
set -uo pipefail
GC="$HOME/google-cloud-sdk/bin/gcloud"; PROJECT="longfu-protein-gpu"; GCP_VM="${GCP_VM:-protein-gpu}"
HERE="$(cd "$(dirname "$0")" && pwd)"
ZONE="${GCP_ZONE:-$(cat "$HERE/.umbrella_zone" 2>/dev/null || echo us-central1-a)}"
DEST="$HERE/../../outputs/php_cycle/C3_umbrella"; mkdir -p "$DEST"
INTERVAL="${INTERVAL:-600}"
S(){ $GC compute ssh "$GCP_VM" --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap --command "$1" 2>/dev/null; }
CP(){ $GC compute scp --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap --recurse "$@" 2>/dev/null; }
ts(){ date -u +%H:%M:%S; }

reason="timeout"
for cycle in $(seq 1 "${MAXCYCLES:-120}"); do
  alive=$($GC compute instances list --project "$PROJECT" --filter="name=$GCP_VM" --format='value(name,status)' 2>/dev/null || true)
  if [ -z "$alive" ]; then echo "[$(ts)] VM GONE (Spot preemption or deleted) after cycle $cycle — salvaging local snapshot"; reason="gone"; break; fi
  lastlog=$(S 'tail -n 1 ~/umb.log 2>/dev/null')
  nwin=$(S '~/be/bin/python -c "import json,os;print(len(json.load(open(os.path.expanduser(\"~/out/umb/window_data.json\")))[\"windows\"]))" 2>/dev/null || echo 0')
  CP "$GCP_VM:~/out/umb" "$DEST/" >/dev/null 2>&1 || true
  CP "$GCP_VM:~/umb.log" "$DEST/" >/dev/null 2>&1 || true
  echo "[$(ts)] cycle $cycle | $alive | windows=$nwin | log: $lastlog"
  if S 'grep -q "UMBRELLA DONE" ~/umb.log 2>/dev/null'; then echo "[$(ts)] UMBRELLA DONE detected"; reason="done"; break; fi
  # crash check: [t] bracket keeps pgrep from matching its own invoking shell (whose cmdline holds the pattern)
  running=$(S 'pgrep -c -f "[t]md_umbrella.py" 2>/dev/null || echo 0')
  if [ "${running:-0}" = "0" ]; then echo "[$(ts)] python NOT running and not DONE — crash? salvaging + exiting"; reason="crash"; break; fi
  sleep "$INTERVAL"
done
echo "=== final download ==="
CP "$GCP_VM:~/out/umb" "$DEST/" >/dev/null 2>&1 || true
CP "$GCP_VM:~/umb.log" "$DEST/" >/dev/null 2>&1 || true
ls -la "$DEST/umb/" 2>/dev/null || echo "(no umb dir downloaded)"
# Auto-teardown ONLY on confirmed success (data already downloaded above). On crash/timeout leave the
# VM up for investigation (still on cap=1 so it must be dealt with); on 'gone' there's nothing to delete.
if [ "$reason" = "done" ]; then
  echo "=== [$(ts)] run complete -> TEARDOWN $GCP_VM (cost safety) ==="
  GCP_ZONE="$ZONE" GCP_VM="$GCP_VM" python3 "$HERE/../../gcp_pipeline/gcp_gpu.py" down \
    || $GC compute instances delete "$GCP_VM" --project "$PROJECT" --zone "$ZONE" --quiet || true
  echo "--- instances now (expect 0) ---"
  $GC compute instances list --project "$PROJECT" --format='value(name,status)' 2>/dev/null || true
else
  echo "=== [$(ts)] reason=$reason -> NOT auto-deleting; VM may still be up (cap=1) — investigate/teardown manually ==="
fi
echo "WATCHER_EXIT reason=$reason"
