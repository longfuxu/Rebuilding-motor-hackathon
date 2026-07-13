#!/usr/bin/env bash
# Poll the detached ratchet run (both modes); download ratchet_concerted/ + ratchet_sequential/; on
# RATCHET_ALL_DONE (or crash / VM-gone) final-download + TEAR DOWN the VM.
set -uo pipefail
GC="$HOME/google-cloud-sdk/bin/gcloud"; PROJECT="longfu-protein-gpu"; GCP_VM="${GCP_VM:-protein-gpu}"
HERE="$(cd "$(dirname "$0")" && pwd)"
ZONE="${GCP_ZONE:-$(cat "$HERE/.ratchet_zone" 2>/dev/null || echo us-central1-a)}"
DEST="$HERE/../../outputs/php_cycle/ratchet"; mkdir -p "$DEST"; INTERVAL="${INTERVAL:-180}"
S(){ $GC compute ssh "$GCP_VM" --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap --command "$1" 2>/dev/null; }
CP(){ $GC compute scp --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap --recurse "$@" 2>/dev/null; }
ts(){ date -u +%H:%M:%S; }
reason="timeout"
for cycle in $(seq 1 80); do
  alive=$($GC compute instances list --project "$PROJECT" --filter="name=$GCP_VM" --format='value(name)' 2>/dev/null || true)
  [ -z "$alive" ] && { echo "[$(ts)] VM gone"; reason="gone"; break; }
  echo "[$(ts)] cycle $cycle | $(S 'tail -n 1 ~/ratchet.log 2>/dev/null')"
  CP "$GCP_VM:~/out/ratchet_*" "$DEST/" >/dev/null 2>&1 || true
  if S 'grep -q RATCHET_ALL_DONE ~/ratchet.log 2>/dev/null'; then echo "[$(ts)] RATCHET_ALL_DONE"; reason="done"; break; fi
  r=$(S 'pgrep -c -f "[t]md_ratchet.py" 2>/dev/null || echo 0'); b=$(S 'pgrep -c -f "[r]atchet.sh" 2>/dev/null || echo 0')
  [ "${r:-0}" = "0" ] && [ "${b:-0}" = "0" ] && { echo "[$(ts)] crash?"; reason="crash"; break; }
  sleep "$INTERVAL"
done
CP "$GCP_VM:~/out/ratchet_*" "$DEST/" >/dev/null 2>&1 || true
CP "$GCP_VM:~/ratchet.log" "$DEST/" >/dev/null 2>&1 || true
ls -d "$DEST"/ratchet_* 2>/dev/null
if [ "$reason" = "done" ] || [ "$reason" = "crash" ]; then
  echo "=== TEARDOWN $GCP_VM ==="
  GCP_ZONE="$ZONE" GCP_VM="$GCP_VM" python3 "$HERE/../../gcp_pipeline/gcp_gpu.py" down \
    || $GC compute instances delete "$GCP_VM" --project "$PROJECT" --zone "$ZONE" --quiet || true
  $GC compute instances list --project "$PROJECT" --format='value(name,status)' 2>/dev/null || true
fi
echo "RATCHET_POLL_EXIT reason=$reason"
