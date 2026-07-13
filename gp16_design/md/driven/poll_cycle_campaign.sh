#!/usr/bin/env bash
# Poll the detached P->H->P cycle campaign; incrementally download all fwd_s*/rev_s* runs; on CAMPAIGN_DONE
# (or crash / VM-gone) do a final download and TEAR DOWN the VM. Interim: Read this task's output file.
set -uo pipefail
GC="$HOME/google-cloud-sdk/bin/gcloud"; PROJECT="longfu-protein-gpu"; GCP_VM="${GCP_VM:-protein-gpu}"
HERE="$(cd "$(dirname "$0")" && pwd)"
ZONE="${GCP_ZONE:-$(cat "$HERE/.campaign_zone" 2>/dev/null || echo us-central1-a)}"
DEST="$HERE/../../outputs/php_cycle/cycle_campaign"; mkdir -p "$DEST"
INTERVAL="${INTERVAL:-180}"
S(){ $GC compute ssh "$GCP_VM" --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap --command "$1" 2>/dev/null; }
CP(){ $GC compute scp --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap --recurse "$@" 2>/dev/null; }
ts(){ date -u +%H:%M:%S; }
reason="timeout"
for cycle in $(seq 1 100); do
  alive=$($GC compute instances list --project "$PROJECT" --filter="name=$GCP_VM" --format='value(name)' 2>/dev/null || true)
  [ -z "$alive" ] && { echo "[$(ts)] VM gone"; reason="gone"; break; }
  last=$(S 'tail -n 1 ~/campaign.log 2>/dev/null'); done_dirs=$(S 'ls -d ~/out/*_s* 2>/dev/null | wc -l')
  CP "$GCP_VM:~/out/*_s*" "$DEST/" >/dev/null 2>&1 || true
  CP "$GCP_VM:~/campaign.log" "$DEST/" >/dev/null 2>&1 || true
  echo "[$(ts)] cycle $cycle | run dirs=$done_dirs | log: $last"
  if S 'grep -q CAMPAIGN_DONE ~/campaign.log 2>/dev/null'; then echo "[$(ts)] CAMPAIGN_DONE"; reason="done"; break; fi
  r=$(S 'pgrep -c -f "[t]md_staircase.py" 2>/dev/null || echo 0')
  bash_alive=$(S 'pgrep -c -f "[c]ampaign.sh" 2>/dev/null || echo 0')
  [ "${r:-0}" = "0" ] && [ "${bash_alive:-0}" = "0" ] && { echo "[$(ts)] no python & no campaign.sh — crash?"; reason="crash"; break; }
  sleep "$INTERVAL"
done
echo "=== final download ==="
CP "$GCP_VM:~/out/*_s*" "$DEST/" >/dev/null 2>&1 || true
CP "$GCP_VM:~/campaign.log" "$DEST/" >/dev/null 2>&1 || true
ls -d "$DEST"/*_s* 2>/dev/null
if [ "$reason" = "done" ] || [ "$reason" = "crash" ]; then
  echo "=== TEARDOWN $GCP_VM ==="
  GCP_ZONE="$ZONE" GCP_VM="$GCP_VM" python3 "$HERE/../../gcp_pipeline/gcp_gpu.py" down \
    || $GC compute instances delete "$GCP_VM" --project "$PROJECT" --zone "$ZONE" --quiet || true
  $GC compute instances list --project "$PROJECT" --format='value(name,status)' 2>/dev/null || true
fi
echo "CAMPAIGN_POLL_EXIT reason=$reason"
