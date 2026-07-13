#!/usr/bin/env bash
# Poll the detached ring+DNA animation MD; when done, download traj + final.pdb -> outputs/, then
# TEAR DOWN the VM (cost safety) and report. Exits on done / crash / VM-gone.
set -uo pipefail
GC="$HOME/google-cloud-sdk/bin/gcloud"; PROJECT="longfu-protein-gpu"; GCP_VM="${GCP_VM:-protein-gpu}"
HERE="$(cd "$(dirname "$0")" && pwd)"
ZONE="${GCP_ZONE:-$(cat "$HERE/.dnaanim_zone" 2>/dev/null || echo us-central1-a)}"
DEST="$HERE/../../outputs/php_cycle/dna_translocation"; mkdir -p "$DEST"
OUTNAME="${OUTNAME:-dna_anim}"
INTERVAL="${INTERVAL:-120}"
S(){ $GC compute ssh "$GCP_VM" --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap --command "$1" 2>/dev/null; }
CP(){ $GC compute scp --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap --recurse "$@" 2>/dev/null; }
ts(){ date -u +%H:%M:%S; }
reason="timeout"
for cycle in $(seq 1 60); do
  alive=$($GC compute instances list --project "$PROJECT" --filter="name=$GCP_VM" --format='value(name)' 2>/dev/null || true)
  [ -z "$alive" ] && { echo "[$(ts)] VM gone"; reason="gone"; break; }
  echo "[$(ts)] cycle $cycle | $(S 'tail -n 1 ~/dna.log 2>/dev/null')"
  if S 'grep -q "RESULT_TMD" ~/dna.log 2>/dev/null'; then echo "[$(ts)] DNA-anim DONE"; reason="done"; break; fi
  r=$(S 'pgrep -c -f "[t]md_staircase.py" 2>/dev/null || echo 0'); [ "${r:-0}" = "0" ] && { echo "[$(ts)] crash?"; reason="crash"; break; }
  sleep "$INTERVAL"
done
echo "=== download traj + final.pdb ==="
CP "$GCP_VM:~/out/$OUTNAME" "$DEST/" >/dev/null 2>&1 || true
CP "$GCP_VM:~/dna.log" "$DEST/${OUTNAME}.log" >/dev/null 2>&1 || true
ls -la "$DEST/$OUTNAME/" 2>/dev/null || echo "(no $OUTNAME dir)"
if [ "$reason" = "done" ] || [ "$reason" = "crash" ]; then
  echo "=== TEARDOWN $GCP_VM ==="
  GCP_ZONE="$ZONE" GCP_VM="$GCP_VM" python3 "$HERE/../../gcp_pipeline/gcp_gpu.py" down \
    || $GC compute instances delete "$GCP_VM" --project "$PROJECT" --zone "$ZONE" --quiet || true
  $GC compute instances list --project "$PROJECT" --format='value(name,status)' 2>/dev/null || true
fi
echo "DNA_ANIM_POLL_EXIT reason=$reason"
