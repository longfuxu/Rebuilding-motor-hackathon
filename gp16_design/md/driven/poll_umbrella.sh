#!/usr/bin/env bash
# Poll the detached A100 umbrella run: tail the log, report windows done, and INCREMENTALLY
# download ~/out/umb -> outputs/php_cycle/C3_umbrella/ (Spot-preemption-safe snapshot every call).
set -uo pipefail
GC="$HOME/google-cloud-sdk/bin/gcloud"; PROJECT="longfu-protein-gpu"; GCP_VM="${GCP_VM:-protein-gpu}"
HERE="$(cd "$(dirname "$0")" && pwd)"
ZONE="${GCP_ZONE:-$(cat "$HERE/.umbrella_zone" 2>/dev/null || echo us-central1-a)}"
DEST="$HERE/../../outputs/php_cycle/C3_umbrella"; mkdir -p "$DEST"
S(){ $GC compute ssh "$GCP_VM" --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap --command "$1" 2>/dev/null; }
CP(){ $GC compute scp --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap --recurse "$@" 2>/dev/null; }

echo "=== instance state ==="
$GC compute instances list --project "$PROJECT" --format='value(name,status,zone)' 2>/dev/null || echo "(list failed)"
echo "=== is the python still running? ==="
S 'pgrep -af tmd_umbrella | head -1 || echo "NOT RUNNING (finished or crashed)"'
echo "=== umb.log tail ==="
S 'tail -n 18 ~/umb.log 2>/dev/null || echo "(no log yet)"'
echo "=== windows done so far ==="
S "~/be/bin/python - <<'PY' 2>/dev/null || echo '(no window_data.json yet)'
import json
d=json.load(open('/home/'+__import__('getpass').getuser()+'/out/umb/window_data.json'))
w=d['windows']; print(len(w),'windows;','last:',{k:round(v,2) for k,v in w[-1].items() if isinstance(v,(int,float))})
PY"
echo "=== incremental download -> $DEST ==="
CP "$GCP_VM:~/out/umb" "$DEST/" && echo "downloaded snapshot" || echo "(download skipped/failed)"
CP "$GCP_VM:~/umb.log" "$DEST/" || true
