#!/usr/bin/env bash
# Ring+DNA animation MD: implicit-solvent targeted-MD driving the SINGLE-CHAIN ring planar->helical
# with dsDNA threaded through the pore (~31k atoms), saving a frequent-frame trajectory to render the
# "motor gripping/moving DNA" animation. Reuses the proven staircase CV (tmd_staircase.py) — the CV
# acts only on the ring's 5 subunits; the DNA (chains F,G) responds freely. Short (~200 ps pull ≈ ~1 h).
#
# *** cap=1 guard: run ONLY after the umbrella VM is gone. *** Non-blocking: launches detached, then
# poll_dna_anim (below) downloads + tears down. VM name protein-gpu (default).
set -euo pipefail
export GCP_VM="${GCP_VM:-protein-gpu}"; GPU="a100"
GC="$HOME/google-cloud-sdk/bin/gcloud"; PROJECT="longfu-protein-gpu"
HERE="$(cd "$(dirname "$0")" && pwd)"; GCPPY="$HERE/../../gcp_pipeline/gcp_gpu.py"
INPUT="${INPUT:-$HERE/inputs/C_plus_dna_relaxed.pdb}"   # planar structure (CV geometry)
POS_FROM="${POS_FROM:-}"                                 # optional helical endpoint for a reverse H->P pull
LAM0="${LAM0:-0}"; LAM1="${LAM1:-1}"; OUTNAME="${OUTNAME:-dna_anim}"; TRIANGLE="${TRIANGLE:-0}"
PULL_PS="${PULL_PS:-200}"; EQUIL_PS="${EQUIL_PS:-50}"; REPORT_PS="${REPORT_PS:-2}"; KCV="${KCV:-30000}"

test -f "$INPUT" || { echo "missing $INPUT"; exit 1; }
[ -n "$POS_FROM" ] && { test -f "$POS_FROM" || { echo "missing POS_FROM $POS_FROM"; exit 1; }; }
echo "=== cap=1 guard ==="
running=$($GC compute instances list --project "$PROJECT" --format='value(name)' 2>/dev/null || true)
[ -n "$running" ] && { echo "REFUSING: instance already up (cap=1): $running"; exit 3; }

echo "=== bring up Spot $GPU as $GCP_VM (zone fallback) ==="
ZONE=""
for Z in "${GCP_ZONE:-us-central1-a}" us-central1-b us-central1-c us-central1-f; do
  if GCP_ZONE="$Z" python3 "$GCPPY" up --gpu "$GPU"; then ZONE="$Z"; break; fi
done
: "${ZONE:?all Spot zones stocked out}"; echo "$ZONE" > "$HERE/.dnaanim_zone"
S(){ $GC compute ssh "$GCP_VM" --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap --command "$1"; }
CP(){ $GC compute scp --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap --recurse "$@"; }

echo "=== upload + install ==="
S 'mkdir -p ~/in/inputs ~/out'
CP "$HERE/tmd_staircase.py" "$GCP_VM:~/in/"
CP "$INPUT" "$GCP_VM:~/in/inputs/"
POS_ARG=""
if [ -n "$POS_FROM" ]; then CP "$POS_FROM" "$GCP_VM:~/in/inputs/"; POS_ARG="--pos_from inputs/$(basename "$POS_FROM")"; fi
S 'sudo apt-get -qq update >/dev/null 2>&1; sudo apt-get -qq install -y python3.10-venv python3-pip >/dev/null 2>&1; \
   test -d ~/be || python3 -m venv ~/be; ~/be/bin/pip -q install openmm numpy mdtraj 2>&1 | tail -1; \
   for i in $(seq 1 30); do nvidia-smi -L >/dev/null 2>&1 && break; sleep 15; done; nvidia-smi -L | head -1'

echo "=== launch ring+DNA targeted-MD DETACHED ($OUTNAME, lam $LAM0->$LAM1, report every ${REPORT_PS}ps -> ~$((PULL_PS/REPORT_PS)) frames) ==="
S "cd ~/in && rm -f ~/dna.log && nohup ~/be/bin/python tmd_staircase.py \
     --input inputs/$(basename "$INPUT") $POS_ARG --out ~/out/$OUTNAME \
     --lam0 $LAM0 --lam1 $LAM1 --triangle $TRIANGLE --pull_ps $PULL_PS --equil_ps $EQUIL_PS --report_ps $REPORT_PS --kcv $KCV \
     --gb_cutoff_nm 2.0 --platform auto --seed 1 > ~/dna.log 2>&1 & echo DNA_PID \$!"
echo "launched. zone=$ZONE out=$OUTNAME. Poll/teardown with: OUTNAME=$OUTNAME bash $HERE/poll_dna_anim.sh"
