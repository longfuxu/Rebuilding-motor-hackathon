#!/usr/bin/env bash
# INTEGRATED helical<->planar + descent-grip model on ONE A100: 2 modes (concerted, sequential) x NSEED seeds
# (replicates → reviewer concern ii). Grip anchored on real DNA-contact residues (concern iv). Ring+dsDNA implicit.
# cap=1 guard; detached; poll+teardown with poll_integrated.sh.
set -euo pipefail
export GCP_VM="${GCP_VM:-protein-gpu}"; GPU="a100"
GC="$HOME/google-cloud-sdk/bin/gcloud"; PROJECT="longfu-protein-gpu"
HERE="$(cd "$(dirname "$0")" && pwd)"; GCPPY="$HERE/../../gcp_pipeline/gcp_gpu.py"
INPUT="$HERE/inputs/C_plus_dna_relaxed.pdb"
NSEED="${NSEED:-2}"; NCYCLES="${NCYCLES:-4}"; CYCLE_PS="${CYCLE_PS:-400}"; KGRIP="${KGRIP:-12000}"; REPORT_PS="${REPORT_PS:-8}"

test -f "$INPUT" || { echo "missing $INPUT"; exit 1; }
echo "=== cap=1 guard ==="
running=$($GC compute instances list --project "$PROJECT" --format='value(name)' 2>/dev/null || true)
[ -n "$running" ] && { echo "REFUSING: instance up (cap=1): $running"; exit 3; }
echo "=== bring up Spot $GPU ==="
ZONE=""
for Z in "${GCP_ZONE:-us-central1-a}" us-central1-b us-central1-c us-central1-f; do
  if GCP_ZONE="$Z" python3 "$GCPPY" up --gpu "$GPU"; then ZONE="$Z"; break; fi
done
: "${ZONE:?stocked out}"; echo "$ZONE" > "$HERE/.integrated_zone"
S(){ $GC compute ssh "$GCP_VM" --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap --command "$1"; }
CP(){ $GC compute scp --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap --recurse "$@"; }
echo "=== upload + install ==="
S 'mkdir -p ~/in/inputs ~/out'
CP "$HERE/tmd_integrated.py" "$HERE/tmd_staircase.py" "$GCP_VM:~/in/"
CP "$INPUT" "$GCP_VM:~/in/inputs/"
S 'sudo apt-get -qq update >/dev/null 2>&1; sudo apt-get -qq install -y python3.10-venv python3-pip >/dev/null 2>&1; \
   test -d ~/be || python3 -m venv ~/be; ~/be/bin/pip -q install openmm numpy mdtraj 2>&1 | tail -1; \
   for i in $(seq 1 30); do nvidia-smi -L >/dev/null 2>&1 && break; sleep 15; done; nvidia-smi -L | head -1'
echo "=== launch detached: 2 modes x $NSEED seeds ==="
S "cat > ~/integ.sh <<'EOS'
set -e; PY=~/be/bin/python; cd ~/in
C=\"--input inputs/C_plus_dna_relaxed.pdb --ncycles $NCYCLES --cycle_ps $CYCLE_PS --kgrip $KGRIP --report_ps $REPORT_PS --equil_ps 50 --platform auto\"
for mode in concerted sequential; do for s in \$(seq 1 $NSEED); do
  \$PY tmd_integrated.py --mode \$mode --seed \$s --out ~/out/integ_\${mode}_s\${s} \$C
done; done
echo INTEGRATED_ALL_DONE
EOS
cd ~/in && rm -f ~/integ.log && nohup bash ~/integ.sh > ~/integ.log 2>&1 & echo INTEGRATED_PID \$!"
echo "launched. zone=$ZONE. Poll/teardown: bash $HERE/poll_integrated.sh"
