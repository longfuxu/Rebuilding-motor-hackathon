#!/usr/bin/env bash
# P->H->P CYCLE campaign on ONE A100: N forward (P->H) + N reverse (H->P) implicit targeted-MD pulls of the
# single-chain ring + dsDNA (~31k atoms), so we can characterise BOTH half-strokes with a few seeds:
#   - ascent  (loading): fwd_s* = C_plus_dna_relaxed.pdb driven lam 0->1
#   - descent (power)  : rev_s* = same CV geometry, but START from the helical endpoint (--pos_from helical.pdb),
#                        driven lam 1->0  (the stroke that should translocate DNA)
# Reverse is run FIRST (rev_s1) so the new reverse logic is validated on-hardware within the first run.
# Detached; poll + teardown with poll_cycle_campaign.sh. cap=1 guard. VM protein-gpu.
set -euo pipefail
export GCP_VM="${GCP_VM:-protein-gpu}"; GPU="a100"
GC="$HOME/google-cloud-sdk/bin/gcloud"; PROJECT="longfu-protein-gpu"
HERE="$(cd "$(dirname "$0")" && pwd)"; GCPPY="$HERE/../../gcp_pipeline/gcp_gpu.py"
PLANAR="$HERE/inputs/C_plus_dna_relaxed.pdb"
HELICAL="${HELICAL:-$HERE/../../outputs/php_cycle/dna_translocation/dna_anim/final.pdb}"
NSEED="${NSEED:-3}"; PULL_PS="${PULL_PS:-500}"; EQUIL_PS="${EQUIL_PS:-50}"; REPORT_PS="${REPORT_PS:-2}"; KCV="${KCV:-30000}"

test -f "$PLANAR" || { echo "missing $PLANAR"; exit 1; }
test -f "$HELICAL" || { echo "missing helical endpoint $HELICAL"; exit 1; }
echo "=== cap=1 guard ==="
running=$($GC compute instances list --project "$PROJECT" --format='value(name)' 2>/dev/null || true)
[ -n "$running" ] && { echo "REFUSING: instance already up (cap=1): $running"; exit 3; }

echo "=== bring up Spot $GPU as $GCP_VM ==="
ZONE=""
for Z in "${GCP_ZONE:-us-central1-a}" us-central1-b us-central1-c us-central1-f; do
  if GCP_ZONE="$Z" python3 "$GCPPY" up --gpu "$GPU"; then ZONE="$Z"; break; fi
done
: "${ZONE:?all Spot zones stocked out}"; echo "$ZONE" > "$HERE/.campaign_zone"
S(){ $GC compute ssh "$GCP_VM" --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap --command "$1"; }
CP(){ $GC compute scp --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap --recurse "$@"; }

echo "=== upload + install ==="
S 'mkdir -p ~/in/inputs ~/out'
CP "$HERE/tmd_staircase.py" "$GCP_VM:~/in/"
CP "$PLANAR" "$GCP_VM:~/in/inputs/C_plus_dna_relaxed.pdb"
CP "$HELICAL" "$GCP_VM:~/in/inputs/helical.pdb"
S 'sudo apt-get -qq update >/dev/null 2>&1; sudo apt-get -qq install -y python3.10-venv python3-pip >/dev/null 2>&1; \
   test -d ~/be || python3 -m venv ~/be; ~/be/bin/pip -q install openmm numpy mdtraj 2>&1 | tail -1; \
   for i in $(seq 1 30); do nvidia-smi -L >/dev/null 2>&1 && break; sleep 15; done; nvidia-smi -L | head -1'

echo "=== launch detached CYCLE campaign (rev first, then fwd, then rest of rev; NSEED=$NSEED, ${PULL_PS}ps each) ==="
S "cat > ~/campaign.sh <<'EOS'
set -e
PY=~/be/bin/python
COMMON=\"--pull_ps $PULL_PS --equil_ps $EQUIL_PS --report_ps $REPORT_PS --kcv $KCV --gb_cutoff_nm 2.0 --platform auto\"
run_rev(){ \$PY tmd_staircase.py --input inputs/C_plus_dna_relaxed.pdb --pos_from inputs/helical.pdb --lam0 1 --lam1 0 --out ~/out/rev_s\$1 --seed \$1 \$COMMON; }
run_fwd(){ \$PY tmd_staircase.py --input inputs/C_plus_dna_relaxed.pdb --lam0 0 --lam1 1 --out ~/out/fwd_s\$1 --seed \$1 \$COMMON; }
cd ~/in
run_rev 1                       # reverse FIRST -> validates the H->P logic on-hardware early
for s in \$(seq 1 $NSEED); do run_fwd \$s; done
for s in \$(seq 2 $NSEED); do run_rev \$s; done
echo CAMPAIGN_DONE
EOS
cd ~/in && rm -f ~/campaign.log && nohup bash ~/campaign.sh > ~/campaign.log 2>&1 & echo CAMPAIGN_PID \$!"
echo "launched. zone=$ZONE. Poll/teardown: bash $HERE/poll_cycle_campaign.sh"
