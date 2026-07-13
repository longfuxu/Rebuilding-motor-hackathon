#!/usr/bin/env bash
# Goal-4 Aim C3 — UMBRELLA SAMPLING of the gp16 single-chain ring along the axial-staircase CV
# to get the P<->H free-energy landscape G(xi). Non-blocking: brings up a Spot A100, uploads the
# driver + helper + input, installs CUDA OpenMM in a venv, and launches tmd_umbrella.py DETACHED
# (nohup) so a 6 h run survives SSH disconnect. Poll + download with poll_umbrella.sh.
#
# *** SAFETY / COORDINATION ***  GPU cap = 1 (raise-to-2 DENIED) -> only ONE GPU VM project-wide.
#   This launcher REFUSES to create a VM if ANY compute instance is already up (never races the cap).
#   VM name = "protein-gpu" (matches GOAL4_HANDOFF s4A + gcp_gpu.py default). Delete ONLY it when done.
# Discipline: Spot, delete the VM the instant the run ends, verify 0 instances, report credit
#   (~$1.2/hr A100-40GB) + $0 self-pay.
#
# Usage: bash run_gcp_umbrella.sh              # production params from GOAL4_HANDOFF s4A
#        NWIN=6 EQUIL_PS=50 SAMPLE_PS=100 BOTH=0 bash run_gcp_umbrella.sh   # short smoke
set -euo pipefail
export GCP_VM="${GCP_VM:-protein-gpu}"
GPU="a100"                                  # 40GB (has capacity; 80GB stocks out)
GC="$HOME/google-cloud-sdk/bin/gcloud"; PROJECT="longfu-protein-gpu"
HERE="$(cd "$(dirname "$0")" && pwd)"; GCPPY="$HERE/../../gcp_pipeline/gcp_gpu.py"
INPUT="$HERE/../openmm_validation/trajectories/C/C_start.pdb"
# umbrella params (GOAL4_HANDOFF s4A defaults; override via env for a smoke test)
NWIN="${NWIN:-24}"; EQUIL_PS="${EQUIL_PS:-300}"; SAMPLE_PS="${SAMPLE_PS:-1200}"
BOTH="${BOTH:-1}"; TIMESTEP_FS="${TIMESTEP_FS:-3}"; SEED="${SEED:-1}"

test -f "$INPUT" || { echo "missing input $INPUT"; exit 1; }
test -f "$HERE/tmd_umbrella.py" || { echo "missing tmd_umbrella.py"; exit 1; }

echo "=== [0] cap=1 guard: refuse to launch if ANY instance is already up ==="
running=$($GC compute instances list --project "$PROJECT" --format='value(name,zone)' 2>/dev/null || true)
if [ -n "$running" ]; then echo "REFUSING: an instance is already up (GPU cap=1):"; echo "$running"; exit 3; fi
echo "GPU slot free."

echo "=== [1] bring up Spot $GPU as $GCP_VM (zone fallback) ==="
ZONE=""
for Z in "${GCP_ZONE:-us-central1-a}" us-central1-b us-central1-c us-central1-f; do
  if GCP_ZONE="$Z" python3 "$GCPPY" up --gpu "$GPU"; then ZONE="$Z"; break; fi
  echo "  $Z stocked out, next..."
done
: "${ZONE:?all Spot zones stocked out}"
echo "$ZONE" > "$HERE/.umbrella_zone"          # remembered for poll/teardown
echo "VM up in $ZONE"

S(){ $GC compute ssh "$GCP_VM" --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap --command "$1"; }
CP(){ $GC compute scp --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap --recurse "$@"; }

echo "=== [2] upload driver + helper + input ==="
S 'mkdir -p ~/in ~/out'
CP "$HERE/tmd_umbrella.py" "$HERE/tmd_staircase.py" "$INPUT" "$GCP_VM:~/in/"

echo "=== [3] install CUDA OpenMM in a venv (cu129 image: apt update FIRST, then venv+pip) ==="
S 'sudo apt-get -qq update >/dev/null 2>&1; sudo apt-get -qq install -y python3.10-venv python3-pip >/dev/null 2>&1; \
   rm -rf ~/be; python3 -m venv ~/be && ~/be/bin/pip -q install --upgrade pip 2>&1 | tail -1 && \
   ~/be/bin/pip -q install openmm numpy mdtraj 2>&1 | tail -2 && \
   ~/be/bin/python -c "from openmm import Platform; print(\"OPENMM_OK platforms:\", [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())])"'

echo "=== [3b] wait for NVIDIA driver (so GPU platform is real, not CPU) ==="
S 'for i in $(seq 1 30); do nvidia-smi -L >/dev/null 2>&1 && { echo GPU_READY: $(nvidia-smi -L | head -1); break; }; echo "waiting-for-nvidia-driver $i"; sleep 20; done; \
   nvidia-smi -L >/dev/null 2>&1 || { echo "ERROR: NVIDIA driver never came up"; exit 1; }'

echo "=== [4] launch umbrella DETACHED (nohup) — nwin=$NWIN equil=$EQUIL_PS sample=$SAMPLE_PS both=$BOTH ==="
S "cd ~/in && rm -f ~/umb.log && nohup ~/be/bin/python tmd_umbrella.py \
     --input C_start.pdb --out ~/out/umb \
     --nwin $NWIN --equil_ps $EQUIL_PS --sample_ps $SAMPLE_PS --both $BOTH \
     --timestep_fs $TIMESTEP_FS --seed $SEED --platform auto \
     > ~/umb.log 2>&1 & echo LAUNCHED_PID \$!"
echo ""
echo "=== launched. VM=$GCP_VM zone=$ZONE ==="
echo "Poll with:  bash $HERE/poll_umbrella.sh"
echo "Teardown:   GCP_ZONE=$ZONE GCP_VM=$GCP_VM python3 $GCPPY down"
