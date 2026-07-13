#!/usr/bin/env bash
# Launch the EXPLICIT-solvent SMD campaign (Track A, steps S1/S2 of STEERED_TARGETED_MD_PLAN)
# on a GCP Spot A100. This is the REAL in-silico optical-tweezers run (~250k atoms) that the
# Mac can't do (no CUDA). Local implicit runs are only a proxy.
#
# Discipline (per the project rules): Spot GPU, DELETE the VM the moment the run ends, and
# verify 0 instances after. Report credit spent (~$1.8/hr A100-80GB) and $0 self-pay.
#
# Prereqs: gcloud auth (`gcloud auth login`); A100 Spot quota in us-central1 (already granted).
# Usage:   bash run_gcp_smd.sh [N_SEEDS] [PULL_A] [PULL_PS]
#   N_SEEDS=1  -> S1 (first force number, ~1 A100-hr, ~$1.8)
#   N_SEEDS=12 -> S2 (Jarzynski dG + <F>, ~half a day, ~$15-20 credit)
set -euo pipefail

N_SEEDS="${1:-1}"
PULL_A="${2:-17}"          # ~5 bp
PULL_PS="${3:-500}"
GPU="${GCP_GPU:-a100-80gb}"
GC="$HOME/google-cloud-sdk/bin/gcloud"
PROJECT="longfu-protein-gpu"
VM="protein-gpu-smd"          # DISTINCT name — never touch the C3 job's default protein-gpu VM
export GCP_VM="$VM"
HERE="$(cd "$(dirname "$0")" && pwd)"
GCPPY="$HERE/../../gcp_pipeline/gcp_gpu.py"

echo "=== [S0] verify the threaded complex exists locally ==="
test -f "$HERE/inputs/C_plus_dna.pdb" || { echo "missing inputs/C_plus_dna.pdb — run build_dna_complex.py first"; exit 1; }

echo "=== [1] bring up Spot $GPU (us-central1; GCP_ZONE overrides for stockout fallback) ==="
# gcp_gpu.py up respects GCP_ZONE; scan zones if the default stocks out.
for Z in "${GCP_ZONE:-us-central1-a}" us-central1-b us-central1-c us-central1-f; do
  if GCP_ZONE="$Z" python3 "$GCPPY" up --gpu "$GPU"; then ZONE="$Z"; break; fi
  echo "  zone $Z stocked out, trying next..."
done
: "${ZONE:?all Spot zones stocked out — retry later or use --on-demand where quota exists}"
echo "VM up in $ZONE"

cleanup() {
  echo "=== [X] tearing down VM (Spot; must reach 0 instances) ==="
  python3 "$GCPPY" down || $GC compute instances delete "$VM" --project "$PROJECT" --zone "$ZONE" --quiet || true
  echo "--- instances now (must be empty) ---"
  $GC compute instances list --project "$PROJECT" || true
}
trap cleanup EXIT

echo "=== [2] upload driven/ scripts + threaded complex ==="
$GC compute ssh "$VM" --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap --command 'mkdir -p ~/in/inputs ~/out'
$GC compute scp --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap --recurse \
  "$HERE/smd_pull_dna.py" "$HERE/inputs/C_plus_dna.pdb" "$VM:~/in/"
$GC compute ssh "$VM" --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap \
  --command 'mv ~/in/C_plus_dna.pdb ~/in/inputs/ 2>/dev/null; true'

echo "=== [3] install OpenMM (pip wheel = CUDA build) + run $N_SEEDS explicit SMD pull(s) ==="
$GC compute ssh "$VM" --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap --command "
  set -e
  pip -q install openmm pdbfixer mdtraj 2>/dev/null
  python -c 'from openmm import Platform; print(\"platforms:\", [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())])'
  cd ~/in
  for s in \$(seq 1 $N_SEEDS); do
    echo \"=== explicit SMD seed \$s ===\"
    python smd_pull_dna.py --input inputs/C_plus_dna.pdb --mode explicit \
      --pull_A $PULL_A --pull_ps $PULL_PS --equil_ps 200 --seed \$s \
      --platform CUDA --out ~/out/gcp_explicit_s\$s
  done
"

echo "=== [4] pull results back ==="
mkdir -p "$HERE/runs/gcp_explicit"
$GC compute scp --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap --recurse \
  "$VM:~/out/*" "$HERE/runs/gcp_explicit/" || true

echo "=== [5] summary (force numbers) ==="
python3 - "$HERE/runs/gcp_explicit" <<'PY'
import sys, glob, json, os
for f in sorted(glob.glob(os.path.join(sys.argv[1], '*', 'series.json'))):
    r = json.load(open(f))['result']
    print(f"{os.path.basename(os.path.dirname(f))}: <F>={r['mean_force_pN']:.0f}pN "
          f"plateau={r['plateau_force_pN']:.0f}pN W={r['W_total_kcal']:.1f}kcal/mol "
          f"n_eng {r['n_eng_start']}->{r['n_eng_final']}")
PY
# cleanup() runs on EXIT (deletes VM, verifies 0 instances)
echo "REMINDER: log credit spent (~\$1.8/A100-80GB-hr), confirm \$0 self-pay, VM deleted above."
