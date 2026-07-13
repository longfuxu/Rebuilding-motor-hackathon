#!/usr/bin/env bash
# ONE Spot A100 session running the heavy driven-MD the Mac's OpenCL can't do fast:
#   Track B (C3-style): targeted MD with the RMSD (difference-to-7JQQ) CV C2 recommends
#     -- RMSDForce is GPU-native/fast on CUDA (it is what cripples the Mac's OpenCL).
#   Track A (S1):        explicit-solvent SMD -- pull dsDNA, first real force number (~<=250k atoms).
#
# *** SAFETY / COORDINATION ***
#   - GPU cap is 1 (GPUS-ALL-REGIONS=1, raise-to-2 was DENIED) -> only ONE GPU VM at a time.
#     This launcher STAGGERS: it waits until the C3 job's "protein-gpu" VM is gone before creating ours.
#   - Uses a DISTINCT VM name "protein-gpu-smd"; teardown deletes ONLY that, NEVER protein-gpu.
#   - Uses 40GB (a100): fits the <=250k-atom explicit SMD AND has capacity (80GB stocks out).
#
# Discipline: Spot GPU, DELETE the VM the instant the run ends, verify 0 instances, report credit
# (~$1.2/hr A100-40GB) + $0 self-pay. Prereq: gcloud auth; A100-40GB Spot quota us-central1.
#
# Usage: bash run_gcp_driven.sh [NSEED_B] [NSEED_A]   (default 2 1)
set -euo pipefail
NSEED_B="${1:-2}"; NSEED_A="${2:-1}"
export GCP_VM="protein-gpu-smd"            # distinct VM -> never collides with the C3 job's protein-gpu
GPU="a100"                                 # 40GB (has capacity; 80GB stocks out)
GC="$HOME/google-cloud-sdk/bin/gcloud"; PROJECT="longfu-protein-gpu"
HERE="$(cd "$(dirname "$0")" && pwd)"; GCPPY="$HERE/../../gcp_pipeline/gcp_gpu.py"
CMIN="$HERE/../openmm_validation/trajectories/C/C_min.pdb"

test -f "$HERE/inputs/C_plus_dna.pdb" || { echo "run build_dna_complex.py first"; exit 1; }
test -f "$CMIN" || { echo "missing C_min.pdb"; exit 1; }

echo "=== [0] STAGGER: wait for the C3 job's protein-gpu VM to free the single GPU slot (cap=1) ==="
for i in $(seq 1 120); do   # up to ~60 min
  running=$($GC compute instances list --project "$PROJECT" --filter='name=protein-gpu' --format='value(name)' 2>/dev/null || true)
  [ -z "$running" ] && { echo "GPU slot free (protein-gpu gone)."; break; }
  echo "  [$i] protein-gpu still up (C3-SBM) — waiting 30s (never touching it)..."; sleep 30
done

echo "=== [1] bring up Spot $GPU as $GCP_VM (zone fallback) ==="
for Z in "${GCP_ZONE:-us-central1-a}" us-central1-b us-central1-c us-central1-f; do
  if GCP_ZONE="$Z" python3 "$GCPPY" up --gpu "$GPU"; then ZONE="$Z"; break; fi
  echo "  $Z stocked out, next..."
done
: "${ZONE:?all Spot zones stocked out}"
cleanup(){ echo "=== teardown (delete ONLY $GCP_VM) ===";
  GCP_VM="$GCP_VM" python3 "$GCPPY" down || \
    $GC compute instances delete "$GCP_VM" --project "$PROJECT" --zone "$ZONE" --quiet || true
  echo "--- instances now ---"; $GC compute instances list --project "$PROJECT" || true; }
trap cleanup EXIT INT TERM
S(){ $GC compute ssh "$GCP_VM" --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap --command "$1"; }
CP(){ $GC compute scp --project "$PROJECT" --zone "$ZONE" --tunnel-through-iap --recurse "$@"; }

echo "=== [2] upload scripts + inputs ==="
S 'mkdir -p ~/in/inputs ~/out'
CP "$HERE/tmd_openclose.py" "$HERE/tmd_staircase.py" "$HERE/smd_pull_dna.py" "$GCP_VM:~/in/"
CP "$CMIN" "$HERE/inputs/C_plus_dna_relaxed.pdb" "$GCP_VM:~/in/inputs/"   # relaxed S0 (raw rigid complex has a 0.43A overlap)

echo "=== [3] install CUDA OpenMM in a venv (image's system python3 can't pip-import openmm; venv works -- proven by the C3 SBM job) ==="
S "sudo apt-get -qq update >/dev/null 2>&1; sudo apt-get -qq install -y python3.10-venv python3-pip >/dev/null 2>&1; \
   rm -rf ~/venv; python3 -m venv ~/venv && ~/venv/bin/pip -q install --upgrade pip 2>&1 | tail -1 && \
   ~/venv/bin/pip -q install openmm pdbfixer mdtraj 2>&1 | tail -2 && \
   ~/venv/bin/python -c 'from openmm import Platform; print(\"OPENMM_OK platforms:\", [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())])'"

echo "=== [3b] wait for the NVIDIA driver so OpenMM's GPU platform (CUDA or GPU-OpenCL) is real, not CPU ==="
# NOTE: on this image the openmm PIP wheel exposes GPU-OpenCL, not CUDA (the C3-SBM job ran on OpenCL).
# So we wait for nvidia-smi (driver loaded), then run with --platform auto (picks CUDA if present, else the
# GPU OpenCL). Track B uses tmd_staircase.py (CustomCentroidBondForce, no RMSDForce) so OpenCL is fast.
S "for i in \$(seq 1 30); do nvidia-smi -L >/dev/null 2>&1 && { echo GPU_READY: \$(nvidia-smi -L | head -1); break; }; echo \"waiting-for-nvidia-driver \$i\"; sleep 20; done; \
   nvidia-smi -L >/dev/null 2>&1 || { echo 'ERROR: NVIDIA driver never came up'; exit 1; }; \
   ~/venv/bin/python -c \"from openmm import Platform; print('OPENMM platforms:', [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())])\""

echo "=== [4] Track B: targeted MD (staircase CV, GPU-native, robust on CUDA or OpenCL), $NSEED_B seeds ==="
S "cd ~/in && for s in \$(seq 1 $NSEED_B); do echo TMD-staircase seed \$s; \
   ~/venv/bin/python tmd_staircase.py --input inputs/C_min.pdb --minimize 0 --equil_ps 20 \
     --pull_ps 200 --kcv 30000 --gb_cutoff_nm 2.0 --report_ps 4 --platform auto --seed \$s --out ~/out/tmdB_stair_s\$s; done"

echo "=== [4b] download Track B NOW (a Track A crash must not lose it) ==="
mkdir -p "$HERE/runs/gcp"
CP "$GCP_VM:~/out/tmdB_stair_s*" "$HERE/runs/gcp/" || true

echo "=== [5] Track A: explicit-solvent SMD, $NSEED_A seed(s), --platform auto ==="
echo "    NOTE: uses the RELAXED S0 (C_plus_dna_relaxed.pdb). Explicit box ~1.16M atoms. A single fast"
echo "    pull yields a DRAG force, NOT the OT stall -- for a 57pN-comparable number use near-quasi-static"
echo "    / multi-seed Jarzynski (much slower). This is a first mechanistic pull, not the stall number."
S "cd ~/in && for s in \$(seq 1 $NSEED_A); do echo SMD-explicit seed \$s; \
   ~/venv/bin/python smd_pull_dna.py --input inputs/C_plus_dna_relaxed.pdb --mode explicit --platform auto \
     --pull_A ${SMD_PULL_A:-20} --pull_ps ${SMD_PULL_PS:-1500} --equil_ps ${SMD_EQUIL_PS:-150} \
     --report_ps 5 --seed \$s --out ~/out/smd_explicit_s\$s; done"

echo "=== [6] download + summarize ==="
mkdir -p "$HERE/runs/gcp"
CP "$GCP_VM:~/out/*" "$HERE/runs/gcp/" || true
python3 - "$HERE/runs/gcp" <<'PY'
import sys, glob, json, os
for f in sorted(glob.glob(os.path.join(sys.argv[1], '*', 'series.json'))):
    r = json.load(open(f))['result']; n = os.path.basename(os.path.dirname(f))
    if 'mean_force_pN' in r:   # SMD (Track A)
        print(f"{n}: <F>={r['mean_force_pN']:.0f}pN plateau={r['plateau_force_pN']:.0f}pN "
              f"W={r['W_total_kcal']:.1f}kcal n_eng {r['n_eng_start']}->{r['n_eng_final']}")
    else:                       # tmd_staircase (Track B)
        print(f"{n}: W={r.get('W_total_kcal',0):.1f}kcal planarity {r.get('planarity_start','?')}->"
              f"{r.get('planarity_final','?')}A axspan {r.get('axial_span_start','?')}->"
              f"{r.get('axial_span_final','?')}A n_eng {r.get('n_eng_start','?')}->{r.get('n_eng_final','?')}")
PY
echo "REMINDER: log credit (~\$1.2/A100-40GB-hr), \$0 self-pay, $GCP_VM deleted above (verify 0 of it)."
