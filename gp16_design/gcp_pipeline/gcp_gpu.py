#!/usr/bin/env python3
"""Reusable GCP GPU pipeline — spins up a Spot GPU VM, runs a protein-design task,
downloads the results, and DELETES the VM (cost-safe). Tasks: rfdiffusion, mpnn,
md (OpenMM), boltz (structure prediction, optional MSA).

Why GCP over Colab: SSH is stable (no websocket hangs / idle self-termination).
Cost: Spot GPUs + auto-delete; usage paid by the account's credits first (see README).

Examples:
  python gcp_gpu.py md        --input struct.pdb                 --out ./out --gpu l4 --ns 5
  python gcp_gpu.py rfdiffusion --input motif.pdb --contig 'A1-100/20-30/B1-100' --num 20 --out ./out --gpu l4
  python gcp_gpu.py mpnn      --input backbone.pdb --fixed '10 11 12' --nseq 8 --out ./out --gpu t4
  python gcp_gpu.py boltz     --input seqs.fasta   --msa core.a3m   --out ./out --gpu l4
  python gcp_gpu.py up --gpu a100          # just create a VM (manual use); then `ssh`, then `down`
  python gcp_gpu.py down                   # delete the VM
"""
import argparse, subprocess, os, sys, time

GCLOUD = os.environ.get("GCLOUD", "/Users/longfu/google-cloud-sdk/bin/gcloud")
PROJECT = os.environ.get("GCP_PROJECT", "longfu-protein-gpu")
VM = os.environ.get("GCP_VM", "protein-gpu")
# GPU -> (machine type, accelerator, default zone). L4 = g2; A100 = a2; T4 = n1.
GPUS = {
    "l4":        ("g2-standard-8",  "type=nvidia-l4,count=1",          "us-central1-a"),
    "a100":      ("a2-highgpu-1g",  "type=nvidia-tesla-a100,count=1",  "us-central1-a"),
    "a100-80gb": ("a2-ultragpu-1g", "type=nvidia-a100-80gb,count=1",   "us-central1-a"),
    "t4":        ("n1-standard-8",  "type=nvidia-tesla-t4,count=1",     "us-central1-a"),
}
# CUDA deep-learning image (has CUDA + conda). PyTorch installs per-recipe.
IMAGE = ["--image-family=common-cu129-ubuntu-2204-nvidia-580", "--image-project=deeplearning-platform-release"]


def run(cmd, timeout=None, check=True, quiet=False):
    if not quiet:
        print("+ " + " ".join(cmd[:6]) + (" ..." if len(cmd) > 6 else ""), flush=True)
    r = subprocess.run(cmd, timeout=timeout)
    if check and r.returncode != 0:
        sys.exit(f"command failed ({r.returncode})")
    return r.returncode


def g(*args, **kw):
    return run([GCLOUD] + list(args), **kw)


def zone_of(gpu):
    return os.environ.get("GCP_ZONE", GPUS[gpu][2])


def vm_up(gpu, spot=True, keep_note=True):
    mt, accel, _z = GPUS[gpu]
    zone = zone_of(gpu)  # respect GCP_ZONE override (capacity/stockout fallback across zones)
    args = ["compute", "instances", "create", VM, "--project", PROJECT, "--zone", zone,
            "--machine-type", mt, "--accelerator", accel, "--maintenance-policy", "TERMINATE",
            "--boot-disk-size", "200GB", "--metadata", "install-nvidia-driver=True"] + IMAGE
    if spot:
        args += ["--provisioning-model=SPOT", "--instance-termination-action=DELETE"]
    print(f"creating {gpu.upper()} VM '{VM}' in {zone} (spot={spot}) — DELETE it when done to stop cost", flush=True)
    g(*args)
    # wait for SSH
    for _ in range(30):
        if g("compute", "ssh", VM, "--project", PROJECT, "--zone", zone, "--command=echo ready",
             "--tunnel-through-iap", check=False, quiet=True) == 0:
            print("VM reachable via SSH", flush=True); return zone
        time.sleep(10)
    sys.exit("VM created but SSH not reachable")


def vm_down():
    for gpu in GPUS:
        g("compute", "instances", "delete", VM, "--project", PROJECT, "--zone", zone_of(gpu),
          "-q", check=False, quiet=True)
    print(f"deleted VM '{VM}' (if it existed)", flush=True)


def ssh(zone, script):
    return g("compute", "ssh", VM, "--project", PROJECT, "--zone", zone, "--tunnel-through-iap",
             "--command", script)


def scp_up(zone, local, remote):
    g("compute", "scp", "--project", PROJECT, "--zone", zone, "--tunnel-through-iap",
      "--recurse", local, f"{VM}:{remote}")


def scp_down(zone, remote, local):
    g("compute", "scp", "--project", PROJECT, "--zone", zone, "--tunnel-through-iap",
      "--recurse", f"{VM}:{remote}", local, check=False)


# ---- recipes: (install_script, run_script_template). {inp}=remote input basename, {out}=/home/out ----
RECIPES = {
    "rfdiffusion": {
        "install": (
            "set -e; test -d ~/RFdiffusion || (git clone -q https://github.com/sokrypton/RFdiffusion.git; "
            "pip -q install jedi omegaconf hydra-core icecream pyrsistent pynvml decorator dllogger 2>/dev/null; "
            "pip -q install --no-dependencies dgl -f https://data.dgl.ai/wheels/torch-2.4/cu124/repo.html 2>/dev/null; "
            "pip -q install --no-dependencies e3nn==0.5.5 opt_einsum_fx 2>/dev/null; "
            "cd ~/RFdiffusion/env/SE3Transformer && pip -q install . ); "
            "test -d ~/RFdiffusion/models || (mkdir -p ~/RFdiffusion/models; cd ~/RFdiffusion/models; "
            "for m in Base_ckpt Complex_base_ckpt; do wget -qnc http://files.ipd.uw.edu/pub/RFdiffusion/$( [ $m = Base_ckpt ] && echo 6f5902ac237024bdd0c176cb93063dc4 || echo e29311f6f1bf1af907f9ef9f44b8328b )/$m.pt; done); "
            "cd ~/RFdiffusion && wget -qnc https://files.ipd.uw.edu/krypton/schedules.zip && unzip -qo schedules.zip 2>/dev/null || true"),
        # run: uses RFdiffusion run_inference directly with a motif + contig
        "run": ("cd ~/RFdiffusion && DGLBACKEND=pytorch python scripts/run_inference.py "
                "inference.input_pdb=~/in/{inp} inference.output_prefix=~/out/design "
                "'contigmap.contigs=[{contig}]' inference.num_designs={num} diffuser.T=50"),
    },
    "mpnn": {
        "install": "test -d ~/ProteinMPNN || git clone -q https://github.com/dauparas/ProteinMPNN.git",
        "run": ("cd ~/ProteinMPNN && python protein_mpnn_run.py --pdb_path ~/in/{inp} "
                "--out_folder ~/out --num_seq_per_target {nseq} --sampling_temp 0.1 "
                "--use_soluble_model --model_name v_48_020 --omit_AAs C {fixed_arg}"),
    },
    "md": {
        "install": "pip -q install openmm pdbfixer mdtraj 2>/dev/null",
        # md_run.py is uploaded alongside the input
        "run": "cd ~/in && python md_run.py ~/in/{inp} run {nsteps} && mv ~/in/*.md.json ~/in/*.dcd ~/in/*_start.pdb ~/out/ 2>/dev/null; true",
    },
    "boltz": {
        "install": "pip -q install boltz 2>/dev/null",
        # boltz predict with optional MSA (enables tiled-MSA folding, unlike single-seq NIM)
        "run": "cd ~/in && boltz predict {inp} --out_dir ~/out --use_msa_server {msa_arg}",
    },
    # BioEmu: sequence->equilibrium ensemble of a MONOMER. Pass a TILED a3m as {inp}
    # (first row = query) so the concatenated single-chain repeats keep clean per-repeat
    # co-evolution (default MMseqs2 mis-aligns them). Large rings need A100 memory.
    # WORKING VERSION SET (verified 2026-07-10 on cu129 image): a clean venv + numpy<2
    # + tensorflow-cpu==2.18.0 (MUST match the jax==0.4.35 era; TF 2.21 -> protobuf/abseil
    # double-registration abort). For long sequences pass --batch_size_100 350 (else the
    # auto batch_size = int(bs100*(100/L)^2) underflows to 0 -> "range() arg 3 must not be
    # zero"). NOTE: long-running (~L^2); drive detached (nohup) + poll for >600 aa.
    "bioemu": {
        "install": ("sudo apt-get -qq install -y python3.10-venv >/dev/null 2>&1; "
                    "python3 -m venv ~/be && ~/be/bin/pip -q install --upgrade pip && "
                    "~/be/bin/pip -q install bioemu 'numpy<2' 'tensorflow-cpu==2.18.0'"),
        "run": ("~/be/bin/python -m bioemu.sample --sequence ~/in/{inp} --num_samples {num} "
                "--batch_size_100 350 --output_dir ~/out/bioemu"),
    },
}


def run_task(task, args):
    if task not in RECIPES:
        sys.exit(f"unknown task {task}; choose from {list(RECIPES)}")
    zone = vm_up(args.gpu, spot=not args.on_demand)
    try:
        ssh(zone, "mkdir -p ~/in ~/out")
        scp_up(zone, args.input, "~/in/")
        for extra in (args.upload or []):
            scp_up(zone, extra, "~/in/")
        inp = os.path.basename(args.input)
        rec = RECIPES[task]
        fmt = dict(inp=inp, out="~/out", num=getattr(args, "num", 1), nseq=getattr(args, "nseq", 8),
                   nsteps=getattr(args, "nsteps", 100000), contig=getattr(args, "contig", ""),
                   fixed_arg=(f"--fixed_positions_jsonl {args.fixed}" if getattr(args, "fixed", None) else ""),
                   msa_arg=(f"--msa {os.path.basename(args.msa)}" if getattr(args, "msa", None) else ""))
        print("=== installing " + task + " on VM ===", flush=True)
        ssh(zone, rec["install"])
        print("=== running " + task + " ===", flush=True)
        ssh(zone, rec["run"].format(**fmt))
        os.makedirs(args.out, exist_ok=True)
        scp_down(zone, "~/out/*", args.out)
        print(f"=== results in {args.out} ===", flush=True)
    finally:
        if not args.keep:
            vm_down()
        else:
            print(f"--keep set: VM '{VM}' left running (delete with: python gcp_gpu.py down)", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("task", choices=list(RECIPES) + ["up", "down", "ssh"])
    ap.add_argument("--input"); ap.add_argument("--out", default="./gcp_out")
    ap.add_argument("--gpu", choices=list(GPUS), default="l4")
    ap.add_argument("--on-demand", action="store_true", help="use on-demand instead of Spot (pricier, less preemption)")
    ap.add_argument("--keep", action="store_true", help="don't delete the VM after the run (remember to `down`!)")
    ap.add_argument("--upload", nargs="*", help="extra local files to upload to ~/in (e.g. md_run.py, msa.a3m)")
    ap.add_argument("--num", type=int, default=8); ap.add_argument("--nseq", type=int, default=8)
    ap.add_argument("--nsteps", type=int, default=100000); ap.add_argument("--contig", default="")
    ap.add_argument("--fixed"); ap.add_argument("--msa")
    a = ap.parse_args()
    if a.task == "down": vm_down(); return
    if a.task == "up": vm_up(a.gpu, spot=not a.on_demand); print("VM up; ssh with: gcloud compute ssh " + VM); return
    if a.task == "ssh": ssh(zone_of(a.gpu), "bash -l"); return
    if not a.input: sys.exit("--input required for task " + a.task)
    run_task(a.task, a)


if __name__ == "__main__":
    main()
