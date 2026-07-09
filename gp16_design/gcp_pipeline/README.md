# gcp_pipeline — reusable GCP GPU tool (RFdiffusion / ProteinMPNN / OpenMM MD / Boltz)

One tool to run protein-design GPU jobs on Google Cloud instead of flaky Colab.
It **spins up a Spot GPU VM → installs the tool → runs your job → downloads results → DELETES the VM**.

## Why GCP (vs Colab)
SSH is stable — no websocket hangs, no idle self-termination. You control the VM.

## Prereqs (already set up 2026-07-09; see system memory `gcp-gpu-account`)
- `gcloud` at `~/google-cloud-sdk/bin/gcloud`, authed as longfu2.xu@gmail.com.
- Project `longfu-protein-gpu` on the **credit** billing account (`0161E2-E9BA4A-8F5A6A`, "Google Gemeni").
- GPU quota: L4 & T4 (us-central1) granted; **A100 requested (pending approval)**.
- First SSH may need: `gcloud compute config-ssh` (keys) — the tool uses `--tunnel-through-iap` (no public IP needed).

## Usage
```bash
GC=gp16_design/gcp_pipeline/gcp_gpu.py
# OpenMM MD (uploads md_run.py automatically-ish; pass it with --upload)
python $GC md --input struct.pdb --upload gp16_design/gcp_pipeline/md_run.py --nsteps 100000 --gpu l4 --out ./out
# RFdiffusion de-novo connector
python $GC rfdiffusion --input motif.pdb --contig 'A1-330/20-30/B1-330' --num 20 --gpu l4 --out ./out
# ProteinMPNN
python $GC mpnn --input backbone.pdb --nseq 8 --gpu t4 --out ./out
# Boltz-2 structure prediction WITH an MSA (fixes the single-seq confound; pass a tiled a3m)
python $GC boltz --input seqs.fasta --msa core.a3m --upload core.a3m --gpu l4 --out ./out
# Manual VM control
python $GC up --gpu a100      # create; python $GC ssh --gpu a100 ; python $GC down   # DELETE
```
GPUs: `--gpu l4|a100|t4` (A100 once quota approves). Spot by default (`--on-demand` to disable). `--keep` leaves the VM up (then run `down` yourself).

## 💰 Cost — no personal money (the important part)
- **Google Cloud applies your CREDITS FIRST, automatically** — every dollar of usage is paid from the ~$45/month Gen-AI-&-Cloud credits (you have several months' worth, >$135, + a $500 cert bonus) **before any card is ever charged**.
- Our jobs are a few GPU-hours (~$1–5 each on Spot). Nowhere near exhausting the credits.
- **The `gp16-gpu-guard` budget ($135) is only an ALERT (email), not a cap** — it doesn't stop or limit spending. It's set near your credit so you'd be warned only if usage ever approached credit exhaustion.
- **The only way to pay out-of-pocket** = spend MORE than ALL your credits in a period. Won't happen with our usage + the two rules below.
- **Two rules that keep it safe (built into this tool):** (1) **Spot GPUs** (~3× cheaper), (2) **the VM is DELETED after every run** — an idle GPU is the only real way to waste credit (same lesson as Colab's idle A100 burn).

## Tasks / recipes
Each task in `gcp_gpu.py` `RECIPES` has an install + run step: `rfdiffusion` (ColabDesign fork + weights), `mpnn` (dauparas/ProteinMPNN), `md` (OpenMM+pdbfixer+mdtraj, runs `md_run.py`), `boltz` (`pip install boltz`, supports `--msa`). Add new tools by adding a recipe.
