import os, sys, json, glob, subprocess, traceback, time
os.chdir("/content")
# Args: progress_json (from a diffusion run), tag for this downstream batch, num_seqs
PROG = sys.argv[1]
TAG = sys.argv[2]
NUM_SEQS = int(sys.argv[3]) if len(sys.argv) > 3 else 8
prog = json.load(open(PROG))
t0 = time.time()
results = {}
for name, info in prog.items():
    if not info.get("ok"):
        results[name] = {"skip": "diffusion failed"}; continue
    pdb0 = f"outputs/{name}_0.pdb"
    if not os.path.isfile(pdb0):
        results[name] = {"skip": "no backbone pdb"}; continue
    contigs = info["contigs"]; copies = info["copies"]
    nd = len(glob.glob(f"outputs/{name}_*.pdb"))
    contig_str = ":".join(contigs)
    opts = [f"--pdb=outputs/{name}_0.pdb", f"--loc=outputs/{name}", f"--contig={contig_str}",
            f"--copies={copies}", f"--num_seqs={NUM_SEQS}", "--num_recycles=1", "--rm_aa=C",
            "--mpnn_sampling_temp=0.1", f"--num_designs={nd}", "--use_soluble"]
    cmd = "python colabdesign/rf/designability_test.py " + " ".join(opts)
    try:
        with open(f"/content/{TAG}.dt.log", "a") as lg:
            lg.write("\n### " + cmd + "\n"); lg.flush()
            p = subprocess.run(cmd, shell=True, stdout=lg, stderr=subprocess.STDOUT, timeout=5400)
        best = f"outputs/{name}/best.pdb"
        info_line = open(best).readline().strip() if os.path.isfile(best) else ""
        results[name] = {"ok": True, "rc": p.returncode, "nd": nd, "best_remark": info_line,
                         "elapsed": round(time.time()-t0, 1)}
    except Exception:
        results[name] = {"ok": False, "err": traceback.format_exc()[-2000:]}
    json.dump(results, open(f"/content/{TAG}.dt.progress.json", "w"), indent=2)
open(f"/content/{TAG}.dt.DONE", "w").write("done\n")
