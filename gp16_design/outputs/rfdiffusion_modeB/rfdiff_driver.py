import subprocess, os, json, traceback, sys, time
os.chdir("/content")
JOBSPEC = sys.argv[1]
spec = json.load(open(JOBSPEC))
tag = spec["tag"]
t0 = time.time()

# Singleton guard: only one driver per tag may run (prevents double-launch races
# from flaky colab-exec retries corrupting shared output files).
LOCK = f"/content/{tag}.lock"
try:
    _fd = os.open(LOCK, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    os.write(_fd, str(os.getpid()).encode()); os.close(_fd)
except FileExistsError:
    print(f"[{tag}] another driver holds {LOCK}; exiting"); sys.exit(0)
import atexit
atexit.register(lambda: (os.path.isfile(LOCK) and os.remove(LOCK)))

# Re-run setup in a fresh namespace: install dirs already exist, so this only
# re-imports and re-defines run_diffusion (fast). No reinstall.
ns = {}
exec(open("/content/rfdiff_setup.py").read(), ns)

# Headless replacement for the notebook's widget-driven run(): just launch
# run_inference.py and block until it finishes (all num_designs done in one call).
def headless_run(command, steps, num_designs=1, visual="none"):
    with open(f"/content/{tag}.rfinf.log", "a") as lg:
        lg.write("\n### CMD: " + command + "\n"); lg.flush()
        p = subprocess.Popen(command, shell=True, stdout=lg, stderr=subprocess.STDOUT)
        p.wait()
        lg.write(f"### run_inference rc={p.returncode}\n"); lg.flush()
ns["run"] = headless_run
run_diffusion = ns["run_diffusion"]

results = {}
for job in spec["jobs"]:
    name, contig, nd, pT, pdb, chains = job
    js = {"contig": contig, "num_designs": nd, "partial_T": pT, "started": round(time.time()-t0,1)}
    try:
        c, cp = run_diffusion(contigs=contig, pdb=pdb, path=name,
                              iterations=spec.get("iterations", 50), symmetry="none", order=1,
                              hotspot="", chains=chains, add_potential=False, partial_T=pT,
                              num_designs=nd, use_beta_model=spec.get("use_beta_model", False),
                              visual="none")
        import glob
        pdbs = sorted(glob.glob(f"outputs/{name}_*.pdb"))
        js.update(ok=True, contigs=c, copies=cp, n_pdbs=len(pdbs),
                  elapsed=round(time.time()-t0, 1))
    except Exception:
        js.update(ok=False, err=traceback.format_exc()[-2500:], elapsed=round(time.time()-t0, 1))
    results[name] = js
    json.dump(results, open(f"/content/{tag}.progress.json", "w"), indent=2)

open(f"/content/{tag}.DONE", "w").write(f"done in {round(time.time()-t0,1)}s\n")
