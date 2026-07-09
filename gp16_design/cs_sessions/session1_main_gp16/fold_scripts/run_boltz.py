
import os, subprocess, glob, shutil, json, sys
os.makedirs("out", exist_ok=True)
JOBS = json.load(open("jobs.json"))
summary={}
for name, chains in JOBS.items():
    d=f"work_{name}"; os.makedirs(d, exist_ok=True)
    fa=os.path.join(d, name+".fasta")
    with open(fa,"w") as f:
        for cid, seq in chains:
            f.write(f">{cid}|protein|\n{seq}\n")
    cmd=["boltz","predict",fa,"--out_dir",d,"--use_msa_server","--override",
         "--output_format","pdb"]
    print("RUN", name, [c[0] for c in chains], "lens", [len(c[1]) for c in chains], flush=True)
    rc=subprocess.run(cmd, capture_output=True, text=True)
    print(name,"boltz rc",rc.returncode, flush=True)
    if rc.returncode!=0:
        print("STDERR tail:\n", rc.stderr[-2500:], flush=True)
    # collect predicted structures + confidence
    for f in glob.glob(f"{d}/**/*.pdb", recursive=True)+glob.glob(f"{d}/**/*.cif", recursive=True):
        shutil.copy(f, f"out/{name}__boltz__"+os.path.basename(f))
    for f in glob.glob(f"{d}/**/confidence*.json", recursive=True)+glob.glob(f"{d}/**/*plddt*.npz", recursive=True)+glob.glob(f"{d}/**/pae*.npz", recursive=True):
        shutil.copy(f, f"out/{name}__boltz__"+os.path.basename(f))
    confs=glob.glob(f"{d}/**/confidence*.json", recursive=True)
    if confs: summary[name]=json.load(open(confs[0]))
json.dump(summary, open("out/boltz_confidence_summary.json","w"), indent=2)
print("BOLTZ DONE", list(summary.keys()), flush=True)
