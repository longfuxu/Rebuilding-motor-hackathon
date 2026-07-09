
import os,json,glob,shutil,subprocess,threading,time as _t,sys
os.makedirs("out",exist_ok=True)
threading.Thread(target=lambda:[print(f"[hb]{int(_t.time())}",flush=True) or _t.sleep(30) for _ in iter(int,1)],daemon=True).start()
SPEC=json.load(open("job_spec.json"))
seeds=SPEC.get("seeds",[1]); rec=SPEC.get("recycling",3); samp=SPEC.get("sampling",200); ds=SPEC.get("diffusion_samples",1)
def run(cmd):
    print("CMD"," ".join(cmd),flush=True); r=subprocess.run(cmd,capture_output=True,text=True)
    if r.returncode: print("ERR",r.stderr[-1800:],flush=True)
    return r.returncode
summary={}
for name,spec in SPEC["constructs"].items():
    chains=spec["chains"]
    for seed in seeds:
        d=f"w_{name}_s{seed}"; os.makedirs(d,exist_ok=True)
        fa=f"{d}/{name}.fasta"
        with open(fa,"w") as fo:
            for cid,seq in chains.items(): fo.write(f">{cid}|protein|\n{seq}\n")
        print("RUN",name,"seed",seed,"chains",list(chains),flush=True)
        rc=run(["boltz","predict",fa,"--out_dir",d,"--use_msa_server","--output_format","pdb",
                "--recycling_steps",str(rec),"--sampling_steps",str(samp),
                "--diffusion_samples",str(ds),"--seed",str(seed),"--override","--write_full_pae"])
        for f in glob.glob(f"{d}/**/*.pdb",recursive=True): shutil.copy(f,f"out/{name}__boltz__s{seed}__"+os.path.basename(f))
        for f in glob.glob(f"{d}/**/confidence*.json",recursive=True): shutil.copy(f,f"out/{name}__boltz__s{seed}__"+os.path.basename(f))
        confs=glob.glob(f"{d}/**/confidence*.json",recursive=True)
        if confs: summary[f"{name}_s{seed}"]=json.load(open(confs[0]))
        json.dump(summary,open("out/summary.json","w"),indent=2)  # save incrementally
print("JOB DONE",list(summary),flush=True)
