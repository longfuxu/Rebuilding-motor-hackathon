
import os, json, glob, shutil, subprocess, threading, time as _t
os.makedirs("out",exist_ok=True)
def _hb():
    while True: print(f"[hb] {int(_t.time())}",flush=True); _t.sleep(30)
threading.Thread(target=_hb,daemon=True).start()
S=json.load(open("stageC_seqs.json"))
def run(cmd):
    print("CMD"," ".join(cmd),flush=True); r=subprocess.run(cmd,capture_output=True,text=True)
    if r.returncode!=0: print("ERR",r.stderr[-1800:],flush=True)
    return r.returncode
summary={}
for name,seq in S.items():
    d=f"w_{name}"; os.makedirs(d,exist_ok=True)
    open(f"{d}/{name}.fasta","w").write(f">A|protein|\n{seq}\n")
    print("RUN",name,len(seq),flush=True)
    rc=run(["boltz","predict",f"{d}/{name}.fasta","--out_dir",d,"--use_msa_server",
            "--output_format","pdb","--recycling_steps","3","--sampling_steps","200",
            "--diffusion_samples","1","--override","--write_full_pae"])
    print(name,"rc",rc,flush=True)
    for f in glob.glob(f"{d}/**/*.pdb",recursive=True): shutil.copy(f,f"out/{name}__boltz__"+os.path.basename(f))
    for f in glob.glob(f"{d}/**/confidence*.json",recursive=True): shutil.copy(f,f"out/{name}__boltz__"+os.path.basename(f))
    for f in glob.glob(f"{d}/**/*.npz",recursive=True): shutil.copy(f,f"out/{name}__boltz__"+os.path.basename(f))
    confs=glob.glob(f"{d}/**/confidence*.json",recursive=True)
    if confs: summary[name]=json.load(open(confs[0]))
json.dump(summary,open("out/stageC_summary.json","w"),indent=2)
print("STAGE C DONE",list(summary),flush=True)
