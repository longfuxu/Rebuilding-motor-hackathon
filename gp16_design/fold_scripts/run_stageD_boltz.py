
import os,json,glob,shutil,subprocess,threading,time as _t
os.makedirs("out",exist_ok=True)
threading.Thread(target=lambda:[print(f"[hb]{int(_t.time())}",flush=True) or _t.sleep(30) for _ in iter(int,1)],daemon=True).start()
S=json.load(open("stageD_seq.json")); name,seq=list(S.items())[0]
open(f"{name}.fasta","w").write(f">A|protein|\n{seq}\n")
r=subprocess.run(["boltz","predict",f"{name}.fasta","--out_dir","bz","--use_msa_server",
    "--output_format","pdb","--recycling_steps","3","--sampling_steps","200",
    "--diffusion_samples","5","--override","--write_full_pae"],capture_output=True,text=True)
print("boltz rc",r.returncode,flush=True)
if r.returncode: print(r.stderr[-2000:],flush=True)
for f in glob.glob("bz/**/*.pdb",recursive=True): shutil.copy(f,"out/"+name+"__boltz__"+os.path.basename(f))
for f in glob.glob("bz/**/confidence*.json",recursive=True): shutil.copy(f,"out/"+name+"__boltz__"+os.path.basename(f))
print("BOLTZ D DONE",len(glob.glob("out/*boltz*model*.pdb")),flush=True)
