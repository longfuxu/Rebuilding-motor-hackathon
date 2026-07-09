
import os,json,glob,shutil,threading,time as _t
from pathlib import Path
threading.Thread(target=lambda:[print(f"[hb]{int(_t.time())}",flush=True) or _t.sleep(30) for _ in iter(int,1)],daemon=True).start()
from chai_lab.chai1 import run_inference
S=json.load(open("stageD_seq.json")); name,seq=list(S.items())[0]
os.makedirs("out",exist_ok=True)
for seed in [42,7]:
    d=Path(f"chai_{seed}"); d.mkdir(exist_ok=True)
    open(d/"in.fasta","w").write(f">protein|name=A\n{seq}\n")
    try:
        run_inference(fasta_file=d/"in.fasta",output_dir=d/"pred",num_trunk_recycles=3,
            num_diffn_timesteps=200,seed=seed,device="cuda",use_esm_embeddings=True)
        for f in glob.glob(str(d/"pred/**/*.cif"),recursive=True): shutil.copy(f,f"out/{name}__chai__seed{seed}__"+os.path.basename(f))
        for f in glob.glob(str(d/"pred/**/*.npz"),recursive=True): shutil.copy(f,f"out/{name}__chai__seed{seed}__"+os.path.basename(f))
        print("chai seed",seed,"ok",flush=True)
    except Exception as e:
        print("chai seed",seed,"FAIL",repr(e)[:400],flush=True)
print("CHAI D DONE",flush=True)
