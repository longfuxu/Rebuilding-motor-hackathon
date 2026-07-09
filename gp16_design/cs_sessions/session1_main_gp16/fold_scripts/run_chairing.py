
import os,json,glob,shutil,threading,time as _t
from pathlib import Path
threading.Thread(target=lambda:[print(f"[hb]{int(_t.time())}",flush=True) or _t.sleep(30) for _ in iter(int,1)],daemon=True).start()
from chai_lab.chai1 import run_inference
SPEC=json.load(open("job_spec.json")); seeds=SPEC.get("seeds",[42])
os.makedirs("out",exist_ok=True)
for name,spec in SPEC["constructs"].items():
    chains=spec["chains"]
    for seed in seeds:
        d=Path(f"c_{name}_s{seed}"); d.mkdir(exist_ok=True)
        with open(d/"in.fasta","w") as fo:
            for cid,seq in chains.items(): fo.write(f">protein|name={cid}\n{seq}\n")
        try:
            run_inference(fasta_file=d/"in.fasta",output_dir=d/"pred",num_trunk_recycles=3,
                num_diffn_timesteps=200,seed=seed,device="cuda",use_esm_embeddings=True)
            for f in glob.glob(str(d/"pred/**/*.cif"),recursive=True): shutil.copy(f,f"out/{name}__chai__s{seed}__"+os.path.basename(f))
            print(name,"seed",seed,"ok",flush=True)
        except Exception as e:
            print(name,"seed",seed,"FAIL",repr(e)[:300],flush=True)
print("CHAI RING DONE",flush=True)
