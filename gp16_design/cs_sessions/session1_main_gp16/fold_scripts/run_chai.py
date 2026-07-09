
import os, json, glob, shutil, torch
from pathlib import Path
from chai_lab.chai1 import run_inference
os.makedirs("out", exist_ok=True)
JOBS=json.load(open("jobs.json"))
summary={}
for name, chains in JOBS.items():
    d=Path(f"work_{name}"); d.mkdir(exist_ok=True)
    fa=d/(name+".fasta")
    with open(fa,"w") as f:
        for cid,seq in chains:
            f.write(f">protein|name={cid}\n{seq}\n")
    print("RUN",name,[c[0] for c in chains],flush=True)
    try:
        out=run_inference(fasta_file=fa, output_dir=d/"pred",
            num_trunk_recycles=3, num_diffn_timesteps=200, seed=42,
            device="cuda", use_esm_embeddings=True)
        # chai writes ranked structures + npz scores
        for f in glob.glob(str(d/"pred/**/*.cif"), recursive=True)+glob.glob(str(d/"pred/**/*.pdb"), recursive=True):
            shutil.copy(f, f"out/{name}__chai__"+os.path.basename(f))
        for f in glob.glob(str(d/"pred/**/*.npz"), recursive=True):
            shutil.copy(f, f"out/{name}__chai__"+os.path.basename(f))
        # aggregate scores if available
        try:
            import numpy as np
            scores=sorted(glob.glob(str(d/"pred/**/*scores*.npz"), recursive=True))
            if scores:
                z=np.load(scores[0]); summary[name]={k:float(np.mean(z[k])) for k in z.files if z[k].dtype.kind=="f"}
        except Exception as e:
            print("score agg err",name,e,flush=True)
        print("CHAI OK",name,flush=True)
    except Exception as e:
        print("CHAI FAIL",name,repr(e)[:500],flush=True)
json.dump(summary, open("out/chai_confidence_summary.json","w"), indent=2)
print("CHAI DONE",list(summary.keys()),flush=True)
