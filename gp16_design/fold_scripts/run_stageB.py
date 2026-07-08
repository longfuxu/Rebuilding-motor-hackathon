
import os, json, glob, shutil, subprocess, sys
os.makedirs("out", exist_ok=True)

import threading, time as _t
def _hb():
    while True:
        print(f"[heartbeat] alive t={int(_t.time())}",flush=True); _t.sleep(30)
threading.Thread(target=_hb,daemon=True).start()

V=json.load(open("stageB_variants.json"))
NATIVE=V["native_monomer"]["seq"]

def run(cmd):
    print("CMD:"," ".join(cmd),flush=True)
    r=subprocess.run(cmd,capture_output=True,text=True)
    if r.returncode!=0: print("STDERR:",r.stderr[-2000:],flush=True)
    return r.returncode

# ---- 1. fold native monomer WITH msa server; capture a3m ----
os.makedirs("nat",exist_ok=True)
open("nat/native.fasta","w").write(f">A|protein|\n{NATIVE}\n")
rc=run(["boltz","predict","nat/native.fasta","--out_dir","nat","--use_msa_server",
        "--output_format","pdb","--recycling_steps","3","--sampling_steps","200",
        "--diffusion_samples","1","--override","--write_full_pae"])
print("native fold rc",rc,flush=True)
# find the a3m boltz generated
a3ms=glob.glob("nat/**/*.a3m",recursive=True)+glob.glob("nat/**/msa/*.a3m",recursive=True)
print("a3m found:",a3ms,flush=True)
native_a3m=None
for f in a3ms:
    txt=open(f).read()
    if txt.count(">")>=1: native_a3m=f; break
assert native_a3m, "no native a3m"
# parse a3m -> match-state strings (drop lowercase insertions); index by native resnum
lines=[l.rstrip("\n") for l in open(native_a3m) if l.strip()]
recs=[]  # (header, matchstr)
i=0
while i<len(lines):
    if lines[i].startswith(">"):
        hdr=lines[i]; seq=lines[i+1]; i+=2
        match="".join(ch for ch in seq if ch=="-" or ch.isupper())
        recs.append((hdr,match))
    else: i+=1
qlen=len(recs[0][1])
print("MSA records:",len(recs),"| query match len:",qlen,flush=True)
# collect Boltz native outputs
for f in glob.glob("nat/**/*.pdb",recursive=True):
    shutil.copy(f,"out/native_monomer__boltz__"+os.path.basename(f))
for f in glob.glob("nat/**/confidence*.json",recursive=True):
    shutil.copy(f,"out/native_monomer__boltz__"+os.path.basename(f))

# ---- 2. build permuted a3m per CP variant, fold with cached MSA, reduced settings ----
summary={"native_monomer":{"a3m":os.path.basename(native_a3m),"msa_depth":len(recs)}}
for name,v in V.items():
    if name=="native_monomer": continue
    posmap=v["posmap"]; seq=v["seq"]
    d=f"work_{name}"; os.makedirs(d,exist_ok=True)
    # permuted a3m: query = variant seq; others = permuted match chars ('-' at linker)
    out_a3m=os.path.join(d,name+".a3m")
    with open(out_a3m,"w") as fo:
        fo.write(f">query\n{seq}\n")
        for hdr,match in recs[1:]:
            row="".join(match[r-1] if (r is not None and 1<=r<=qlen) else "-" for r in posmap)
            # skip all-gap rows
            if row.count("-")<len(row):
                fo.write(f"{hdr}\n{row}\n")
    # boltz fasta referencing precomputed msa
    fa=os.path.join(d,name+".fasta")
    open(fa,"w").write(f">A|protein|{os.path.abspath(out_a3m)}\n{seq}\n")
    print(f"RUN {name} cut={v['cut']} (cached MSA depth {len(recs)}, reduced settings)",flush=True)
    rc=run(["boltz","predict",fa,"--out_dir",d,"--output_format","pdb",
            "--recycling_steps","1","--sampling_steps","50","--diffusion_samples","1",
            "--override"])
    print(name,"rc",rc,flush=True)
    for f in glob.glob(f"{d}/**/*.pdb",recursive=True):
        shutil.copy(f,f"out/{name}__boltz__"+os.path.basename(f))
    confs=glob.glob(f"{d}/**/confidence*.json",recursive=True)
    for f in confs: shutil.copy(f,f"out/{name}__boltz__"+os.path.basename(f))
    if confs: summary[name]=json.load(open(confs[0]))
json.dump(summary,open("out/stageB_summary.json","w"),indent=2)
print("STAGE B DONE",list(summary.keys()),flush=True)
