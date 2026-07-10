#!/usr/bin/env python3
"""ProteinMPNN-design the connector of each Rho salami backbone (motor motifs frozen).

Backbone layout (contig A175-414 / connector / B175-414): output chain A residues 1..N,
motorA = first 230 modeled native residues, connector = middle, motorB = last 230.
We thread the native Rho motor sequence onto the two motifs, fix them, and design ONLY
the connector, then emit full designed dimer sequences (motorA + connector + motorB).
Catalytic check: R366 and Walker-A K184 must be preserved in the two motif copies.
"""
import os, sys, glob, subprocess, json
import numpy as np
MPNN="/Users/longfu/.claude/jobs/48b7a543/tmp/ProteinMPNN"
PY="/Users/longfu/miniforge3/bin/python3.13"
BASE="/Users/longfu/Developer/claude-science-hackthon/.claude/worktrees/export-gp16-results/gp16_design/outputs/rho_generative"
BB=f"{BASE}/rfdiffusion"
WORK=f"{BASE}/mpnn_work"; os.makedirs(WORK, exist_ok=True)

# native Rho motor res175-414 sequence (== monomer_ref in the manifest, 240 aa)
MAN=json.load(open("/Users/longfu/Developer/claude-science-hackthon/.claude/worktrees/export-gp16-results/gp16_design/pipelines/tiled_msa_fold/manifests/rho_direct_tiledMSA.json"))
RHO_MOTOR="".join(MAN["monomer_ref"].split())   # native res175..414
# the RFdiffusion contig kept these native residues (loops dropped); we thread by ORDINAL onto modeled motif residues
A_segs=["175-185","187-204","206-218","220-244","246-326","328-340","342-379","381-389","391-395","397-404","406-414"]
def native_list(segs):
    out=[]
    for s in segs:
        a,b=s.split('-'); out+=list(range(int(a),int(b)+1))
    return out
MOTIF_NATIVE=native_list(A_segs)          # native residue numbers present in each motif, in order
NA=len(MOTIF_NATIVE)
AA1TO3={"A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","Q":"GLN","E":"GLU","G":"GLY","H":"HIS",
        "I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE","P":"PRO","S":"SER","T":"THR","W":"TRP",
        "Y":"TYR","V":"VAL"}
# native residue -> aa (offset from 175)
def nat_aa(rn): return RHO_MOTOR[rn-175]

def read_ca_res(path):
    rns=[]
    for L in open(path):
        if L[:4]=="ATOM" and L[12:16].strip()=="CA":
            rns.append(int(L[22:26]))
    return rns

def thread(pdb_in, pdb_out):
    """Relabel the first NA and last NA modeled residues to native Rho motor names (by ordinal);
    connector residues left unchanged. Returns (N_modeled_output_max, connector_output_resnums)."""
    rns=read_ca_res(pdb_in)
    # modeled residues in order; motifA = first NA, motifB = last NA
    motifA_out=rns[:NA]; motifB_out=rns[-NA:]
    conn_out=rns[NA:len(rns)-NA]
    a_map={out:MOTIF_NATIVE[i] for i,out in enumerate(motifA_out)}
    b_map={out:MOTIF_NATIVE[i] for i,out in enumerate(motifB_out)}
    out=[]
    for L in open(pdb_in):
        if L[:4]=="ATOM":
            rn=int(L[22:26]); nm=None
            if rn in a_map: nm=AA1TO3[nat_aa(a_map[rn])]
            elif rn in b_map: nm=AA1TO3[nat_aa(b_map[rn])]
            if nm: L=L[:17]+nm+L[20:]
        out.append(L)
    open(pdb_out,"w").write("".join(out))
    return conn_out, motifA_out, motifB_out

def run(cmd,cwd=None):
    env=dict(os.environ); env['PYTHONPATH']=MPNN
    r=subprocess.run(cmd,capture_output=True,text=True,cwd=cwd,env=env)
    if r.returncode!=0: print("ERR",cmd[1].split('/')[-1], r.stderr[-300:])
    return r

def main():
    open(f"{WORK}/all_rho_designs.fasta","w").close()
    results={}
    for pdb in sorted(glob.glob(f"{BB}/rhosal_*.pdb")):
        tag=os.path.splitext(os.path.basename(pdb))[0]
        d=f"{WORK}/{tag}"; os.makedirs(d,exist_ok=True)
        threaded=f"{d}/threaded.pdb"
        conn_out,mA,mB=thread(pdb,threaded)
        run([PY,f"{MPNN}/helper_scripts/parse_multiple_chains.py","--input_path",d,"--output_path",f"{d}/parsed.jsonl"])
        # design ONLY connector positions (chain A), everything else fixed. positions are 1-based indices along chain A residues as ordered
        rns=read_ca_res(threaded)
        pos_index=[str(i+1) for i,rn in enumerate(rns) if rn in set(conn_out)]  # 1-based ordinal along chain
        run([PY,f"{MPNN}/helper_scripts/make_fixed_positions_dict.py","--input_path",f"{d}/parsed.jsonl",
             "--output_path",f"{d}/fixed.jsonl","--chain_list","A","--position_list"," ".join(pos_index),"--specify_non_fixed"])
        run([PY,f"{MPNN}/protein_mpnn_run.py","--jsonl_path",f"{d}/parsed.jsonl","--fixed_positions_jsonl",f"{d}/fixed.jsonl",
             "--out_folder",d,"--num_seq_per_target","4","--sampling_temp","0.1","--seed","37","--batch_size","1",
             "--omit_AAs","C","--use_soluble_model","--model_name","v_48_020"],cwd=MPNN)
        fa=f"{d}/seqs/threaded.fa"
        if os.path.isfile(fa):
            recs=[]; 
            for line in open(fa):
                if line.startswith(">"): recs.append([line.strip(),""])
                elif recs: recs[-1][1]+=line.strip()
            designs=recs[1:]
            n=len(rns)
            # verify motif catalytic residues in designed seq: motifA R366 at ordinal, K184
            posA={rn:i for i,rn in enumerate(rns[:NA])}
            r366_idx=posA.get(366); k184_idx=posA.get(184)
            ok=[]
            for hdr,s in designs:
                good = len(s)==n and (r366_idx is None or s[r366_idx]=="R") and (k184_idx is None or s[k184_idx]=="K")
                ok.append(good)
            results[tag]={"N":n,"Lc":len(conn_out),"n_designs":len(designs),"cat_ok":sum(ok)}
            with open(f"{WORK}/all_rho_designs.fasta","a") as o:
                for i,(hdr,s) in enumerate(designs):
                    o.write(f">{tag}_d{i}\n{s}\n")
        else:
            results[tag]={"error":"no mpnn out"}
        print(tag,results[tag])
    json.dump(results,open(f"{WORK}/rho_mpnn_summary.json","w"),indent=2)
    print("done",len(results),"backbones")

if __name__=="__main__": main()
