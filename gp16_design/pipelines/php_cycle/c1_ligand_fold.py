# C1 — ligand-conditioned endpoints via Boltz-2 NIM. Does ATP open the ring (planar->helical)?
import os, json, time, urllib.request, urllib.error, math
import numpy as np

NIM="https://health.api.nvidia.com/v1/biology/mit/boltz2/predict"
REPO="/Users/longfu/Developer/claude-science-hackthon/.claude/worktrees/export-gp16-results/gp16_design"
OUT=f"{REPO}/outputs/php_cycle/C1_boltz"; os.makedirs(OUT, exist_ok=True)
KEY=[l.split("=",1)[1].strip().strip('"') for l in open("/Users/longfu/Developer/claude-science-hackthon/.env")
     if l.startswith("NVIDIA_API_KEY")][0]

NATIVE=("SLFYNPQKMLSYDRILNFVIGARGIGKSYAMKVYPINRFIKYGEQFIYVRRYKPELAKVSNYFNDVAQEFPDHELVVKGRRFYIDGKLAGWAIPLSVWQSEKS"
        "NAYPNVSTIVFDEFIREKDNSNYIPNEVSALLNLMDTVFRNRERVRCICLSNAVSVVNPYFLFFNLVPDVNKRFNVYDDALIEIPDSLDFSSERRKTRFGRLID"
        "GTEYGEMSLDNQFIGDSQVFIEKRSKDSKFVFSIVYNGFTLGVWVDVNQGLMYIDTAHDPSTKNVYTLTTDDLNENMMLITNYKNNYHLRKLASAFMNGYLRFD"
        "NQVIRNIAYELFRKMR")
CORE=open(f"{REPO}/pipelines/tiled_msa_fold/gp16_core.a3m").read()
CP233=json.load(open(f"{REPO}/pipelines/tiled_msa_fold/manifests/cp233_WT.json"))["sequence"]
TILED=open(f"{REPO}/outputs/directed_evolution/folds/cp233_WTrep.tiled.a3m").read()
# single-chain 5 copies: residue blocks (1-based, inclusive)
SC_BLOCKS=[(1,342),(353,694),(705,1046),(1057,1398),(1409,1750)]

def poly(cid, seq, a3m): return {"id":cid,"molecule_type":"protein","sequence":seq,
                                 "msa":{"default":{"a3m":{"alignment":a3m,"format":"a3m"}}}}
def lig(cid, ccd): return {"id":cid,"ccd":ccd}
def atp_mg(n=5): return [lig(f"L{i+1}","ATP") for i in range(n)]+[lig(f"G{i+1}","MG") for i in range(n)]
def ags_mg(n=5): return [lig(f"L{i+1}","AGS") for i in range(n)]+[lig(f"G{i+1}","MG") for i in range(n)]

CONDS=[
 ("native_apo","native",[poly(c,NATIVE,CORE) for c in "ABCDE"],[]),
 ("native_atp","native",[poly(c,NATIVE,CORE) for c in "ABCDE"],atp_mg()),
 ("native_ags","native",[poly(c,NATIVE,CORE) for c in "ABCDE"],ags_mg()),
 ("sc_apo","sc",[poly("A",CP233,TILED)],[]),
 ("sc_atp","sc",[poly("A",CP233,TILED)],atp_mg()),
]

def call(polymers, ligands, tries=6):
    body=json.dumps({"polymers":polymers,"ligands":ligands,
                     "recycling_steps":3,"sampling_steps":50,"diffusion_samples":1}).encode()
    for a in range(tries):
        try:
            r=urllib.request.Request(NIM,data=body,headers={"Authorization":f"Bearer {KEY}",
                                     "Content-Type":"application/json"})
            t0=time.time(); resp=json.load(urllib.request.urlopen(r,timeout=1800))
            return resp, round(time.time()-t0,1)
        except urllib.error.HTTPError as e:
            b=e.read()[:400].decode(errors="replace")
            if e.code in (429,500,502,503,504): print(f"  HTTP {e.code} backoff",flush=True); time.sleep(20*(a+1)); continue
            return {"__err":f"HTTP {e.code}: {b}"},0
        except Exception as e:
            print(f"  {repr(e)[:120]} retry",flush=True); time.sleep(15)
    return {"__err":"retries exhausted"},0

def ca_by_chain(cif):
    """parse CA coords; returns dict chain-> list[(seq_id,xyz)] for protein CA."""
    hdr=[]; data=[]; inloop=False
    for l in open(cif):
        if l.startswith("_atom_site."): hdr.append(l.strip()); inloop=True
        elif inloop and (l.startswith("ATOM") or l.startswith("HETATM")): data.append(l.split())
        elif inloop and l.startswith("#"): break
    def ci(n):
        for i,h in enumerate(hdr):
            if h.endswith("."+n): return i
        return None
    a=ci("label_atom_id"); ch=ci("label_asym_id"); sq=ci("label_seq_id")
    x,y,z=ci("Cartn_x"),ci("Cartn_y"),ci("Cartn_z"); comp=ci("label_comp_id")
    out={}
    for r in data:
        if len(r)<=max(a,ch,sq,x,y,z): continue
        if r[a]!="CA": continue
        try: xyz=np.array([float(r[x]),float(r[y]),float(r[z])])
        except: continue
        try: s=int(r[sq])
        except: s=len(out.get(r[ch],[]))
        out.setdefault(r[ch],[]).append((s,xyz))
    return out

def subunit_centroids(cif, kind):
    cc=ca_by_chain(cif)
    cents=[]
    if kind=="native":
        for c in sorted(cc):
            pts=np.array([p for _,p in cc[c]])
            if len(pts)>50: cents.append(pts.mean(0))   # protein chains only
    else:  # single-chain: split the one chain by residue blocks
        # pick the chain with the most CA (the protein)
        c=max(cc, key=lambda k:len(cc[k]))
        resmap={s:p for s,p in cc[c]}
        for (a,b) in SC_BLOCKS:
            pts=np.array([resmap[s] for s in range(a,b+1) if s in resmap])
            if len(pts)>50: cents.append(pts.mean(0))
    return np.array(cents)

def geom(cents):
    """planarity + axial spread (P vs H discriminator) + per-subunit axial z + radius."""
    if len(cents)<3: return None
    c0=cents.mean(0); X=cents-c0
    # best-fit plane normal = smallest-singular-vector
    U,S,Vt=np.linalg.svd(X); axis=Vt[2]
    z=X@axis                      # axial coord per subunit (Å)
    inplane=X-np.outer(z,axis)
    radius=np.linalg.norm(inplane,axis=1).mean()
    return dict(n=len(cents), axial_spread=float(z.max()-z.min()),
                planarity_rms=float(np.sqrt((z**2).mean())),
                radius=float(radius), z=[round(float(v),1) for v in np.sort(z)])

def main():
    rows=[]
    for name,kind,polys,ligs in CONDS:
        print(f"\n=== {name} ({kind}, {len(polys)} chains, {len(ligs)} ligands) ===",flush=True)
        resp,dt=call(polys,ligs)
        if "__err" in resp:
            print(f"  FAILED {resp['__err'][:160]}",flush=True); rows.append(dict(name=name,ok=False,err=resp['__err'][:200])); continue
        cif=f"{OUT}/{name}.cif"; open(cif,"w").write(resp["structures"][0]["structure"])
        conf={k:resp.get(k) for k in ("confidence_scores","ptm_scores","iptm_scores") if k in resp}
        g=geom(subunit_centroids(cif,kind))
        print(f"  ok {dt}s  geom={g}  conf={conf}",flush=True)
        rows.append(dict(name=name,ok=True,dt=dt,**(g or {}),conf=conf))
        json.dump(rows,open(f"{OUT}/C1_results.json","w"),indent=1)
    print("\n===== SUMMARY (axial_spread: planar~small, helical~large) =====",flush=True)
    for r in rows:
        if r.get("ok"): print(f"  {r['name']:12s} spread {r.get('axial_spread','?'):>6} Å  planarity_rms {r.get('planarity_rms','?')}  radius {r.get('radius','?')}  z={r.get('z')}",flush=True)
        else: print(f"  {r['name']:12s} FAILED {r.get('err','')[:80]}",flush=True)
    print("\nDONE",flush=True)

if __name__=="__main__":
    main()
