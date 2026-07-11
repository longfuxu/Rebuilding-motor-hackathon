import sys, json, os
sys.path.insert(0,"/Users/longfu/.claude/jobs/48b7a543/tmp")
from c1_ligand_fold import (call, poly, lig, atp_mg, ags_mg, geom, subunit_centroids,
                            NATIVE, CORE, CP233, TILED, OUT)
def dna(cid, seq): return {"id":cid,"molecule_type":"dna","sequence":seq}
D1="CGCGAATTCGCGATCGATCG"; D2="CGATCGATCGCGAATTCGCG"
DNA=[dna("D1",D1), dna("D2",D2)]
CONDS=[
 ("native_apo_dna","native",[poly(c,NATIVE,CORE) for c in "ABCDE"]+DNA,[]),
 ("native_atp_dna","native",[poly(c,NATIVE,CORE) for c in "ABCDE"]+DNA,atp_mg()),
 ("native_ags_dna","native",[poly(c,NATIVE,CORE) for c in "ABCDE"]+DNA,ags_mg()),
 ("sc_atp_dna","sc",[poly("A",CP233,TILED)]+DNA,atp_mg()),
]
rows=[]
for name,kind,polys,ligs in CONDS:
    print(f"\n=== {name} ({sum(1 for p in polys if p['molecule_type']=='protein')} prot + 2 DNA, {len(ligs)} lig) ===",flush=True)
    resp,dt=call(polys,ligs)
    if "__err" in resp: print(f"  FAILED {resp['__err'][:180]}",flush=True); rows.append((name,None,resp['__err'][:120])); continue
    cif=f"{OUT}/{name}.cif"; open(cif,"w").write(resp["structures"][0]["structure"])
    g=geom(subunit_centroids(cif,kind))
    print(f"  ok {dt}s  spread={g['axial_spread']:.2f}A planarity={g['planarity_rms']:.2f} radius={g['radius']:.1f} z={g['z']}",flush=True)
    rows.append((name,g,None))
print("\n===== +DNA SUMMARY (planar~0.1A, 7JQQ helical~4.8A) =====",flush=True)
for name,g,err in rows:
    print(f"  {name:16s} spread {g['axial_spread']:.2f} A" if g else f"  {name:16s} FAIL {err}",flush=True)
print("DONE",flush=True)
