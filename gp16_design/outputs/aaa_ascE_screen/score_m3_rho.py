#!/usr/bin/env python3
"""M3 channel gate for single-chain Rho ring folds (RNA translocase).

The atlas cautionary lesson: an apo single-chain ring can CLOSE the coupler interface
(M2) yet still be non-functional if the fusion/linker fouls the translocation channel.
M3 therefore checks, on the folded ring, GEOMETRY ONLY (never pLDDT/pTM):

  M3a  pore radius profile along the ring axis: the minimum CA-based pore radius must
       admit ssRNA (~5 A radius clearance). Reported as min_pore_radius and admits_ssRNA.
  M3b  functional-residue integrity: the 6 copies' Walker-A K184 and arginine-finger R366
       must sit at native-like radial positions (they line the pore / interface). We compare
       the median radial position of K184 and R366 CA across copies to the native 3ICE ring;
       a large deviation means the linker distorted the active site.
  M3c  linker intrusion: no designed-linker CA may sit within the substrate lumen
       (radial position < min_pore_radius + 2 A) — a linker in the lumen blocks RNA.

Copies are copy ranges within chain A (from the fold manifest). Rho monomer = res175-414
mapped to within-copy offsets; native landmarks K184 (Walker-A), R366 (arginine finger),
first modeled native res = 175.
Run: score_m3_rho.py <fold.cif> --copies "A:1-280,..." [--linker-spans ...]
"""
import sys, argparse, json
import numpy as np

def parse_atoms(path):
    out=[]
    if path.lower().endswith((".cif",".mmcif")):
        cols={}; order=[]; inloop=False
        for L in open(path):
            s=L.strip()
            if s.startswith("_atom_site."):
                cols[s]=len(order); order.append(s); inloop=True; continue
            if inloop and (s.startswith("_") or s.startswith("loop_") or s=="#"):
                if not s.startswith("_atom_site."): inloop=False
                continue
            if inloop and s and not s.startswith("#"):
                f=s.split()
                if len(f)<len(order): continue
                try:
                    ch=f[cols.get("_atom_site.auth_asym_id",cols.get("_atom_site.label_asym_id"))]
                    rn=int(f[cols["_atom_site.auth_seq_id"]] if "_atom_site.auth_seq_id" in cols else f[cols["_atom_site.label_seq_id"]])
                    an=f[cols["_atom_site.label_atom_id"]].strip('"')
                    x=float(f[cols["_atom_site.Cartn_x"]]); y=float(f[cols["_atom_site.Cartn_y"]]); z=float(f[cols["_atom_site.Cartn_z"]])
                    out.append((ch,rn,an,np.array([x,y,z])))
                except: pass
    else:
        for L in open(path):
            if L[:4]=="ATOM":
                ch=L[21]; rn=int(L[22:26]); an=L[12:16].strip()
                out.append((ch,rn,an,np.array([float(L[30:38]),float(L[38:46]),float(L[46:54])])))
    return out

def ring_axis(cas):
    C=np.array(cas); ctr=C.mean(0); X=C-ctr
    _,_,Vt=np.linalg.svd(X); return ctr, Vt[2]

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("cif")
    ap.add_argument("--copies",required=True)      # A:1-280,A:281-560,...
    ap.add_argument("--first-native",type=int,default=175)
    ap.add_argument("--k184",type=int,default=184)
    ap.add_argument("--r366",type=int,default=366)
    ap.add_argument("--linker-len",type=int,default=None, help="if set, treat last <linker-len> res of each copy span as linker")
    args=ap.parse_args()
    atoms=parse_atoms(args.cif)
    # copy spans
    spans=[]
    for tok in args.copies.split(","):
        ch,rng=tok.split(":"); a,b=rng.split("-"); spans.append((ch,int(a),int(b)))
    # all CA
    ca={} 
    for ch,rn,an,xyz in atoms:
        if an=="CA": ca[(ch,rn)]=xyz
    # per-copy CA lists + landmark coords
    copy_cas=[]; k184s=[]; r366s=[]; linker_cas=[]
    for (ch,a,b) in spans:
        cc=[ca[(ch,r)] for r in range(a,b+1) if (ch,r) in ca]
        copy_cas.append(cc)
        # within-copy offset of native residue = (native - first_native); res number in construct = a + offset
        off_k=args.k184-args.first_native; off_r=args.r366-args.first_native
        if (ch,a+off_k) in ca: k184s.append(ca[(ch,a+off_k)])
        if (ch,a+off_r) in ca: r366s.append(ca[(ch,a+off_r)])
        if args.linker_len:
            for r in range(b-args.linker_len+1,b+1):
                if (ch,r) in ca: linker_cas.append(ca[(ch,r)])
    allca=[c for cc in copy_cas for c in cc]
    ctr,axis=ring_axis([np.mean(cc,0) for cc in copy_cas])
    def radial(p):
        v=np.array(p)-ctr; return float(np.linalg.norm(v-np.dot(v,axis)*axis))
    # M3a pore radius profile: bin residues by axial coord, min radial per bin = pore wall
    axial=np.array([np.dot(np.array(p)-ctr,axis) for p in allca])
    rad=np.array([radial(p) for p in allca])
    nb=20; edges=np.linspace(axial.min(),axial.max(),nb+1)
    pore=[]
    for i in range(nb):
        sel=(axial>=edges[i])&(axial<edges[i+1])
        if sel.sum()>=3: pore.append(rad[sel].min())
    min_pore=float(min(pore)) if pore else None
    admits=bool(min_pore is not None and min_pore>=5.0)
    # M3b functional residue radial positions
    k184_rad=[radial(p) for p in k184s]; r366_rad=[radial(p) for p in r366s]
    # M3c linker intrusion
    linker_min_rad=float(min(radial(p) for p in linker_cas)) if linker_cas else None
    linker_in_lumen = bool(linker_min_rad is not None and min_pore is not None and linker_min_rad < min_pore+2.0)
    res={"cif":args.cif.split("/")[-1],
         "M3a_min_pore_radius":round(min_pore,1) if min_pore else None,
         "M3a_admits_ssRNA":admits,
         "M3b_K184_med_radial":round(float(np.median(k184_rad)),1) if k184_rad else None,
         "M3b_R366_med_radial":round(float(np.median(r366_rad)),1) if r366_rad else None,
         "M3b_n_K184":len(k184s),"M3b_n_R366":len(r366s),
         "M3c_linker_min_radial":round(linker_min_rad,1) if linker_min_rad else None,
         "M3c_linker_in_lumen":linker_in_lumen}
    print(json.dumps(res))

if __name__=="__main__": main()
