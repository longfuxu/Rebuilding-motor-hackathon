#!/usr/bin/env python3
"""
M3/M4/M5 functional scoring of gp16 single-chain ring designs via the "7JQQ subunit-swap".

7JQQ = complete phi29 packaging motor: gp16 ATPase pentamer A-E (A,B,C carry AGS = ATP-gamma-S
+ Mg2+), coaxial dsDNA F,G, and pRNA pentamer K-O. We place each design into that frame and
score by GEOMETRY ONLY (never pLDDT/pTM).

KEY STRUCTURAL FACT (measured here): gp16 is a two-domain motor.
  - N-terminal ASCE ATPase domain  ~ native 4-205  (Walker-A 24-31, Walker-B D118/E119,
    arginine finger R146, and DNA contacts 55-130 / pRNA contacts 10-182)
  - C-terminal domain              ~ native 206-330 (DNA contacts 292/293/297/330,
    pRNA contacts 237-273)
A single rigid whole-monomer fit is therefore meaningless (inter-domain hinge). We anchor on
each domain separately.

Grounded functional contacts (7JQQ 4.5 A, native gp16 numbering):
  DNA (20): 55-60,82-83,98-100,125-130,292-293,297,330
  pRNA(29): 10-17,37-45,148-152,182,237-240,254,267-273

Constructs: native reference (A_apo apo pentamer + 7JQQ gp16 in-frame), cp233 (=C_design),
cp285, cp297. Design copies are circular permutations: per copy = seg1[native cp+1..330] +
GS-linker + seg2[native 4..cp]. Cut severs the native (cp)-(cp+1) peptide bond.

Run: /Users/longfu/miniforge3/bin/python3.13 score_m3_m4_m5.py
"""
import warnings; warnings.filterwarnings("ignore")
import json, re
import numpy as np
from scipy.spatial import cKDTree
import biotite.structure.io.pdbx as pdbx
import biotite.structure.io.pdb as pdb
import biotite.structure as struc
import biotite.sequence as bseq
import biotite.sequence.align as balign

BASE = "/Users/longfu/Developer/claude-science-hackthon/.claude/worktrees/export-gp16-results/gp16_design/"
OUT  = BASE + "outputs/m3_m4_m5/"

NATSEQ = ("MDKSLFYNPQKMLSYDRILNFVIGARGIGKSYAMKVYPINRFIKYGEQFIYVRRYKPELAKVSNYFNDVAQEFPDHELVVKGRRFYIDGKLAGWAIPLS"
          "VWQSEKSNAYPNVSTIVFDEFIREKDNSNYIPNEVSALLNLMDTVFRNRERVRCICLSNAVSVVNPYFLFFNLVPDVNKRFNVYDDALIEIPDSLDFSS"
          "ERRKTRFGRLIDGTEYGEMSLDNQFIGDSQVFIEKRSKDSKFVFSIVYNGFTLGVWVDVNQGLMYIDTAHDPSTKNVYTLTTDDLNENMMLITNYKNNY"
          "HLRKLASAFMNGYLRFDNQVIRNIAYELFRKMRIQ")

DNA_CONTACTS  = [55,56,57,59,60,82,83,98,99,100,125,126,127,128,129,130,292,293,297,330]
PRNA_CONTACTS = [10,11,13,14,15,16,17,37,40,41,44,45,148,149,150,152,182,237,238,239,240,254,267,268,269,270,271,272,273]
NDOM = lambda r: r <= 205
CDOM = lambda r: r >= 206
NDOM_ANCHOR = (20, 200)     # rigid ATPase core used to place each subunit
CDOM_ANCHOR = (234, 330)    # C-domain structural core (holds all C-domain contacts)
WALKER_A = list(range(24,32)); K_LYS = 30; WALKER_B = [118,119]; RFINGER = 146
POCKET = WALKER_A + WALKER_B + [RFINGER]

THREE2ONE = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
 'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
 'THR':'T','TRP':'W','TYR':'Y','VAL':'V','MSE':'M'}
GP16 = ['A','B','C','D','E']

def load(path):
    if path.endswith(".cif"):
        return pdbx.get_structure(pdbx.CIFFile.read(path), model=1)
    return pdb.get_structure(pdb.PDBFile.read(path), model=1)

# ------------------------------------------------------------------ design -> native map
MATRIX = balign.SubstitutionMatrix.std_protein_matrix()
def map_design_positions(seqstr):
    """position(0based) -> (native_resid or -1, copy 0-4, is_linker)."""
    n = len(seqstr); linker = np.zeros(n, bool)
    for m in re.finditer(r'(?:GGGGS){2,}', seqstr): linker[m.start():m.end()] = True
    segs=[]; i=0
    while i<n:
        if not linker[i]:
            j=i
            while j<n and not linker[j]: j+=1
            segs.append((i,j)); i=j
        else: i+=1
    natmap=np.full(n,-1); copy=np.full(n,-1); nat=bseq.ProteinSequence(NATSEQ)
    for k,(a0,a1) in enumerate(segs):
        al=balign.align_optimal(bseq.ProteinSequence(seqstr[a0:a1].replace('X','A')),nat,MATRIX,
                                local=True,gap_penalty=(-10,-1))[0]
        for r in range(al.trace.shape[0]):
            qi,ni=al.trace[r]
            if qi>=0 and ni>=0: natmap[a0+qi]=ni+1; copy[a0+qi]=k//2
    last=0
    for i in range(n):
        if copy[i]>=0: last=copy[i]
        else: copy[i]=last
    return natmap, copy, linker

def build_construct(name, kind, path):
    """Return dict with heavy AtomArray + per-atom arrays native_resid, copy, islink, and
    per-copy CA dict {copy: {native_resid: coord}}. Coords are the raw (own-frame) coords."""
    if kind == "inframe_native":
        a = load(BASE+"data/raw/7JQQ.cif")
        heavy = a[struc.filter_amino_acids(a) & np.isin(a.chain_id, GP16) & (a.element!="H")]
        native_resid = heavy.res_id.copy()
        copy = np.array([GP16.index(c) for c in heavy.chain_id])
        islink = np.zeros(heavy.array_length(), bool)
    elif kind == "apo":
        a = load(path)
        heavy = a[struc.filter_amino_acids(a) & (a.element!="H")]
        native_resid = heavy.res_id.copy()
        copy = np.array([GP16.index(c) for c in heavy.chain_id])
        islink = np.zeros(heavy.array_length(), bool)
    else:  # design single chain A
        a = load(path)
        heavy = a[struc.filter_amino_acids(a) & (a.chain_id=="A") & (a.element!="H")]
        ca = heavy[heavy.atom_name=="CA"]; order=np.argsort(ca.res_id); ca=ca[order]
        s="".join(THREE2ONE.get(r,'X') for r in ca.res_name)
        nm,cp,lk = map_design_positions(s)
        resid2 = {int(ca.res_id[i]):(int(nm[i]),int(cp[i]),bool(lk[i])) for i in range(len(s))}
        native_resid = np.array([resid2[int(r)][0] for r in heavy.res_id])
        copy         = np.array([resid2[int(r)][1] for r in heavy.res_id])
        islink       = np.array([resid2[int(r)][2] for r in heavy.res_id])
    # per-copy CA dicts (exclude linkers, require valid native id)
    cabyc=[dict() for _ in range(5)]
    ca_mask = heavy.atom_name=="CA"
    for x,nr,c,l in zip(heavy.coord[ca_mask], native_resid[ca_mask], copy[ca_mask], islink[ca_mask]):
        if nr>0 and not l: cabyc[c][int(nr)]=x
    return dict(name=name, kind=kind, heavy=heavy, coord=heavy.coord.copy(),
                native_resid=native_resid, copy=copy, islink=islink,
                atom_name=heavy.atom_name.copy(), cabyc=cabyc)

# ------------------------------------------------------------------ 7JQQ context
qc = build_construct("native_7JQQ","inframe_native",None)
q  = load(BASE+"data/raw/7JQQ.cif")
Q_DNA  = q[struc.filter_nucleotides(q) & np.isin(q.chain_id,['F','G']) & (q.element!="H")]
Q_PRNA = q[struc.filter_nucleotides(q) & np.isin(q.chain_id,['K','L','M','N','O']) & (q.element!="H")]
Q_AGS  = {ch: q[(q.res_name=="AGS")&(q.chain_id==ch)] for ch in ['A','B','C']}
Q_MG   = {ch: q[(q.res_name=="MG")&(q.chain_id==ch)]  for ch in ['A','B','C']}
DNA_tree  = cKDTree(Q_DNA.coord)
PRNA_tree = cKDTree(Q_PRNA.coord)
# 7JQQ per-chain CA dict (for anchoring)
q_ca = qc['cabyc']  # copy index == GP16 chain index

# ------------------------------------------------------------------ geometry helpers
def kabsch(fixed, mobile):
    fitted, tf = struc.superimpose(fixed, mobile)
    return tf, float(struc.rmsd(fixed, fitted))

def domain_ca(cabyc_copy, lo, hi):
    return {r:x for r,x in cabyc_copy.items() if lo<=r<=hi}

def anchor_transform(mob_ca, chain_idx, lo, hi):
    """Best-fit transform mapping mob copy's domain [lo,hi] onto 7JQQ chain[chain_idx]."""
    fd = domain_ca(q_ca[chain_idx], lo, hi); md = domain_ca(mob_ca, lo, hi)
    common = sorted(set(fd)&set(md))
    if len(common)<5: return None, None, 0
    F=np.array([fd[r] for r in common]); M=np.array([md[r] for r in common])
    tf, rmsd = kabsch(F, M)
    return tf, rmsd, len(common)

def best_assignment(cabyc, lo, hi):
    """Choose (direction, offset) mapping design copies -> 7JQQ chains by global N-domain ring fit."""
    best=None
    for direction in (+1,-1):
        for offset in range(5):
            F=[]; M=[]
            for c in range(5):
                ci=(offset+direction*c)%5
                fd=domain_ca(q_ca[ci],lo,hi); md=domain_ca(cabyc[c],lo,hi)
                for r in sorted(set(fd)&set(md)): F.append(fd[r]); M.append(md[r])
            F=np.array(F); M=np.array(M); tf,rmsd=kabsch(F,M)
            if best is None or rmsd<best[1]: best=((direction,offset),rmsd)
    return best[0], best[1]

def ring_axis(cabyc, lo=20, hi=200):
    # axis from the native-like ATPase (N) domain -> stable even when C-domain is displaced/broken
    cents=[np.array(list(domain_ca(cabyc[c],lo,hi).values())).mean(0) for c in range(5)]
    cents=np.array(cents); center=cents.mean(0)
    _,_,vt=np.linalg.svd(cents-center,full_matrices=False)
    normal=vt[2]
    return center, normal

def radial_of(coord, center, axis):
    d=coord-center; t=d@axis; return np.linalg.norm(d-np.outer(t,axis),axis=1)

def pore_profile(ca_coord, center, axis, step=2.0, halfwin=3.0, lo_pct=8, hi_pct=92):
    d=ca_coord-center; t=d@axis; perp=d-np.outer(t,axis); r=np.linalg.norm(perp,axis=1)
    tlo,thi=np.percentile(t,lo_pct),np.percentile(t,hi_pct)
    prof=[]
    for tt in np.arange(tlo,thi+step,step):
        sel=np.abs(t-tt)<halfwin
        if sel.sum()>=3: prof.append((float(tt),float(r[sel].min())))
    return prof

def residue_atoms(con, nr, cp, sidechain=True):
    m=(con['native_resid']==nr)&(con['copy']==cp)&(~con['islink'])
    if not m.any(): return None
    coord=con['coord'][m]; names=con['atom_name'][m]
    if sidechain:
        sc=coord[~np.isin(names,["N","C","O"])]
        return sc if len(sc)>0 else coord
    return coord

# ------------------------------------------------------------------ per-construct scoring
CONSTRUCTS = [
    ("native_7JQQ","inframe_native",None),
    ("A_apo","apo",BASE+"md/openmm_validation/inputs/A_apo.pdb"),
    ("cp233","design",BASE+"md/openmm_validation/inputs/C_design.cif"),
    ("cp285","design",BASE+"outputs/cp_grid_screen/cp285_int15_inter10.cif"),
    ("cp297","design",BASE+"outputs/cp_grid_screen/cp297_int15_inter10.cif"),
]
CUT = {"cp233":233,"cp285":285,"cp297":297}   # severed native peptide bond cp..cp+1

results={}
placements={}   # name -> transformed heavy coords (global N-domain ring fit) for viz/clash
for name,kind,path in CONSTRUCTS:
    con = qc if name=="native_7JQQ" else build_construct(name,kind,path)
    rec={"name":name,"kind":kind}

    # ---- M3a intrinsic pore (own frame, CA, N-domain-defined ring axis)
    center,axis = ring_axis(con['cabyc'])
    ca_mask = con['atom_name']=="CA"
    ca_all  = con['coord'][ca_mask]
    ca_nres = con['native_resid'][ca_mask]; ca_link = con['islink'][ca_mask]
    # all-CA pore = whole channel (a collapsed/displaced C-domain occludes it)
    prof = pore_profile(ca_all, center, axis)
    pr = np.array([r for _,r in prof])
    rec["M3a_min_pore_radius_CA"]=round(float(pr.min()),1)
    rec["M3a_pore_diam"]=round(float(pr.min())*2,1)
    rec["M3a_admits_dsDNA"]=bool(pr.min()>=9.0)   # ~apo channel; 20A dsDNA gripped
    # N-domain-only pore = ATPase channel alone (isolates the primary translocation channel)
    ndm = ca_mask.copy()
    nsel = (ca_nres>=20)&(ca_nres<=205)&(~ca_link)
    profN = pore_profile(ca_all[nsel], center, axis)
    rec["M3a_min_pore_radius_Ndom"]=round(float(np.array([r for _,r in profN]).min()),1)
    rec["_pore_profile"]=prof

    # ---- assignment + domain-integrity + placement
    if kind=="inframe_native":
        assign=("native",0); rec["assign_rmsd_Ndom"]=0.0
    else:
        assign, arms = best_assignment(con['cabyc'], *NDOM_ANCHOR)
        rec["assign_rmsd_Ndom"]=round(arms,2)
    direction,offset = (1,0) if kind=="inframe_native" else assign
    def chain_of(c): return c if kind=="inframe_native" else (offset+direction*c)%5

    # domain-integrity: internal-fold RMSD of each subunit's N/C domain vs its 7JQQ chain
    nd=[]; cd=[]; ndc_rmsd=[]; cdc_rmsd=[]
    for c in range(5):
        ci=chain_of(c)
        _,rn,_=anchor_transform(con['cabyc'][c],ci,*NDOM_ANCHOR); nd.append(rn)
        _,rcd,_=anchor_transform(con['cabyc'][c],ci,*CDOM_ANCHOR); cd.append(rcd)
    nd=[x for x in nd if x is not None]; cd=[x for x in cd if x is not None]
    rec["Ndom_intfold_rmsd"]=round(float(np.mean(nd)),2) if nd else None
    rec["Cdom_intfold_rmsd"]=round(float(np.mean(cd)),2) if cd else None

    # ---- global N-domain ring placement (single rigid transform) for clash + viz
    F=[];M=[]
    for c in range(5):
        ci=chain_of(c); fd=domain_ca(q_ca[ci],*NDOM_ANCHOR); md=domain_ca(con['cabyc'][c],*NDOM_ANCHOR)
        for r in sorted(set(fd)&set(md)): F.append(fd[r]); M.append(md[r])
    if kind=="inframe_native":
        placed=con['coord'].copy(); rec["global_place_rmsd"]=0.0
    else:
        tf,grms=kabsch(np.array(F),np.array(M)); placed=tf.apply(con['coord']); rec["global_place_rmsd"]=round(grms,2)
    placements[name]=dict(coord=placed, con=con)

    # ---- M3b / M4 : do the contact residues LINE the channel? (symmetric, own-frame radial test)
    # Native gp16 grips DNA on the channel wall, so DNA-contact sidechains sit at small radial
    # distance from the axis and point INWARD (sidechain tip closer to axis than its CA). This is
    # frame-independent (avoids 7JQQ's asymmetric partial-DNA engagement). Also record distance to
    # the actual 7JQQ DNA/pRNA (domain-anchored) in the detail table.
    def lining_scores(contact_list, tree):
        per=[]  # (copy, nat, domain, radial_tip, inward(bool), dist_to_NA)
        for c in range(5):
            ci=chain_of(c)
            tfN,_,_=anchor_transform(con['cabyc'][c],ci,*NDOM_ANCHOR)
            tfC,_,_=anchor_transform(con['cabyc'][c],ci,*CDOM_ANCHOR)
            for nr in contact_list:
                sc=residue_atoms(con,nr,c,sidechain=True)
                if sc is None: continue
                r_sc=radial_of(sc,center,axis)
                mca=(con['native_resid']==nr)&(con['copy']==c)&(~con['islink'])&(con['atom_name']=="CA")
                rca=float(radial_of(con['coord'][mca],center,axis)[0]) if mca.any() else float(r_sc.mean())
                rtip=float(r_sc.min())
                inward = rtip < rca
                tf = tfN if NDOM(nr) else tfC
                dna_d = float(tree.query(tf.apply(sc),k=1)[0].min()) if tf is not None else None
                per.append((c,nr,'N' if NDOM(nr) else 'C',round(rtip,2),bool(inward),
                            round(dna_d,2) if dna_d is not None else None))
        return per
    dna_per  = lining_scores(DNA_CONTACTS, DNA_tree)
    prna_per = lining_scores(PRNA_CONTACTS, PRNA_tree)
    def summ(per):
        out={}
        for dom in ('N','C'):
            L=[p for p in per if p[2]==dom]
            if not L: out[dom]=(None,None,0); continue
            inward=sum(1 for p in L if p[4])
            out[dom]=(round(float(np.median([p[3] for p in L])),1), inward, len(L))
        return out
    ds=summ(dna_per); ps=summ(prna_per)
    rec["M3b_DNA_med_radial_Ndom"]=ds['N'][0]; rec["M3b_DNA_inward_Ndom"]=f"{ds['N'][1]}/{ds['N'][2]}"
    rec["M3b_DNA_med_radial_Cdom"]=ds['C'][0]; rec["M3b_DNA_inward_Cdom"]=f"{ds['C'][1]}/{ds['C'][2]}"
    rec["M4_pRNA_med_radial_Ndom"]=ps['N'][0]; rec["M4_pRNA_inward_Ndom"]=f"{ps['N'][1]}/{ps['N'][2]}"
    rec["M4_pRNA_med_radial_Cdom"]=ps['C'][0]; rec["M4_pRNA_inward_Cdom"]=f"{ps['C'][1]}/{ps['C'][2]}"
    rec["_dna_per"]=dna_per; rec["_prna_per"]=prna_per

    # ---- M3c clash vs DNA (global N-placement); linker atoms specifically
    dprot=DNA_tree.query(placed,k=1)[0]
    rec["M3c_clash_atoms_vs_DNA(<2.5A)"]=int((dprot<2.5).sum())
    rec["M3c_clash_linker_atoms"]=int(((dprot<2.5)&con['islink']).sum())
    rec["M3c_min_dist_prot_DNA"]=round(float(dprot.min()),2)

    # ---- M3d cut disruption: does the CP sever/relocate a DNA contact?
    if name in CUT:
        cutres=CUT[name]
        # native peptide bond cutres..cutres+1 is broken; cutres becomes a new C-terminus of seg2,
        # cutres+1 the new N-terminus (seg1) of the copy; both are linker-adjacent termini.
        rec["M3d_cut_bond"]=f"{cutres}-{cutres+1}"
        rec["M3d_cut_hits_DNA_contact"]=bool(cutres in DNA_CONTACTS or (cutres+1) in DNA_CONTACTS)
        rec["M3d_cut_hits_pRNA_contact"]=bool(cutres in PRNA_CONTACTS or (cutres+1) in PRNA_CONTACTS)
    else:
        rec["M3d_cut_bond"]="none(native)"; rec["M3d_cut_hits_DNA_contact"]=False; rec["M3d_cut_hits_pRNA_contact"]=False

    # ---- M5 ATP pocket: anchor N-domain onto AGS-bearing 7JQQ chain A; measure catalytic geometry
    #      vs that chain's AGS/Mg. Reference = 7JQQ chain A's own catalytic residues.
    agsA=Q_AGS['A']; mgA=Q_MG['A']
    p_beta = agsA[np.isin(agsA.atom_name,["PB","O1B","O2B","O3B"])].coord
    p_gamma= agsA[np.isin(agsA.atom_name,["PG","O1G","O2G","O3G","S1G"])].coord
    mg_xyz = mgA.coord[0] if mgA.array_length() else None
    k30_d=[]; d118_d=[]; e119_d=[]; r146_d=[]; pkt_rmsd=[]
    for c in range(5):
        # anchor this subunit's N-domain onto chain A (AGS present)
        tf,rms,_=anchor_transform(con['cabyc'][c],0,*NDOM_ANCHOR)
        if tf is None: continue
        def far(nr,atoms):
            m=(con['native_resid']==nr)&(con['copy']==c)&(~con['islink'])&np.isin(con['atom_name'],atoms)
            if not m.any(): return None
            return tf.apply(con['coord'][m])
        nz=far(K_LYS,["NZ"])
        if nz is not None and len(p_gamma):
            k30_d.append(float(cKDTree(np.vstack([p_beta,p_gamma])).query(nz,k=1)[0].min()))
        dcg=far(118,["OD1","OD2","CG"]);
        if dcg is not None and mg_xyz is not None: d118_d.append(float(np.linalg.norm(dcg-mg_xyz,axis=1).min()))
        ecg=far(119,["OE1","OE2","CD"]);
        if ecg is not None and mg_xyz is not None: e119_d.append(float(np.linalg.norm(ecg-mg_xyz,axis=1).min()))
        rcz=far(146,["CZ","NH1","NH2","NE"])
        if rcz is not None: r146_d.append(float(cKDTree(agsA.coord).query(rcz,k=1)[0].min()))
        # pocket CA rmsd
        pf=[];pm=[]
        for nr in POCKET:
            if nr in q_ca[0] and nr in con['cabyc'][c]:
                pf.append(q_ca[0][nr]); pm.append(con['cabyc'][c][nr])
        if len(pf)>=4:
            pff=np.array(pf); pmm=tf.apply(np.array(pm)); pkt_rmsd.append(float(np.sqrt(((pff-pmm)**2).sum(1).mean())))
    mean=lambda L: round(float(np.mean(L)),2) if L else None
    rec["M5_K30_to_ATP_phosphate"]=mean(k30_d)
    rec["M5_D118_to_Mg"]=mean(d118_d)
    rec["M5_E119_to_Mg"]=mean(e119_d)
    rec["M5_R146_to_ATP"]=mean(r146_d)
    rec["M5_pocket_CA_rmsd_vs_native"]=mean(pkt_rmsd)

    results[name]=rec
    print(f"[{name:12s}] pore all/Ndom {rec['M3a_min_pore_radius_CA']:>4}/{rec['M3a_min_pore_radius_Ndom']:>4}  "
          f"intfold N/C {rec['Ndom_intfold_rmsd']}/{rec['Cdom_intfold_rmsd']}  "
          f"DNA inward N/C {rec['M3b_DNA_inward_Ndom']}|{rec['M3b_DNA_inward_Cdom']}  "
          f"linkerClash {rec['M3c_clash_linker_atoms']}  D118-Mg {rec['M5_D118_to_Mg']}  pktRMSD {rec['M5_pocket_CA_rmsd_vs_native']}")

# ------------------------------------------------------------------ write outputs
clean={k:{kk:vv for kk,vv in v.items() if not kk.startswith('_')} for k,v in results.items()}
with open(OUT+"m3_m4_m5_scores.json","w") as f: json.dump(clean,f,indent=2)

cols=["name","kind","M3a_min_pore_radius_CA","M3a_min_pore_radius_Ndom","M3a_pore_diam","M3a_admits_dsDNA",
      "Ndom_intfold_rmsd","Cdom_intfold_rmsd",
      "M3b_DNA_med_radial_Ndom","M3b_DNA_inward_Ndom","M3b_DNA_med_radial_Cdom","M3b_DNA_inward_Cdom",
      "M4_pRNA_med_radial_Ndom","M4_pRNA_inward_Ndom","M4_pRNA_med_radial_Cdom","M4_pRNA_inward_Cdom",
      "M3c_clash_atoms_vs_DNA(<2.5A)","M3c_clash_linker_atoms","M3c_min_dist_prot_DNA","M3d_cut_bond",
      "M3d_cut_hits_DNA_contact","M3d_cut_hits_pRNA_contact",
      "M5_D118_to_Mg","M5_E119_to_Mg","M5_R146_to_ATP","M5_pocket_CA_rmsd_vs_native",
      "assign_rmsd_Ndom","global_place_rmsd"]
import csv
with open(OUT+"m3_m4_m5_scores.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(cols)
    for name in results:
        w.writerow([results[name].get(c,"") for c in cols])

# per-contact detail table (for M3d / audit)
with open(OUT+"contact_detail.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["construct","nucleic","copy","native_resid","domain","radial_tip_A","inward","dist_to_NA_A"])
    for name in results:
        for tag,key in [("DNA","_dna_per"),("pRNA","_prna_per")]:
            for (c,nr,dom,rtip,inward,d) in results[name][key]:
                w.writerow([name,tag,c,nr,dom,rtip,inward,d])

import pickle
dump={n:{'coord':placements[n]['coord'],
         'native_resid':placements[n]['con']['native_resid'],
         'copy':placements[n]['con']['copy'],
         'islink':placements[n]['con']['islink'],
         'atom_name':placements[n]['con']['atom_name'],
         'res_id':placements[n]['con']['heavy'].res_id.copy()} for n in placements}
dump['_DNA']=Q_DNA.coord; dump['_DNA_atom']=Q_DNA.atom_name.copy(); dump['_DNA_chain']=Q_DNA.chain_id.copy(); dump['_DNA_resid']=Q_DNA.res_id.copy()
dump['_PRNA']=Q_PRNA.coord; dump['_PRNA_atom']=Q_PRNA.atom_name.copy(); dump['_PRNA_chain']=Q_PRNA.chain_id.copy(); dump['_PRNA_resid']=Q_PRNA.res_id.copy()
dump['_DNA_CONTACTS']=DNA_CONTACTS; dump['_PRNA_CONTACTS']=PRNA_CONTACTS
with open(OUT+"_placements.pkl","wb") as f: pickle.dump(dump,f)
print("\nWROTE:", OUT+"m3_m4_m5_scores.{json,csv}, contact_detail.csv, _placements.pkl")
