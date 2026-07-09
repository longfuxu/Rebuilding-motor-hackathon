"""ClpX hexamer M1/M2 scorer — retrospective validation of fold-in-context framework.
M2: trans arginine-finger R307 (in-copy 246) guanidinium -> neighbor Walker-A (in-copy 64-73)
    min heavy-atom distance, engaged <8A, sequential cyclic over 6 copies, floor >=5/6.
    Calibrated on functional substrate-bound hexamer 6PP5 (5/6 engaged ~5.5A, seam 9.1A).
M1: ring geometry from 6 copy centroids: radius_CV<0.35 compact, planarity RMS,
    sequential_consistent = designed cyclic order == spatial angular order.
Global pTM NOT used (per campaign rule). Native ClpX numbering per UniProt P0A6H1 (E. coli ClpX; sequence fetched + catalytic residues verified in-session); ClpXdN construct res 62-424.
"""
import numpy as np, gemmi
def score_clpx(cif_path, spans, monomer_len, argfinger_incopy, walker_incopy_range,
               engaged_thr=8.0, radius_cv_thr=0.35, guan=("NE","CZ","NH1","NH2")):
    """
    ClpX hexamer M1/M2 scorer (single chain A, `spans` = list of 6 (lo,hi) 1-based copy spans).
    M2: trans R307(copy i) guanidinium -> Walker-A(copy i+1, cyclic) min heavy-atom dist; engaged < thr.
    M1: ring geometry from copy centroids (radius_CV, planarity) + sequential_consistent
        (designed cyclic order == spatial angular order).
    """
    st=gemmi.read_structure(cif_path); st.setup_entities()
    m=st[0]
    ch=m[0]  # single chain
    # residue -> heavy atoms
    at={}
    for r in ch:
        if r.name=="HOH": continue
        at[r.seqid.num]={a.name:np.array([a.pos.x,a.pos.y,a.pos.z]) for a in r if a.element.name!="H"}
    ncopy=len(spans)
    def construct_pos(ci, incopy): return spans[ci][0]+incopy-1
    # ---- M2: sequential trans engagement ----
    wlo,whi=walker_incopy_range
    dists=[]; engaged=[]
    for i in range(ncopy):
        j=(i+1)%ncopy
        rpos=construct_pos(i, argfinger_incopy)
        gu=[at[rpos][k] for k in guan if rpos in at and k in at[rpos]]
        wa=[]
        for ic in range(wlo,whi+1):
            wp=construct_pos(j,ic)
            if wp in at: wa.extend(at[wp].values())
        if not gu or not wa:
            dists.append(None); engaged.append(False); continue
        G=np.array(gu); W=np.array(wa)
        d=float(np.min(np.sqrt(((G[:,None,:]-W[None,:,:])**2).sum(-1))))
        dists.append(round(d,2)); engaged.append(d<engaged_thr)
    m2=sum(engaged)
    # ---- M1: ring geometry ----
    cents=[]
    for (lo,hi) in spans:
        pts=[at[p]["CA"] for p in range(lo,hi+1) if p in at and "CA" in at[p]]
        cents.append(np.mean(pts,axis=0))
    cents=np.array(cents); ring_c=cents.mean(0)
    radii=np.linalg.norm(cents-ring_c,axis=1)
    radius_cv=float(np.std(radii)/np.mean(radii))
    compact=radius_cv<radius_cv_thr
    # planarity: fit plane, rms of out-of-plane
    u,s,vt=np.linalg.svd(cents-ring_c); normal=vt[2]
    planarity_rms=float(np.sqrt(np.mean((( cents-ring_c)@normal)**2)))
    # sequential consistency: spatial angular order around normal vs designed order 0..5
    x=vt[0]; y=vt[1]
    ang=np.array([np.arctan2((c-ring_c)@y,(c-ring_c)@x) for c in cents])
    order=list(np.argsort(ang))
    # designed order is 0,1,2,3,4,5; spatial order should be a cyclic rotation (either direction)
    def is_cyclic(o):
        n=len(o)
        for d in (1,-1):
            for start in range(n):
                if o==[ (start+d*k)%n for k in range(n)]: return True
        return False
    seq_consistent=is_cyclic(order)
    return {"m2_engaged":m2,"m2_total":ncopy,"m2_dists":dists,"m2_engaged_bool":engaged,
            "radius_cv":round(radius_cv,3),"compact_ring":compact,
            "planarity_rms":round(planarity_rms,2),"spatial_order":order,
            "sequential_consistent":bool(seq_consistent)}
