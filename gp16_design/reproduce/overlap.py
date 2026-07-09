import sys, numpy as np, scipy.linalg as sla
sys.path.insert(0,"gp16_design/reproduce")
from score_m2 import parse_atoms
def ca_by_chain(path):
    d={}
    for c,rn,an,xyz,b in parse_atoms(path):
        if an=="CA": d.setdefault(c,{})[rn]=np.array(xyz)
    return d
apo=ca_by_chain("gp16_design/outputs/structures/cycle3/native_ring__boltz__s1__native_ring_model_0.pdb")
jqq=ca_by_chain("gp16_design/data/raw/7JQQ.cif")
chs=list("ABCDE")
res=sorted(set.intersection(*[set(apo[c]) for c in chs], *[set(jqq[c]) for c in chs]))  # common residues
print(f"common residues per subunit: {len(res)} (res {res[0]}-{res[-1]}); ring = {5*len(res)} CA")
Xa=np.array([apo[c][r] for c in chs for r in res])          # apo ring, fixed order
def kabsch(P,Q):    # superpose P onto Q, return aligned P and RMSD
    Pc=P-P.mean(0); Qc=Q-Q.mean(0)
    V,S,Wt=np.linalg.svd(Pc.T@Qc); d=np.sign(np.linalg.det(V@Wt)); D=np.diag([1,1,d])
    R=V@D@Wt; Pal=Pc@R+Q.mean(0); return Pal, np.sqrt(((Pal-Q)**2).sum(1).mean())
# try cyclic rotations + reflection of 7JQQ chain assignment
best=None
for order in [chs[i:]+chs[:i] for i in range(5)]+[ (chs[::-1])[i:]+(chs[::-1])[:i] for i in range(5)]:
    Xj=np.array([jqq[order[k]][r] for k in range(5) for r in res])
    Xal,rmsd=kabsch(Xj,Xa)
    if best is None or rmsd<best[0]: best=(rmsd,order,Xal)
rmsd,order,Xj_al=best
print(f"best chain assignment apo(ABCDE)->7JQQ{order}, ring RMSD apo<->7JQQ = {rmsd:.1f} A (they differ = the transition)")
Delta=(Xa - Xj_al).ravel(); Delta/=np.linalg.norm(Delta)   # planar->helical difference vector (unit)
# ANM on apo (common residues)
N=len(Xa); cut=13.0; H=np.zeros((3*N,3*N))
for i in range(N):
    dd=Xa-Xa[i]; r2=(dd*dd).sum(1); r2[i]=1e9
    for j in np.where(r2<cut*cut)[0]:
        e=(Xa[j]-Xa[i]); e/=np.linalg.norm(e); k=np.outer(e,e)
        H[3*i:3*i+3,3*i:3*i+3]+=k; H[3*i:3*i+3,3*j:3*j+3]-=k
w,V=sla.eigh(H, subset_by_index=[0,25])
ov=np.array([abs(V[:,m]@Delta) for m in range(len(w))])    # overlap of each mode with the transition
print(f"\nsingle-mode overlaps with apo->7JQQ transition (modes 7-15):")
for m in range(6,16):
    print(f"  mode {m+1:>2}: overlap {ov[m]:.2f}" + ("   <== soft mode explains the transition" if ov[m]>0.3 else ""))
for K in (2,5,10,20):
    cum=np.sqrt((ov[6:6+K]**2).sum())
    print(f"cumulative overlap, {K} softest modes: {cum:.2f}  ({cum**2*100:.0f}% of the transition captured)")
