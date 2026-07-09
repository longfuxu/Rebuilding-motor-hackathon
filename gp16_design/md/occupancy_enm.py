"""(b) per-subunit apo->7JQQ opening vs the 3-ATP occupancy pattern, and
(a) occupancy-dependent elastic-network linear-response: does adding ATP one-by-one open the ring?
Free, local, numpy/scipy. apo = closed planar (Boltz native ring); 7JQQ = 3-ATP partial-occupancy helical."""
import sys, numpy as np, scipy.linalg as sla
sys.path.insert(0,"gp16_design/reproduce")
from score_m2 import parse_atoms
def ca(path):
    d={}
    for c,rn,an,xyz,b in parse_atoms(path):
        if an=="CA": d.setdefault(c,{})[rn]=np.array(xyz)
    return d
apo=ca("gp16_design/outputs/structures/cycle3/native_ring__boltz__s1__native_ring_model_0.pdb")
jqq=ca("gp16_design/data/raw/7JQQ.cif"); chs=list("ABCDE")
res=sorted(set.intersection(*[set(apo[c]) for c in chs],*[set(jqq[c]) for c in chs]))
Xa=np.array([apo[c][r] for c in chs for r in res]); n=len(res)
def kabsch(P,Q):
    Pc=P-P.mean(0);Qc=Q-Q.mean(0);V,S,Wt=np.linalg.svd(Pc.T@Qc)
    d=np.sign(np.linalg.det(V@Wt));R=V@np.diag([1,1,d])@Wt;return Pc@R+Q.mean(0)
best=None
for order in [chs[i:]+chs[:i] for i in range(5)]+[(chs[::-1])[i:]+(chs[::-1])[:i] for i in range(5)]:
    Xj=np.array([jqq[order[k]][r] for k in range(5) for r in res]);Xal=kabsch(Xj,Xa)
    rmsd=np.sqrt(((Xal-Xa)**2).sum(1).mean())
    if best is None or rmsd<best[0]: best=(rmsd,order,Xal)
rmsd,order,Xj=best
# ring normal (apo)
cent=np.array([Xa[i*n:(i+1)*n].mean(0) for i in range(5)]);ctr=cent.mean(0)
_,_,Vt=np.linalg.svd(cent-ctr);normal=Vt[2]
print(f"# apo<->7JQQ RMSD {rmsd:.1f}A (7JQQ=3-ATP intermediate). chain map apo(ABCDE)->7JQQ{order}")
# (b) per-subunit displacement + 7JQQ per-interface R146
def gua(X,i): 
    return np.array([X[i*n+res.index(146)]]) if 146 in res else np.zeros((0,3))
print("\n(b) per-subunit apo->7JQQ CA displacement + 7JQQ trans-R146(donor->next):")
for i,c in enumerate(chs):
    disp=np.sqrt(((Xj[i*n:(i+1)*n]-Xa[i*n:(i+1)*n])**2).sum(1).mean())
    oop=abs((Xj[i*n:(i+1)*n].mean(0)-Xa[i*n:(i+1)*n].mean(0))@normal)
    # R146 of subunit i -> Walker-A (res24-31) of next subunit, in 7JQQ
    j=(i+1)%5
    r146=np.array([Xj[i*n+res.index(146)]]) if 146 in res else np.zeros((0,3))
    wa=np.array([Xj[j*n+res.index(r)] for r in range(24,32) if r in res])
    d146=np.sqrt(((r146[:,None]-wa[None])**2).sum(-1)).min() if len(r146) and len(wa) else np.nan
    print(f"  subunit {c}: displacement {disp:5.1f}A (out-of-plane {oop:5.1f}A) | 7JQQ R146->{chs[j]} = {d146:5.1f}A {'(ATP-engaged)' if d146<8 else '(OPEN seam/apo)'}")
# (a) occupancy-dependent linear response
N=len(Xa);cut=13.0;H=np.zeros((3*N,3*N))
for i in range(N):
    dd=Xa-Xa[i];r2=(dd*dd).sum(1);r2[i]=1e9
    for j in np.where(r2<cut*cut)[0]:
        e=Xa[j]-Xa[i];e/=np.linalg.norm(e);k=np.outer(e,e)
        H[3*i:3*i+3,3*i:3*i+3]+=k;H[3*i:3*i+3,3*j:3*j+3]-=k
w,V=sla.eigh(H,subset_by_index=[0,150])
Hinv=sum(np.outer(V[:,m],V[:,m])/w[m] for m in range(6,151))   # pseudo-inverse (soft modes)
Dtar=(Xj-Xa)                                                    # per-CA target displacement toward 7JQQ
def opening(subset):
    f=np.zeros(3*N)
    for i in subset:                                           # apply target-pull force to ATP-bound subunits
        for a in range(i*n,(i+1)*n): f[3*a:3*a+3]=Dtar[a]
    dx=(Hinv@f).reshape(N,3)
    cdisp=np.array([dx[i*n:(i+1)*n].mean(0) for i in range(5)])
    return np.sqrt((np.array([abs(cdisp[i]@normal) for i in range(5)])**2).mean())  # out-of-plane (helicity)
print("\n(a) occupancy-dependent ENM linear response (ring out-of-plane opening vs # ATP):")
print("  #ATP (adjacent):", "  ".join(f"{k}:{opening(list(range(k))):.2f}" for k in range(6)))
# pattern test at 3 ATP: adjacent vs spread (sequential vs concerted signature)
import itertools
print("  3-ATP patterns (which 3 subunits):")
for pat in [(0,1,2),(0,2,4),(0,1,3)]:
    print(f"     {pat}: opening {opening(list(pat)):.2f}")

# --- transition path (NEB-lite): ENM energy along apo->7JQQ interpolation ---
def enm_energy(X):
    E=0.0
    for i in range(N):
        dd=X-X[i]; r2=(dd*dd).sum(1); r2[i]=1e9
        for j in np.where(r2<cut*cut)[0]:
            if j>i:
                r0=np.linalg.norm(Xa[j]-Xa[i]); r=np.linalg.norm(X[j]-X[i]); E+=0.5*(r-r0)**2
    return E
print("\n(NEB-lite) ENM energy along the apo->7JQQ linear path (barrier shape):")
prof=[enm_energy(Xa+t*(Xj-Xa)) for t in np.linspace(0,1,9)]
prof=[p-prof[0] for p in prof]
print("  t: " + "  ".join(f"{t:.2f}:{e:6.0f}" for t,e in zip(np.linspace(0,1,9),prof)))
print(f"  barrier ~ {max(prof):.0f} (arb. ENM units); monotone-uphill={all(prof[i]<=prof[i+1]+1 for i in range(len(prof)-1))}")
