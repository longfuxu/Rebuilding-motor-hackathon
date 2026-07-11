"""C3 (refined) — Cα structure-based (Go) model of the gp16 ring P->H transition, OpenMM.
Screen: is the planar(P)->helical(H) rearrangement concerted or sequential?

Refinements vs v1 (which was inconclusive: Q_H started at 0.60 because P&H share ~60% of
contacts, and the union dual-basin was frustrated / stuck at 0.68):
  - readout uses ONLY H-UNIQUE contacts (present in H, absent in P) -> starts ~0 in P, ->1 in H,
    so it actually reports the inter-subunit rearrangement, not shared intra-domain contacts.
  - single-basin-on-P + strong slow steer toward H (RMSD CV) -> drives cleanly to H (no union
    frustration); the ORDER subunits form their H-unique contacts = concerted vs sequential.
  - per-subunit crossing window (Q=0.5) reported -> concerted: all near same window; sequential: spread.

Usage:
  python c3_sbm.py test
  python c3_sbm.py steer --steps 6000000 --ksteer 8000
"""
import sys, argparse, os, json, numpy as np
import openmm as mm
from openmm import unit
HERE=os.path.dirname(os.path.abspath(__file__))
SBM=os.environ.get("SBM_DIR", os.path.join(HERE,"..","..","outputs","php_cycle","C3_sbm"))
KJ=unit.kilojoule_per_mole; NM=unit.nanometer

def load_ca(pdb):
    xyz=[]; ch=[]
    for l in open(pdb):
        if l.startswith("ATOM") and l[12:16].strip()=="CA":
            xyz.append([float(l[30:38]),float(l[38:46]),float(l[46:54])]); ch.append(l[21])
    return np.array(xyz)/10.0, ch

def contacts(xyz, chain, cut=1.0, sep=4):
    from scipy.spatial import cKDTree
    ch=np.array(chain); pairs={}
    for i,j in cKDTree(xyz).query_pairs(cut):
        a,b=(i,j) if i<j else (j,i)
        if ch[a]==ch[b] and (b-a)<sep: continue
        pairs[(a,b)]=float(np.linalg.norm(xyz[a]-xyz[b]))
    return pairs

def build(xyzP, chain, contact_map, kbond=20000., kang=40., eps_c=1.0):
    n=len(xyzP); s=mm.System()
    for _ in range(n): s.addParticle(110.0)
    bond=mm.HarmonicBondForce(); bonded=set()
    for i in range(n-1):
        if chain[i]==chain[i+1]:
            r0=np.linalg.norm(xyzP[i]-xyzP[i+1]); bond.addBond(i,i+1,r0*NM,kbond*KJ/NM**2); bonded.add((i,i+1))
    s.addForce(bond)
    ang=mm.HarmonicAngleForce()
    for i in range(n-2):
        if chain[i]==chain[i+1]==chain[i+2]:
            a=xyzP[i]-xyzP[i+1]; b=xyzP[i+2]-xyzP[i+1]
            th=np.arccos(np.clip(a.dot(b)/(np.linalg.norm(a)*np.linalg.norm(b)),-1,1))
            ang.addAngle(i,i+1,i+2,th*unit.radian,kang*KJ/unit.radian**2)
    s.addForce(ang)
    go=mm.CustomBondForce("eps*(5*(r0/r)^12-6*(r0/r)^10)")
    go.addPerBondParameter("r0"); go.addPerBondParameter("eps"); cset=set()
    for (i,j),r0 in contact_map.items(): go.addBond(i,j,[r0*NM,eps_c*KJ]); cset.add((i,j))
    s.addForce(go)
    ev=mm.CustomNonbondedForce("epsev*(sig/r)^12")
    ev.addGlobalParameter("epsev",1.0); ev.addGlobalParameter("sig",0.4)
    ev.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic); ev.setCutoffDistance(1.5*NM)
    for _ in range(n): ev.addParticle([])
    for (i,j) in bonded|cset: ev.addExclusion(i,j)
    s.addForce(ev)
    return s

def per_sub_Q(pos, cset, chain, tol=1.2):
    per={c:[0,0] for c in sorted(set(chain))}
    for (i,j),r0 in cset.items():
        per[chain[i]][1]+=1
        if np.linalg.norm(pos[i]-pos[j])<tol*r0: per[chain[i]][0]+=1
    return {c:(v[0]/v[1] if v[1] else 0.) for c,v in per.items()}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("mode",choices=["test","steer"])
    ap.add_argument("--steps",type=int,default=6000000); ap.add_argument("--ksteer",type=float,default=8000.)
    ap.add_argument("--Tred",type=float,default=100.); a=ap.parse_args()
    xyzP,chain=load_ca(os.path.join(SBM,"P_matched.pdb"))
    xyzH,_=load_ca(os.path.join(SBM,"H_matched.pdb"))
    cP=contacts(xyzP,chain); cH=contacts(xyzH,chain)
    Huniq={k:v for k,v in cH.items() if k not in cP}   # H-unique = the rearrangement signal
    print(f"beads {len(xyzP)}  P-contacts {len(cP)}  H-contacts {len(cH)}  H-UNIQUE {len(Huniq)}",flush=True)
    s=build(xyzP,chain,cP)                              # single-basin on P
    rmsd=mm.RMSDForce((xyzH*NM),list(range(len(xyzH))))
    cv=mm.CustomCVForce("0.5*ksteer*(rmsd-r0)^2"); cv.addCollectiveVariable("rmsd",rmsd)
    cv.addGlobalParameter("ksteer",a.ksteer); cv.addGlobalParameter("r0",100.); s.addForce(cv)
    integ=mm.LangevinMiddleIntegrator(a.Tred*unit.kelvin,1./unit.picosecond,0.010*unit.picoseconds)
    names=[mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
    plat=mm.Platform.getPlatformByName(next(p for p in ["CUDA","OpenCL","CPU"] if p in names))
    ctx=mm.Context(s,integ,plat); ctx.setPositions(xyzP*NM); print("platform",plat.getName(),flush=True)
    mm.LocalEnergyMinimizer.minimize(ctx,maxIterations=2000)
    if a.mode=="test":
        ctx.setParameter("r0",4.0); integ.step(20000)
        st=ctx.getState(getPositions=True,getEnergy=True)
        pos=st.getPositions(asNumpy=True).value_in_unit(NM)
        print("test 20k: E=%.0f finite=%s  Q_Huniq=%s"%(st.getPotentialEnergy().value_in_unit(KJ),
              np.isfinite(pos).all(), {k:round(v,2) for k,v in per_sub_Q(pos,Huniq,chain).items()}),flush=True)
        return
    # steer: ramp r0 4.3->0.2 nm, log per-subunit Q_Huniq + crossing window
    nwin=60; nstep=a.steps//nwin; log=[]; cross={c:None for c in sorted(set(chain))}
    for w in range(nwin+1):
        r0=4.3*(1-w/nwin)+0.2*(w/nwin); ctx.setParameter("r0",r0); integ.step(nstep)
        pos=ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(NM)
        q=per_sub_Q(pos,Huniq,chain)
        for c in q:
            if cross[c] is None and q[c]>=0.5: cross[c]=w
        log.append({"win":w,"r0":round(r0,2),"Q":{k:round(v,2) for k,v in q.items()},
                    "Qmean":round(float(np.mean(list(q.values()))),3)})
        if w%5==0: print(f"win {w:2d} r0={r0:.2f} Qmean={log[-1]['Qmean']} Q={log[-1]['Q']}",flush=True)
        json.dump({"log":log,"cross":cross},open(os.path.join(SBM,"steer_refined.json"),"w"),indent=1)
    order=sorted(cross.items(), key=lambda kv:(kv[1] is None, kv[1]))
    spread=[w for w in cross.values() if w is not None]
    print("CROSSING WINDOW per subunit (Q_Huniq>=0.5):",cross,flush=True)
    print("ORDER:",[c for c,_ in order],"  window-spread=",(max(spread)-min(spread)) if len(spread)>1 else "n/a",flush=True)
    print("=> spread small => CONCERTED ; spread large / staggered => SEQUENTIAL",flush=True)
    print("STEER DONE",flush=True)

if __name__=="__main__": main()
