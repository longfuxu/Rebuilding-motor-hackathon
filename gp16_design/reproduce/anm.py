import numpy as np, scipy.linalg as sla
PDB="gp16_design/outputs/structures/cycle3/native_ring__boltz__s1__native_ring_model_0.pdb"
xyz=[]; ch=[]
for L in open(PDB):
    if L[:4]=="ATOM" and L[12:16].strip()=="CA":
        xyz.append([float(L[30:38]),float(L[38:46]),float(L[46:54])]); ch.append(L[21])
X=np.array(xyz); ch=np.array(ch); N=len(X)
# --- Anisotropic Network Model (elastic network): springs between CA within cutoff ---
cut=13.0; H=np.zeros((3*N,3*N))
for i in range(N):
    d=X-X[i]; r2=(d*d).sum(1); r2[i]=1e9
    js=np.where(r2<cut*cut)[0]
    for j in js:
        e=(X[j]-X[i]); e/= np.linalg.norm(e); k=np.outer(e,e)
        H[3*i:3*i+3,3*i:3*i+3]+=k; H[3*i:3*i+3,3*j:3*j+3]-=k
w,V=sla.eigh(H, subset_by_index=[0,19])   # lowest 20 eigenpairs
# ring frame from subunit centroids
chs=sorted(set(ch)); cent=np.array([X[ch==c].mean(0) for c in chs]); ctr=cent.mean(0)
_,_,Vt=np.linalg.svd(cent-ctr); normal=Vt[2]
rad=[(p-ctr)-((p-ctr)@normal)*normal for p in cent]; rad=[r/np.linalg.norm(r) for r in rad]
print(f"# ANM on native gp16 ring ({N} CA, cutoff {cut}A). Modes 1-6 = rigid body (skip).")
print(f"{'mode':>5} {'softness(1/eig)':>15} {'helical/out-of-plane':>21} {'radial-breathing':>17} {'localization':>13}  character")
soft0=None
for m in range(6,20):
    v=V[:,m].reshape(N,3)
    dc=np.array([v[ch==c].mean(0) for c in chs])           # per-subunit centroid displacement
    mag=np.linalg.norm(dc,axis=1)
    oop=np.array([abs(dc[i]@normal) for i in range(5)])
    radc=np.array([ (dc[i]-(dc[i]@normal)*normal)@rad[i] for i in range(5)])
    helical=np.sqrt((oop**2).mean())/ (np.sqrt((mag**2).mean())+1e-9)
    breathing=abs(radc.mean())/(np.sqrt((mag**2).mean())+1e-9)
    local=mag.max()/(mag.mean()+1e-9)
    eig=w[m]; soft=1.0/eig if eig>1e-6 else 9e9
    if soft0 is None: soft0=soft
    lab=("HELICAL/opening" if helical>0.5 else ("breathing(radial)" if breathing>0.5 else "in-plane/other"))
    if local>2.2: lab+=" +localized(seam)"
    print(f"{m+1:>5} {soft/soft0:>14.2f}x {helical:>21.2f} {breathing:>17.2f} {local:>13.2f}  {lab}")
