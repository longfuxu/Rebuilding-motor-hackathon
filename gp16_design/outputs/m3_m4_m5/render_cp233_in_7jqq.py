#!/usr/bin/env python3
"""Render cp233 design (Ca cartoon) placed in the 7JQQ frame with the 7JQQ dsDNA threading its
pore and the pRNA below, DNA-contact residues highlighted facing the channel.
matplotlib 3D fallback (no PyMOL/ChimeraX available). -> cp233_in_7jqq.png"""
import warnings; warnings.filterwarnings("ignore")
import pickle, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection

OUT="/Users/longfu/Developer/claude-science-hackthon/.claude/worktrees/export-gp16-results/gp16_design/outputs/m3_m4_m5/"
D=pickle.load(open(OUT+"_placements.pkl","rb"))
P=D["cp233"]
coord=P["coord"]; nres=P["native_resid"]; copy=P["copy"]; link=P["islink"]; anm=P["atom_name"]; rid=P["res_id"]
DNA=D["_DNA"]; DNA_atom=D["_DNA_atom"]; DNA_chain=D["_DNA_chain"]
PRNA=D["_PRNA"]; PRNA_atom=D["_PRNA_atom"]; PRNA_chain=D["_PRNA_chain"]
DNAC=set(D["_DNA_CONTACTS"])

# ---- orient: DNA principal axis -> z, DNA channel center -> origin
dc=DNA.mean(0); u,s,vt=np.linalg.svd(DNA-dc,full_matrices=False); axis=vt[0]
if axis[2]<0: axis=-axis
# build rotation mapping axis->z
z=np.array([0,0,1.0]); v=np.cross(axis,z); c=axis@z
if np.linalg.norm(v)<1e-6: R=np.eye(3)
else:
    vx=np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
    R=np.eye(3)+vx+vx@vx*(1/(1+c))
def T(x): return (x-dc)@R.T
coordT=T(coord); DNAT=T(DNA); PRNAT=T(PRNA)
# focus on the gp16 ring z-band: use design Ca z-range
caz=coordT[anm=="CA"][:,2]
zlo,zhi=np.percentile(caz,1)-8, np.percentile(caz,99)+8
def band(x): return (x[:,2]>zlo)&(x[:,2]<zhi)

COL=["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd"]  # 5 subunits

fig=plt.figure(figsize=(15,7.2))
plt.rcParams.update({"font.size":9})

def draw_backbone(ax, subunit_alpha=0.9, lw=1.4):
    # per copy, per contiguous native segment, connect Ca in res_id order
    ca_mask=(anm=="CA")&(~link)
    for c in range(5):
        m=ca_mask&(copy==c)
        r=rid[m]; xyz=coordT[m]; order=np.argsort(r); r=r[order]; xyz=xyz[order]
        # split where res_id jumps (linker gaps) -> separate polyline segments
        segs=[]; cur=[0]
        for i in range(1,len(r)):
            if r[i]-r[i-1]>1: segs.append(cur); cur=[i]
            else: cur.append(i)
        segs.append(cur)
        for seg in segs:
            if len(seg)<2: continue
            p=xyz[seg]
            ax.plot(p[:,0],p[:,1],p[:,2],color=COL[c],lw=lw,alpha=subunit_alpha,solid_capstyle="round")

def draw_dna(ax, lw=3.5):
    for ch,col in [("F","#111111"),("G","#555555")]:
        m=(DNA_chain==ch)&np.isin(DNA_atom,["P","C1'","O5'"])
        # order along axis (z)
        p=DNAT[m]; p=p[np.argsort(p[:,2])]
        ax.plot(p[:,0],p[:,1],p[:,2],color=col,lw=lw,alpha=0.95,solid_capstyle="round")

def draw_prna(ax, lw=1.2):
    for ch in ["K","L","M","N","O"]:
        m=(PRNA_chain==ch)&np.isin(PRNA_atom,["P","C1'"])
        p=PRNAT[m]; p=p[np.argsort(p[:,2])]
        ax.plot(p[:,0],p[:,1],p[:,2],color="#7fb3d5",lw=lw,alpha=0.35)

def draw_contacts(ax, s=26):
    # DNA-contact residue sidechain tips (closest sidechain atom to a DNA atom would be ideal;
    # here plot Ca of contact residues, colored by domain)
    from scipy.spatial import cKDTree
    tree=cKDTree(DNAT)
    for c in range(5):
        for nr in DNAC:
            m=(nres==nr)&(copy==c)&(~link)&(~np.isin(anm,["N","C","O"]))
            if not m.any(): continue
            sc=coordT[m]; tip=sc[np.argmin(tree.query(sc,k=1)[0])]
            dom = "#e41a1c" if nr<=205 else "#ff9500"   # N-dom red, C-dom orange
            ax.scatter(*tip,color=dom,s=s,edgecolor="k",linewidth=0.3,depthshade=True,zorder=5)

# ---------- Panel A: side view (DNA threading the pore)
axA=fig.add_subplot(1,2,1,projection="3d")
draw_prna(axA); draw_backbone(axA); draw_dna(axA); draw_contacts(axA)
axA.set_title("Side view — 7JQQ dsDNA threads the cp233 ring pore\n(pRNA faint blue below; DNA-contact residues: red=ATPase-domain, orange=C-domain)",fontsize=9)
axA.view_init(elev=8,azim=-60)
axA.set_box_aspect((1,1,1.6))

# ---------- Panel B: top view down the channel axis
axB=fig.add_subplot(1,2,2,projection="3d")
# only ring band for clarity
draw_backbone(axB);
# DNA cross section in band
for ch,col in [("F","#111111"),("G","#555555")]:
    m=(DNA_chain==ch)&np.isin(DNA_atom,["P","C1'","O5'"]); p=DNAT[m]; p=p[band(p)]
    axB.scatter(p[:,0],p[:,1],p[:,2],color=col,s=8,alpha=0.9)
draw_contacts(axB)
axB.set_title("Axial view (down channel) — contact residues face the central DNA\nfive subunits color-coded; central dark points = dsDNA backbone",fontsize=9)
axB.view_init(elev=88,azim=-90)
axB.set_box_aspect((1,1,0.4))

for ax in (axA,axB):
    ax.set_xlabel("x (Å)"); ax.set_ylabel("y (Å)"); ax.set_zlabel("channel axis (Å)")
    ax.grid(False)
    try: ax.set_zlim(zlo,zhi)
    except Exception: pass

fig.suptitle("cp233 single-chain gp16 ring placed in the 7JQQ motor frame (ATPase-domain superposition, RMSD 2.1 Å)\n"
             "dsDNA (F,G) threads the design pore; ATPase-domain DNA contacts line the channel — geometry only, no pLDDT",
             fontsize=10.5, y=0.99)
fig.tight_layout(rect=(0,0,1,0.94))
fig.savefig(OUT+"cp233_in_7jqq.png",dpi=170)
print("wrote", OUT+"cp233_in_7jqq.png")
