import mdtraj as md, numpy as np, json, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa

REPO="/Users/longfu/Developer/claude-science-hackthon/.claude/worktrees/export-gp16-results/gp16_design"
OUT="/Users/longfu/.claude/jobs/48b7a543/tmp/bioemu_analysis"; os.makedirs(OUT, exist_ok=True)

t   = md.load(f"{REPO}/outputs/bioemu/cp233_apo/samples.xtc", top=f"{REPO}/outputs/bioemu/cp233_apo/topology.pdb")
ref = md.load(f"{REPO}/md/openmm_validation/trajectories/C/C_start.pdb")   # folded design (MD start)
ca_t   = t.topology.select("name CA")
ca_ref = ref.topology.select("name CA")
print("samples frames", t.n_frames, "nCA", len(ca_t), "| ref nCA", len(ca_ref))

# 5 subunits: residues 1-342, 353-694, 705-1046, 1057-1398, 1409-1750 (0-indexed CA blocks of 342, stride 352)
SUB=[(352*k, 352*k+342) for k in range(5)]
def rg(x):  # x:(n,3) Å -> nm
    c=x.mean(0); return float(np.sqrt(((x-c)**2).sum(1).mean())/10.0)

xyz   = t.xyz[:, ca_t, :]*10.0        # frames,1750,3 in Å
xref  = ref.xyz[0, ca_ref, :]*10.0    # 1750,3 Å
n=min(xyz.shape[1], xref.shape[0]); xyz=xyz[:,:n]; xref=xref[:n]

def kabsch_rmsd(P,Q):
    Pc=P-P.mean(0); Qc=Q-Q.mean(0)
    V,S,Wt=np.linalg.svd(Pc.T@Qc); d=np.sign(np.linalg.det(V@Wt)); D=np.diag([1,1,d])
    U=V@D@Wt; return float(np.sqrt(((Pc@U-Qc)**2).sum(1).mean()))

rows=[]
for f in range(t.n_frames):
    x=xyz[f]
    tot_rg=rg(x)
    sub_rg=[rg(x[a:b]) for a,b in SUB]
    sub_rmsd=[kabsch_rmsd(x[a:b], xref[a:b]) for a,b in SUB]
    whole_rmsd=kabsch_rmsd(x, xref)
    rows.append(dict(frame=f, total_rg=tot_rg, sub_rg=sub_rg, sub_rmsd=sub_rmsd, whole_rmsd=whole_rmsd))

ref_sub_rg=[rg(xref[a:b]) for a,b in SUB]; ref_tot_rg=rg(xref)
print("REF folded: total Rg %.2f nm, per-subunit Rg %s"%(ref_tot_rg,[round(r,2) for r in ref_sub_rg]))
for r in rows:
    print("frame %d: totRg %.2f | subRg %s | subRMSD(Å) %s | wholeRMSD %.1f"%(
        r["frame"], r["total_rg"], [round(v,2) for v in r["sub_rg"]],
        [round(v,1) for v in r["sub_rmsd"]], r["whole_rmsd"]))

summary=dict(ref_total_rg=ref_tot_rg, ref_sub_rg=ref_sub_rg,
             sample_total_rg=[r["total_rg"] for r in rows],
             sample_sub_rg=[r["sub_rg"] for r in rows],
             sample_sub_rmsd=[r["sub_rmsd"] for r in rows],
             sample_whole_rmsd=[r["whole_rmsd"] for r in rows],
             mono_folded_rg=2.29)
json.dump(summary, open(f"{OUT}/metrics.json","w"), indent=1)

# ---- Plot 1: total Rg distribution vs references ----
fig,ax=plt.subplots(figsize=(6,3.4))
srg=[r["total_rg"] for r in rows]
ax.scatter(srg, np.random.RandomState(0).uniform(-.1,.1,len(srg)), s=60, c="#c0392b", zorder=3, label="BioEmu samples (10)")
for val,name,col in [(3.68,"folded design","#27ae60"),(3.54,"native apo","#2980b9"),(3.84,"7JQQ helical","#8e44ad")]:
    ax.axvline(val, color=col, ls="--", lw=1.6, label=f"{name} {val}nm")
ax.set_xlabel("radius of gyration Rg (nm)"); ax.set_yticks([]); ax.set_xlim(3,9.2)
ax.set_title("cp233 BioEmu ensemble: total Rg vs folded references"); ax.legend(fontsize=7, loc="upper right")
plt.tight_layout(); plt.savefig(f"{OUT}/fig_rg.png", dpi=130); plt.close()

# ---- Plot 2: per-subunit Rg (are the 5 domains individually folded?) ----
fig,ax=plt.subplots(figsize=(6,3.4))
arr=np.array([r["sub_rg"] for r in rows])  # 10x5
for k in range(5): ax.scatter([k+1]*len(arr), arr[:,k], s=30, c="#c0392b", alpha=.7)
ax.plot(range(1,6), ref_sub_rg, "o-", c="#27ae60", label="folded design subunits")
ax.axhline(2.29, color="#7f8c8d", ls=":", label="isolated gp16 domain (2.29nm)")
ax.set_xlabel("subunit (of 5)"); ax.set_ylabel("per-subunit Rg (nm)")
ax.set_title("Are the individual domains folded?"); ax.legend(fontsize=7); ax.set_xticks(range(1,6))
plt.tight_layout(); plt.savefig(f"{OUT}/fig_subunit_rg.png", dpi=130); plt.close()

# ---- Plot 3: 3D CA traces colored by subunit — most compact sample vs folded ref ----
order=np.argsort(srg); f_compact=int(order[0]); f_expand=int(order[-1])
cols=["#e74c3c","#f39c12","#27ae60","#2980b9","#8e44ad"]
def plot3d(ax, X, title):
    for k,(a,b) in enumerate(SUB):
        seg=X[a:b]; ax.plot(seg[:,0],seg[:,1],seg[:,2], c=cols[k], lw=1)
    ax.set_title(title, fontsize=9); ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
fig=plt.figure(figsize=(10,3.6))
ax1=fig.add_subplot(131,projection="3d"); plot3d(ax1, xref, "folded design (ref)\nRg 3.68nm")
ax2=fig.add_subplot(132,projection="3d"); plot3d(ax2, xyz[f_compact], f"BioEmu most-compact\nRg {srg[f_compact]:.1f}nm")
ax3=fig.add_subplot(133,projection="3d"); plot3d(ax3, xyz[f_expand], f"BioEmu most-expanded\nRg {srg[f_expand]:.1f}nm")
plt.tight_layout(); plt.savefig(f"{OUT}/fig_3d.png", dpi=130); plt.close()
print("WROTE", OUT)
PY_DONE=1
print("DONE")
