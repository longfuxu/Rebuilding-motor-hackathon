import numpy as np, matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt
mpl.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False,
    "axes.titlesize":9.5,"axes.titleweight":"bold","font.family":"DejaVu Sans","figure.dpi":150})
BLU="#2166ac"; ORN="#d6604d"; GRY="#b8b8b8"; TEAL="#1b7837"; GLD="#e08214"
fig,ax=plt.subplots(2,3,figsize=(11.5,6.6))
fig.suptitle("gp16 ring opening: an ATP-occupancy mechanism that structure predictors miss",
             fontsize=12,fontweight="bold",y=0.99)

# A: opening vs #ATP
a=ax[0,0]; x=[0,1,2,3,4,5]; y=[0,0.99,4.19,7.19,12.17,13.24]
a.plot(x,y,"-o",color=BLU,lw=2.2,ms=6,mfc="white",mec=BLU,mew=1.8)
a.set_title("A  Opening is progressive & cooperative"); a.set_xlabel("# subunits with ATP bound")
a.set_ylabel("ring opening (a.u.)"); a.set_xticks(x)
a.annotate("each ATP opens\nthe ring more",(3,7.19),(0.4,10.5),fontsize=8,color=BLU,
           arrowprops=dict(arrowstyle="->",color=BLU))

# B: 3-ATP pattern
b=ax[0,1]; labs=["adjacent\n(1-2-3)","spread\n(1-3-5)","(1-2-4)"]; vals=[7.19,5.97,9.83]
cols=[ORN,GRY,GRY]; bars=b.bar(labs,vals,color=cols,width=0.62,edgecolor="k",lw=0.6)
b.set_title("B  3-ATP: which subunits sets the state"); b.set_ylabel("ring opening (a.u.)")
for r,v in zip(bars,vals): b.text(r.get_x()+r.get_width()/2,v+0.15,f"{v:.1f}",ha="center",fontsize=8.5)
b.text(0,3.2,"7JQQ\nadjacent\n(sequential)",ha="center",fontsize=7.2,color="white",fontweight="bold",linespacing=1.15)
b.set_ylim(0,11)

# C: per-subunit out-of-plane displacement apo->7JQQ
c=ax[0,2]; su=["A","B","C","D","E"]; oop=[0.3,1.5,3.9,4.5,1.1]
cc=[GRY,GRY,GLD,GLD,GRY]; c.bar(su,oop,color=cc,edgecolor="k",lw=0.6,width=0.66)
c.set_title("C  Opening localizes (apo→7JQQ)"); c.set_ylabel("out-of-plane shift (Å)"); c.set_ylim(0,5.4)
c.set_xlabel("subunit"); c.text(0.05,5.15,"helical distortion at the ATP region (C,D)",ha="left",va="top",fontsize=7.4,color=GLD)

# D: NMA cumulative overlap
d=ax[1,0]; m=[2,5,10,20]; pct=[1,19,33,41]
d.plot(m,pct,"-o",color=TEAL,lw=2.2,ms=6,mfc="white",mec=TEAL,mew=1.8)
d.axhline(0.4,ls="--",color=GRY,lw=1.2); d.text(20,2.5,"random ≈0.4%",ha="right",fontsize=7.5,color=GRY)
d.set_title("D  Opening is soft-mode-enriched"); d.set_xlabel("# softest ANM modes")
d.set_ylabel("% of apo→7JQQ captured"); d.set_xticks(m)
d.text(10,33,"enriched but NOT\none clean mode",fontsize=7.6,color=TEAL,va="top")

# E: AF3 3-state planarity vs experiment
e=ax[1,1]; st=["apo","+ADP","+ATP","7JQQ\n(exp.)"]; pl=[0.08,0.07,0.01,1.80]
ce=[BLU,BLU,BLU,ORN]; bars=e.bar(st,pl,color=ce,edgecolor="k",lw=0.6,width=0.64)
e.set_title("E  Predictors can't open the ring"); e.set_ylabel("ring non-planarity (Å)")
for r,v in zip(bars,pl): e.text(r.get_x()+r.get_width()/2,v+0.04,f"{v:.2f}",ha="center",fontsize=8)
e.text(1,1.0,"AF3 (all templates OFF):\nclosed planar in every\nnucleotide state",ha="center",fontsize=7.6,color=BLU)
e.set_ylim(0,2.05)

# F: text takeaway panel
f=ax[1,2]; f.axis("off")
f.text(0.0,0.98,"The mechanism prediction",fontsize=10,fontweight="bold",va="top")
txt=("• ATP opens the ring progressively &\n  cooperatively (A), and the OPENING\n  DEPENDS ON WHICH subunits fire (B).\n\n"
     "• Sequential firing → an ADJACENT ATP\n  block (7JQQ pattern, opening 7.2);\n  concerted/spread → 6.0 or 9.8.\n\n"
     "• Static predictors give a closed planar\n  ring for apo/ADP/ATP alike (E) — the\n  opening is real physics they miss (D).\n\n"
     "→ A covalent addressable ring + single-\n  molecule tweezers can SET the ATP\n  pattern and read step/force → decide\n  sequential vs concerted.")
f.text(0.0,0.86,txt,fontsize=8.2,va="top",linespacing=1.25)
plt.tight_layout(rect=[0,0,1,0.96])
plt.savefig("gp16_design/outputs/figures/fig_mechanism.png",dpi=300,bbox_inches="tight")
plt.savefig("gp16_design/outputs/figures/fig_mechanism.pdf",bbox_inches="tight")
print("wrote gp16_design/outputs/figures/fig_mechanism.{png,pdf}")
