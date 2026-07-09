#!/usr/bin/env python3
"""Reusable: generate WT / one-seat / all-seat point-mutant sequences of a tiled single-chain
ring, for AF3/Boltz functional-locality tests. Locates the target native residue in every copy
by sequence anchor, mutates, writes FASTAs.

Usage: python make_deadseat_variants.py <construct.txt> <native_ref.pdb> <native_resid> <to_aa> <outdir>
Example (cp233 E119Q): ... cp233_int15_inter10_AF3_sequence.txt A_apo.pdb 119 Q ~/Desktop/out
"""
import sys, os
import biotite.structure.io.pdb as pdb
from biotite.sequence import ProteinSequence
def one(r):
    try: return ProteinSequence.convert_letter_3to1(r)
    except: return "X"
def findall(h,n):
    o=[];i=h.find(n)
    while i!=-1:o.append(i);i=h.find(n,i+1)
    return o
def main(construct,nativepdb,resid,toaa,outdir):
    resid=int(resid)
    a=pdb.PDBFile.read(nativepdb).get_structure(model=1)
    ch=sorted(set(a.chain_id))[0]
    nseq="".join(one(r) for r in a[(a.chain_id==ch)&(a.atom_name=="CA")].res_name)
    dseq=open(construct).read().strip()
    anchor=nseq[resid-6:resid+5]                 # 11-mer around the target
    hits=findall(dseq,anchor); pos=[h+5 for h in hits]  # target within each copy
    assert hits, f"anchor {anchor!r} not found in construct"
    wt=nseq[resid-1]
    print(f"native res {resid}={wt}; found {len(pos)} copies at {[p+1 for p in pos]} (all {[dseq[p] for p in pos]})")
    def mut(seq,ps):
        s=list(seq)
        for p in ps: s[p]=toaa
        return "".join(s)
    os.makedirs(outdir,exist_ok=True)
    variants={f"WT":(dseq,"baseline"),
              f"{wt}{resid}{toaa}_1seat":(mut(dseq,pos[:1]),"one dead seat"),
              f"{wt}{resid}{toaa}_{len(pos)}seat":(mut(dseq,pos),"all seats knockout")}
    for nm,(sq,d) in variants.items():
        open(f"{outdir}/{nm}.fasta","w").write(f">{nm} {d}\n{sq}\n")
        print("wrote",nm)
if __name__=="__main__":
    main(*sys.argv[1:6])
