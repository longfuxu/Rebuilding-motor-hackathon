#!/usr/bin/env python3
"""Build gp17 single-chain constructs (ATPase-domain ring) + write subunit fasta for MSA.

gp17 = T4 large terminase (UniProt P17312, 610 aa). Functional motor ring = N-terminal
ATPase domain (res 1-360; 2O0H covers 1-357). We build the ring from the ATPase domain,
the ClpX(deltaN) analog: the part that carries Walker-A(161-167,K166), Walker-B(251-256,
D255), catalytic E256, arg-finger candidate R162, and forms the DNA channel.
"""
import json, os

# UniProt P17312 (TERL_BPT4) full sequence
FULL = ("MEQPINVLNDFHPLNEAGKILIKHPSLAERKDEDGIHWIKSQWDGKWYPEKFSDYLRLHK"
        "IVKIPNNSDKPELFQTYKDKNNKRSRYMGLPNLKRANIKTQWTREMVEEWKKCRDDIVYF"
        "AETYCAITHIDYGVIKVQLRDYQRDMLKIMSSKRMTVCNLSRQLGKTTVVAIFLAHFVCF"
        "NKDKAVGILAHKGSMSAEVLDRTKQAIELLPDFLQPGIVEWNKGSIELDNGSSIGAYASS"
        "PDAVRGNSFAMIYIDECAFIPNFHDSWLAIQPVISSGRRSKIIITTTPNGLNHFYDIWTA"
        "AVEGKSGFEPYTAIWNSVKERLYNDEDIFDDGWQWSIQTINGSSLAQFRQEHTAAFEGTS"
        "GTLISGMKLAVMDFIEVTPDDHGFHQFKKPEPDRKYIATLDCSEGRGQDYHALHIIDVTD"
        "DVWEQVGVLHSNTISHLILPDIVMRYLVEYNECPVYIELNSTGVSVAKSLYMDLEYEGVI"
        "CDSYTDLGMKQTKRTKAVGCSTLKDLIEKDKLIIHHRATIQEFRTFSEKGVSWAAEEGYH"
        "DDLVMSLVIFGWLSTQSKFIDYADKDDMRLASEVFSKELQDMSDDYAPVIFVDSVHSAEY"
        "VPVSHGMSMV")
assert len(FULL) == 610, len(FULL)

ATP_LO, ATP_HI = 1, 360                      # ATPase domain (1-based inclusive)
SUB = FULL[ATP_LO-1:ATP_HI]                  # 360-aa subunit
assert len(SUB) == 360

# sanity: landmark residues (1-based in SUB == native numbering here)
assert SUB[166-1] == 'K', SUB[166-1]          # Walker-A K166
assert SUB[162-1] == 'R', SUB[162-1]          # arg-finger candidate R162
assert SUB[255-1] == 'D' and SUB[256-1] == 'E'  # Walker-B D255 / catalytic E256
assert SUB[245-1] == 'R'                        # alt arg-finger R245

OUT = os.path.dirname(os.path.abspath(__file__))
LINK20 = "GGGGS" * 4       # 20-aa flexible linker (spans the 34 A C->N gap)

def pentamer(sub, linker):
    return linker.join([sub]*5), [1 + i*(len(sub)+len(linker)) for i in range(5)]

# --- construct 1: direct C->N fusion, 20-aa linker ---
df_seq, df_starts = pentamer(SUB, LINK20)
constructs = {
  "gp17_atpase_monomer": {"seq": SUB, "note": "ATPase-domain monomer (res1-360), positive control for fold+MSA"},
  "gp17_directfusion_L20": {"seq": df_seq, "linker": LINK20, "copy_starts": df_starts,
        "copy_len": len(SUB), "note": "direct C(res360)->N(res1) fusion x5, 20-aa GS linker"},
}
print("subunit len:", len(SUB))
for name, c in constructs.items():
    print(f"{name}: {len(c['seq'])} aa" + (f" | copy_starts {c.get('copy_starts')}" if 'copy_starts' in c else ""))

# write fastas
with open(f"{OUT}/seqs/gp17_subunit.fasta","w") as f:
    f.write(f">gp17_atpase_res1-360\n{SUB}\n")
for name, c in constructs.items():
    with open(f"{OUT}/seqs/{name}.fasta","w") as f:
        f.write(f">{name}\n{c['seq']}\n")
json.dump(constructs, open(f"{OUT}/seqs/constructs.json","w"), indent=2)
print("wrote seqs/ + constructs.json")
