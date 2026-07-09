#!/usr/bin/env python3
"""Geometric pre-filter for RFdiffusion Mode-B connector backbones.

Each output is a single chain (res 1..N): subunit A = 1..327 (gp16 4-330),
connector = 328..(N-327), subunit B = (N-326)..N (gp16 4-330). Motif is 654 res
(327+327) held fixed by RFdiffusion, so connector_len = N - 654.

Backbone-only (GLY), so only CA-based proxies. Ranks by: no clash, connector
compactness, and preservation of the A->B trans interface (R146(A)->WalkerA(B)).
"""
import sys, glob, os
import numpy as np

MOTIF = 327               # residues per subunit (gp16 4-330)
# within-subunit-A output indices (res1=gp16 4): R146 -> 143, WalkerA 24-31 -> 21-28
R146_A = 146 - 4 + 1      # 143
WALKER_OFF = [r - 4 + 1 for r in range(24, 32)]  # 21..28


def ca_coords(path):
    d = {}
    for L in open(path):
        if L[:4] == "ATOM" and L[12:16].strip() == "CA":
            d[int(L[22:26])] = np.array([float(L[30:38]), float(L[38:46]), float(L[46:54])])
    return d


def rg(P):
    P = np.array(P); c = P.mean(0); return float(np.sqrt(((P - c) ** 2).sum(1).mean()))


def analyze(path):
    ca = ca_coords(path)
    N = max(ca)
    clen = N - 2 * MOTIF
    if clen < 1:
        return None
    A = list(range(1, MOTIF + 1))
    con = list(range(MOTIF + 1, MOTIF + clen + 1))
    B = list(range(N - MOTIF + 1, N + 1))
    Bmap = {gp: B[0] + (gp - 1) for gp in range(1, MOTIF + 1)}  # subunitB output idx by within-copy pos
    conP = [ca[r] for r in con]
    motifP = [ca[r] for r in A + B]
    # connector compactness: Rg per residue (lower = more compact/folded)
    rg_con = rg(conP) if len(conP) > 1 else 0.0
    rg_per = rg_con / clen
    # clash: min distance connector-CA to any motif-CA (exclude the 2 junction neighbours)
    conM = np.array(conP); motM = np.array(motifP)
    dcm = np.sqrt(((conM[:, None, :] - motM[None, :, :]) ** 2).sum(-1))
    min_clash = float(dcm.min())
    # A->B trans interface: R146(A) CA -> nearest WalkerA(B) CA
    r146 = ca[R146_A]
    waB = np.array([ca[Bmap[o]] for o in WALKER_OFF])
    iface = float(np.sqrt(((waB - r146) ** 2).sum(1)).min())
    # crude SS proxy: fraction of connector with CA(i)-CA(i+3) ~5-6A (helix) or extended
    ssh = 0
    for i in range(len(conP) - 3):
        d = np.linalg.norm(conP[i] - conP[i + 3])
        if 4.5 < d < 6.5:
            ssh += 1
    helix_frac = ssh / max(1, len(conP) - 3)
    return dict(name=os.path.basename(path), N=N, clen=clen, rg_con=round(rg_con, 1),
                rg_per=round(rg_per, 3), min_clash=round(min_clash, 2),
                iface_CA=round(iface, 2), helix_frac=round(helix_frac, 2))


def main():
    paths = sorted(glob.glob(sys.argv[1]))
    rows = [r for r in (analyze(p) for p in paths) if r]
    # rank: no clash (>=3.5A backbone), then compact connector, then some SS
    def key(r):
        clash_ok = r["min_clash"] >= 3.5
        iface_ok = r["iface_CA"] <= 12.0
        return (not clash_ok, not iface_ok, r["rg_per"], -r["helix_frac"])
    rows.sort(key=key)
    print(f"{'name':<22}{'N':>5}{'clen':>5}{'rg_con':>7}{'rg/res':>7}{'clash':>7}{'iface':>7}{'helix':>6}  rank")
    for i, r in enumerate(rows):
        flag = "OK" if (r["min_clash"] >= 3.5 and r["iface_CA"] <= 12.0) else "--"
        print(f"{r['name']:<22}{r['N']:>5}{r['clen']:>5}{r['rg_con']:>7}{r['rg_per']:>7}"
              f"{r['min_clash']:>7}{r['iface_CA']:>7}{r['helix_frac']:>6}  {flag}")
    import json
    json.dump(rows, open(os.path.join(os.path.dirname(sys.argv[1].rstrip('*')) or '.', 'geom_rank.json'), 'w'), indent=2)


if __name__ == "__main__":
    main()
