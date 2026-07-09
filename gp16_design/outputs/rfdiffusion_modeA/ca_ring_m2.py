#!/usr/bin/env python3
"""CA-based sequential R146->WalkerA for a CP cp233 pentamer (backbone-only OK).
cp233_int15_inter10: copies start [1,353,705,1057,1409], each 342 aa;
within-copy R146=255, Walker-A=133-140. Reports per-sequential-interface
R146(CA)->min WalkerA(CA) distance for the 5 copies (designed cyclic order)."""
import sys, numpy as np

STARTS = [1, 353, 705, 1057, 1409]
R146_INCOPY = 255
WA_INCOPY = list(range(133, 141))


def ca(path):
    d = {}
    if path.lower().endswith((".cif", ".mmcif")):
        order = {}; cols = {}; inl = False
        for L in open(path):
            s = L.strip()
            if s.startswith("_atom_site."):
                order[len(order)] = s.split(".", 1)[1]; cols = {v: k for k, v in order.items()}; inl = True; continue
            if inl and s.startswith("ATOM"):
                f = s.split()
                if len(f) < len(order): continue
                if f[cols["label_atom_id"]].strip('"') != "CA": continue
                rn = int(f[cols.get("auth_seq_id", cols.get("label_seq_id"))])
                d[rn] = np.array([float(f[cols["Cartn_x"]]), float(f[cols["Cartn_y"]]), float(f[cols["Cartn_z"]])])
            elif inl and s == "#": inl = False
    else:
        for L in open(path):
            if L[:4] == "ATOM" and L[12:16].strip() == "CA":
                d[int(L[22:26])] = np.array([float(L[30:38]), float(L[38:46]), float(L[46:54])])
    return d


def analyze(path):
    c = ca(path)
    dists = []
    for i in range(5):
        donor = STARTS[i]; acc = STARTS[(i + 1) % 5]
        r = c.get(donor - 1 + R146_INCOPY)
        wa = [c[acc - 1 + w] for w in WA_INCOPY if (acc - 1 + w) in c]
        if r is None or not wa:
            dists.append(None); continue
        dists.append(float(min(np.linalg.norm(r - w) for w in wa)))
    ok = [d for d in dists if d is not None]
    return dists, (np.mean(ok) if ok else None)


if __name__ == "__main__":
    for p in sys.argv[1:]:
        dists, mean = analyze(p)
        ds = " ".join(f"{d:.1f}" if d else "NA" for d in dists)
        eng = sum(1 for d in dists if d and d < 12.0)  # CA-based looser cutoff
        print(f"{p.split('/')[-1]:<42} R146->WA CA: [{ds}]  mean {mean:.1f}  <12A:{eng}/5" if mean else f"{p}: NA")
