#!/usr/bin/env python3
"""Physics-independent (elastic-network) open/close test: does the covalent
single-chain cp233 retain the native ring's soft ring-opening/breathing motion,
or does the linkage over-rigidify it? No learned priors -- pure ENM (springs).

Metric: cumulative overlap of the ring-breathing collective coordinate with the
lowest-K ENM modes. High = the opening motion lives in the soft (low-energy)
modes = mechanically able to open/close like the real motor."""
import sys, numpy as np, scipy.linalg as sla


def load(path, copies=None):
    """Return CA coords X and a subunit label per CA. copies=[(lo,hi),...] within chain A, else use chain id."""
    rows = []
    if path.endswith(".cif"):
        order = []; cols = {}; inl = False
        for L in open(path):
            s = L.strip()
            if s.startswith("_atom_site."):
                order.append(s.split(".", 1)[1]); cols = {k: i for i, k in enumerate(order)}; inl = True; continue
            if inl and s.startswith("ATOM"):
                f = s.split()
                if len(f) >= len(order) and f[cols["label_atom_id"]].strip('"') == "CA":
                    rows.append((f[cols.get("auth_asym_id", cols.get("label_asym_id"))],
                                 int(f[cols.get("auth_seq_id", cols.get("label_seq_id"))]),
                                 [float(f[cols["Cartn_x"]]), float(f[cols["Cartn_y"]]), float(f[cols["Cartn_z"]])]))
            elif inl and s == "#": inl = False
    else:
        for L in open(path):
            if L[:4] == "ATOM" and L[12:16].strip() == "CA":
                rows.append((L[21], int(L[22:26]), [float(L[30:38]), float(L[38:46]), float(L[46:54])]))
    X = np.array([r[2] for r in rows]); lab = []
    for ch, rn, _ in rows:
        if copies:
            s = None
            for k, (lo, hi) in enumerate(copies):
                if lo <= rn <= hi: s = k; break
            lab.append(s)
        else:
            lab.append(ch)
    return X, np.array(lab, dtype=object)


def enm_modes(X, cut=13.0, nmodes=30):
    N = len(X); H = np.zeros((3 * N, 3 * N))
    for i in range(N):
        d = X - X[i]; r2 = (d * d).sum(1); r2[i] = 1e9
        for j in np.where(r2 < cut * cut)[0]:
            e = (X[j] - X[i]); e /= np.linalg.norm(e); k = np.outer(e, e)
            H[3 * i:3 * i + 3, 3 * i:3 * i + 3] += k; H[3 * i:3 * i + 3, 3 * j:3 * j + 3] -= k
    w, V = sla.eigh(H, subset_by_index=[0, nmodes + 5])
    return w[6:], V[:, 6:]  # drop 6 rigid-body


def breathing_vector(X, lab):
    """Radial (in-plane) unit displacement of every CA = ring opening/breathing coordinate."""
    subs = sorted(set(lab), key=lambda x: (x is None, x))
    cent = np.array([X[lab == s].mean(0) for s in subs]); ctr = cent.mean(0)
    _, _, Vt = np.linalg.svd(cent - ctr); normal = Vt[2]
    B = np.zeros(3 * len(X))
    for i, p in enumerate(X):
        v = (p - ctr) - ((p - ctr) @ normal) * normal
        n = np.linalg.norm(v)
        if n > 1e-6:
            B[3 * i:3 * i + 3] = v / n
    return B / np.linalg.norm(B)


def cum_overlap(B, V, ks=(1, 3, 5, 10, 20)):
    ov = (V.T @ B)  # projection of breathing onto each mode
    cum = np.cumsum(ov ** 2)
    return {k: float(cum[min(k, len(cum)) - 1]) for k in ks}


def analyze(path, copies, name):
    X, lab = load(path, copies)
    w, V = enm_modes(X)
    B = breathing_vector(X, lab)
    co = cum_overlap(B, V)
    print(f"{name:<28} N={len(X):<5} breathing captured by softest modes: "
          + " ".join(f"top{k}={co[k]:.2f}" for k in (1, 3, 5, 10, 20)))


if __name__ == "__main__":
    REPO = "/Users/longfu/Developer/claude-science-hackthon/.claude/worktrees/export-gp16-results/gp16_design"
    # native ring (5 chains, non-covalent)
    analyze(f"{REPO}/outputs/structures/cycle3/native_ring__boltz__s1__native_ring_model_0.pdb", None, "NATIVE ring (non-covalent)")
    # cp233 covalent single chain: copy ranges
    cp = [(1, 342), (353, 694), (705, 1046), (1057, 1398), (1409, 1750)]
    analyze(f"{REPO}/outputs/structures/af3_sweep/cp233/inter10/fold_2026_07_08_cp233_int15_inter10_model_0.cif", cp, "cp233 (covalent single chain)")
