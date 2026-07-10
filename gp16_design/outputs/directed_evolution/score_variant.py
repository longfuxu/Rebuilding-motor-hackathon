#!/usr/bin/env python3
"""score_variant.py -- score ONE folded cp233 variant ring for the directed-evolution loop.

WHY A HEAVY-ATOM-CONTACT-WEIGHTED NETWORK
    All variants fold to nearly the SAME closed ring (that is the point of a robust
    scaffold), so a CA-only elastic network cannot tell point mutants apart -- the
    backbone barely moves. The mechanism of a rigidifying substitution is a change in
    SIDECHAIN packing (a Phe/Trp/Pro adds heavy-atom contacts and stiffens local springs).
    We therefore weight every network spring by the number of heavy-atom contacts between
    the two residues in the folded structure:
        gamma_ij = 1 + BETA * (#heavy-atom pairs of residues i,j within CONTACT A)
    A CA-cutoff backbone keeps the network connected (gamma>=1); sidechain contacts add
    stiffness. Mutating Leu->Phe adds atoms -> more contacts on that residue's springs ->
    the coordination / force / soft-mode proxies register the change deterministically,
    on top of whatever backbone change the fold produced. gamma=1 (uniform) reproduces the
    existing CA-ANM baselines (cross-check).

    The proxies here RANK candidates; they do NOT give absolute force/power (only a
    single-molecule experiment can). Judge by M1/M2 + these geometry/coupling proxies,
    never by pLDDT.

PROXIES (all on the same weighted network, node = native residue 4..330 x 5 subunits, N=1635)
    P1 coordination  : ANM perturbation-response ATP-site -> ring-neighbour ATP-site coupling
    P2 force         : static GNM-cross-correlation network, Y129 node betweenness + ATP->Y129->DNA path
    P3 DNA grip      : own-frame radial inward-lining of the 20 DNA-contact residues + pore radius
    P4 interface     : inter-subunit heavy-atom contact density at the designed k->k+1 coupler
    P5 soft-mode     : cumulative overlap of the ring open/close (breathing) coordinate with the
                       softest weighted-ANM modes  (the "still able to power-stroke" metric)

GATES (read/derived; enforced in rank_variants.py)
    G1 M1  : closed compact sequential ring   (from the fold result.json)
    G2 M2  : >=4/5 handedness-robust arginine fingers engaged (from the fold result.json)
    G3 M5  : ATP-pocket CA geometry preserved (pocket RMSD vs WT fold, N-domain superposed)
    G4 soft-mode floor : breathing overlap retained (computed vs WT in rank_variants.py)
"""
import sys, os, json, argparse, warnings
warnings.filterwarnings("ignore")
import numpy as np
import mdtraj as md
import networkx as nx
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
COORD = os.path.abspath(os.path.join(HERE, "..", "coordination"))
sys.path.insert(0, COORD)
import coord_common as cc

from prody import ANM, GNM, Gamma, calcPerturbResponse, calcCrossCorr, confProDy
confProDy(verbosity="none")

FOLDS = os.path.join(HERE, "folds")
WT_CIF = os.path.join(FOLDS, "cp233_WT.cif")

# ---- network / weighting parameters ----
BETA = 0.25            # stiffness added per heavy-atom contact
CONTACT_A = 4.5        # heavy-atom contact distance (A)
ANM_CUTOFF = 15.0
GNM_CUTOFF = 10.0
SOFT_CUTOFF = 13.0
SOFT_NMODES = 25
FORCE_CONTACT = 9.0    # CA-CA contact for the force-transmission graph
EPS = 1e-6

DNA_CONTACT = cc.DNA_CONTACT           # 20 native DNA-contact residues
# M5 catalytic pocket -- Walker-A + Walker-B + arginine finger, PLUS the wider catalytic set
# the phi29 biochem mutation data flags as essential: glu-switch E58, trans-catalysis K105/S106,
# and N158 (a variant that distorts any of these should fail the M5 geometry gate).
POCKET_RES = list(range(24, 32)) + [118, 119, 146, 58, 105, 106, 158]
NDOM = (20, 205)


class MatrixGamma(Gamma):
    """Per-edge force constant looked up from a precomputed NxN gamma matrix (row order =
    coords order). Plugs into prody buildHessian/buildKirchhoff."""
    def __init__(self, G, cutoff):
        self._G = G
        self.cutoff = float(cutoff)
        self.gamma_ = 1.0
    def gamma(self, dist2, i, j):
        return float(self._G[i, j])


def _heavy_by_resseq(top):
    """design resSeq -> np.array of heavy-atom indices."""
    d = {}
    for a in top.atoms:
        if a.residue.is_protein and a.element.symbol != "H":
            d.setdefault(a.residue.resSeq, []).append(a.index)
    return {k: np.array(v) for k, v in d.items()}


def build_gamma(traj, rows):
    """Heavy-atom-contact matrix HC[i,j] (network-row indices) and gamma = 1 + BETA*HC."""
    top = traj.topology
    heavy_by_res = _heavy_by_resseq(top)
    xyz = traj.xyz[0] * 10.0                      # Angstrom
    N = len(rows)
    # collect heavy atoms per network row, tagged by row index
    atom_idx = []
    atom_row = []
    for (r, k, n, ai) in rows:
        drs = cc.design_resseq(k, n)
        h = heavy_by_res.get(drs)
        if h is None:
            continue
        atom_idx.append(h)
        atom_row.append(np.full(len(h), r))
    atom_idx = np.concatenate(atom_idx)
    atom_row = np.concatenate(atom_row)
    coords = xyz[atom_idx]
    tree = cKDTree(coords)
    pairs = tree.query_pairs(CONTACT_A, output_type="ndarray")
    HC = np.zeros((N, N), dtype=np.float32)
    if len(pairs):
        ri = atom_row[pairs[:, 0]]
        rj = atom_row[pairs[:, 1]]
        m = ri != rj                              # ignore same-residue atom pairs
        ri, rj = ri[m], rj[m]
        np.add.at(HC, (ri, rj), 1.0)
        np.add.at(HC, (rj, ri), 1.0)
    G = 1.0 + BETA * HC
    return G, HC


def setup(cif):
    traj = md.load(cif)
    rows = cc.ca_table_design(traj)
    ai = cc.atom_indices(rows)
    xyz = traj.xyz[0][ai] * 10.0                  # CA coords (reduced order), Angstrom
    lm = cc.landmark_rows(rows)
    order, nxt, prv, cents = cc.ring_neighbor_order(traj, rows, xyz=traj.xyz[0])
    G, HC = build_gamma(traj, rows)
    return dict(traj=traj, rows=rows, xyz=xyz, lm=lm, order=order, nxt=nxt, prv=prv,
                G=G, HC=HC, N=len(rows), ai=ai)


# --------------------------------------------------------------- P1 coordination (PRS)
def coordination(S):
    lm, nxt = S["lm"], S["nxt"]
    anm = ANM("v")
    anm.buildHessian(S["xyz"], cutoff=ANM_CUTOFF, gamma=MatrixGamma(S["G"], ANM_CUTOFF))
    anm.calcModes(n_modes="all")
    prs = calcPerturbResponse(anm)
    prs = np.asarray(prs[0] if isinstance(prs, tuple) else prs)

    def block(mat, a, b):
        return float(mat[np.ix_(a, b)].mean())

    wA = lm["walkerA"]
    intra = np.mean([block(prs, wA[k], wA[k]) for k in range(5)])
    nn = np.mean([0.5 * (block(prs, wA[k], wA[nxt[k]]) + block(prs, wA[nxt[k]], wA[k]))
                  for k in range(5)])
    sub = cc.subunit_ids(S["rows"])
    inter = sub[:, None] != sub[None, :]
    intra_mask = ~inter.copy()
    np.fill_diagonal(intra_mask, False)
    return dict(prs_ATP_NN_coupling=round(float(nn), 6),
                prs_ATP_intra_self=round(float(intra), 6),
                prs_ATP_NN_over_intra=round(float(nn / intra), 5),
                prs_inter_over_intra_all=round(float(prs[inter].mean() / prs[intra_mask].mean()), 5))


# --------------------------------------------------------------- P2 force network (static GNM)
def force_network(S):
    lm = S["lm"]
    rows = S["rows"]
    natres = cc.native_res_ids(rows)
    gnm = GNM("v")
    gnm.buildKirchhoff(S["xyz"], cutoff=GNM_CUTOFF, gamma=MatrixGamma(S["G"], GNM_CUTOFF))
    gnm.calcModes(n_modes="all")
    C = np.asarray(calcCrossCorr(gnm))
    # contact graph, weight = -log|cross-corr| (strong co-motion = short edge)
    D = np.linalg.norm(S["xyz"][:, None, :] - S["xyz"][None, :, :], axis=-1)
    contact = (D < FORCE_CONTACT) & (D > 0)
    W = -np.log(np.clip(np.abs(C), EPS, 1.0))
    Gph = nx.Graph()
    Gph.add_nodes_from(range(S["N"]))
    ii, jj = np.where(np.triu(contact, 1))
    for i, j in zip(ii.tolist(), jj.tolist()):
        Gph.add_edge(i, j, weight=float(W[i, j]))
    # ATP -> Y129 -> nearest DNA contact, per subunit
    per = []
    for k in range(5):
        wA = lm["walkerA"][k]
        y = lm["y129"][k]
        dna = [r for r in lm["dna"][k] if r != y]
        try:
            dl_ay, _ = nx.multi_source_dijkstra(Gph, set(wA), target=y, weight="weight")
            dl_yd, _ = nx.multi_source_dijkstra(Gph, {y}, weight="weight")
        except nx.NetworkXNoPath:
            continue
        yd = min((dl_yd[t] for t in dna if t in dl_yd), default=np.nan)
        per.append(dl_ay + yd)
    y_btw_nodes = [lm["y129"][k] for k in range(5)]
    btw = nx.betweenness_centrality(Gph, weight="weight", normalized=True,
                                    k=min(400, Gph.number_of_nodes()), seed=1)
    y_btw = float(np.mean([btw[n] for n in y_btw_nodes]))
    return dict(Y129_betweenness=round(y_btw, 6),
                len_ATP_Y129_DNA=round(float(np.nanmean(per)), 4) if per else None,
                n_edges=Gph.number_of_edges())


# --------------------------------------------------------------- P3 DNA grip (own-frame)
def ring_axis(xyz, sub, natres):
    lo, hi = NDOM
    cents = []
    for k in range(5):
        m = (sub == k) & (natres >= lo) & (natres <= hi)
        cents.append(xyz[m].mean(0))
    cents = np.array(cents)
    center = cents.mean(0)
    _, _, vt = np.linalg.svd(cents - center)
    return center, vt[2]


def grip(S):
    traj, rows = S["traj"], S["rows"]
    top = traj.topology
    heavy_by_res = _heavy_by_resseq(top)
    xyzA = traj.xyz[0] * 10.0
    sub = cc.subunit_ids(rows)
    natres = cc.native_res_ids(rows)
    center, axis = ring_axis(S["xyz"], sub, natres)

    def radial(P):
        d = P - center
        t = d @ axis
        return np.linalg.norm(d - np.outer(t, axis), axis=1)

    inward_N = inward_C = tot_N = tot_C = 0
    radials = []
    for (r, k, n, ai) in rows:
        if n not in DNA_CONTACT:
            continue
        drs = cc.design_resseq(k, n)
        h = heavy_by_res.get(drs)
        if h is None:
            continue
        names = [top.atom(x).name for x in h]
        sc = np.array([xyzA[x] for x, nm in zip(h, names) if nm not in ("N", "C", "O", "CA")])
        if len(sc) == 0:
            sc = xyzA[h]
        r_tip = float(radial(sc).min())
        r_ca = float(radial(xyzA[ai][None, :])[0])
        inward = r_tip < r_ca
        radials.append(r_tip)
        if n <= 205:
            tot_N += 1; inward_N += int(inward)
        else:
            tot_C += 1; inward_C += int(inward)
    # min pore radius (all-CA, along axis)
    d = S["xyz"] - center
    t = d @ axis
    perp = np.linalg.norm(d - np.outer(t, axis), axis=1)
    tlo, thi = np.percentile(t, 8), np.percentile(t, 92)
    prof = []
    for tt in np.arange(tlo, thi, 2.0):
        sel = np.abs(t - tt) < 3.0
        if sel.sum() >= 3:
            prof.append(perp[sel].min())
    min_pore = float(np.min(prof)) if prof else None
    return dict(DNA_inward_N=f"{inward_N}/{tot_N}", DNA_inward_C=f"{inward_C}/{tot_C}",
                DNA_inward_total=inward_N + inward_C,
                DNA_med_radial_tip=round(float(np.median(radials)), 2) if radials else None,
                min_pore_radius=round(min_pore, 2) if min_pore else None,
                admits_dsDNA=bool(min_pore is not None and min_pore >= 9.0))


# --------------------------------------------------------------- P4 interface density
def interface_density(S):
    traj, rows, nxt = S["traj"], S["rows"], S["nxt"]
    top = traj.topology
    heavy_by_res = _heavy_by_resseq(top)
    xyzA = traj.xyz[0] * 10.0
    # heavy atoms per subunit (network residues only)
    sub_atoms = {k: [] for k in range(5)}
    for (r, k, n, ai) in rows:
        h = heavy_by_res.get(cc.design_resseq(k, n))
        if h is not None:
            sub_atoms[k].append(h)
    sub_atoms = {k: np.concatenate(v) for k, v in sub_atoms.items()}
    counts = []
    for k in range(5):
        m = nxt[k]
        A = xyzA[sub_atoms[k]]
        B = xyzA[sub_atoms[m]]
        tree = cKDTree(B)
        c = sum(len(tree.query_ball_point(a, CONTACT_A)) for a in A)
        counts.append(c)
    return dict(interSub_heavy_contacts=round(float(np.mean(counts)), 1),
                interSub_heavy_contacts_min=int(np.min(counts)))


# --------------------------------------------------------------- P5 soft-mode (breathing)
def breathing_vector(xyz, sub):
    cents = np.array([xyz[sub == k].mean(0) for k in range(5)])
    ctr = cents.mean(0)
    _, _, vt = np.linalg.svd(cents - ctr)
    normal = vt[2]
    B = np.zeros(3 * len(xyz))
    for i, p in enumerate(xyz):
        v = (p - ctr) - ((p - ctr) @ normal) * normal
        nrm = np.linalg.norm(v)
        if nrm > 1e-6:
            B[3 * i:3 * i + 3] = v / nrm
    return B / np.linalg.norm(B)


def soft_mode(S):
    sub = cc.subunit_ids(S["rows"])
    anm = ANM("s")
    anm.buildHessian(S["xyz"], cutoff=SOFT_CUTOFF, gamma=MatrixGamma(S["G"], SOFT_CUTOFF))
    anm.calcModes(n_modes=SOFT_NMODES)
    V = anm.getEigvecs()                          # (3N, nmodes), softest first (trivial removed)
    B = breathing_vector(S["xyz"], sub)
    ov = V.T @ B
    cum = np.cumsum(ov ** 2)
    ks = (1, 3, 5, 10, 20)
    return {f"breathing_top{k}": round(float(cum[min(k, len(cum)) - 1]), 4) for k in ks}


# --------------------------------------------------------------- G3 M5 pocket geometry
def _kabsch(P, Q):
    """rotate/translate Q onto P; return rmsd."""
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    H = Qc.T @ Pc
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    Qr = Qc @ R.T
    return float(np.sqrt(((Qr - Pc) ** 2).sum(1).mean()))


def pocket_rmsd_vs_wt(S, wt_S):
    """Per subunit: superpose N-domain (20-205) CA onto WT subunit-0 N-domain, pocket CA RMSD."""
    def ca_map(St, k):
        rows = St["rows"]; xyz = St["xyz"]
        return {n: xyz[r] for (r, kk, n, ai) in rows if kk == k}
    wt0 = ca_map(wt_S, 0)
    lo, hi = NDOM
    nd_res = [r for r in range(lo, hi + 1) if r in wt0]
    pk_res = [r for r in POCKET_RES if r in wt0]
    rmsds = []
    for k in range(5):
        vk = ca_map(S, k)
        nd_common = [r for r in nd_res if r in vk]
        pk_common = [r for r in pk_res if r in vk]
        if len(nd_common) < 20 or len(pk_common) < 6:
            continue
        P = np.array([wt0[r] for r in nd_common]); Q = np.array([vk[r] for r in nd_common])
        # build transform from N-domain, apply to pocket
        Pc0 = P.mean(0); Qc0 = Q.mean(0)
        H = (Q - Qc0).T @ (P - Pc0)
        U, _, Vt = np.linalg.svd(H)
        dsign = np.sign(np.linalg.det(Vt.T @ U.T))
        R = Vt.T @ np.diag([1, 1, dsign]) @ U.T
        Pp = np.array([wt0[r] for r in pk_common])
        Qp = np.array([vk[r] for r in pk_common])
        Qp_fit = (Qp - Qc0) @ R.T + Pc0
        rmsds.append(float(np.sqrt(((Qp_fit - Pp) ** 2).sum(1).mean())))
    return round(float(np.mean(rmsds)), 3) if rmsds else None


# --------------------------------------------------------------- fold gate (M1/M2)
def read_gate(name):
    rj = os.path.join(FOLDS, f"{name}.result.json")
    if not os.path.exists(rj):
        return dict(M1="", M2="", m1_pass=None, m2_engaged=None, m2_pass=None, fold_ok=None)
    d = json.load(open(rj))
    m1 = d.get("M1", "")
    ring = d.get("ring", {}) or {}
    eng = ring.get("engaged")
    m1_pass = ("compact_ring True" in m1) and ("sequential_consistent: YES" in m1)
    return dict(M1=m1, M2=d.get("M2", ""), m1_pass=bool(m1_pass),
                m2_engaged=eng, m2_pass=(eng is not None and eng >= 4),
                fold_ok=bool(d.get("ok")))


# --------------------------------------------------------------- driver
def score_cif(name, cif, wt_S=None):
    S = setup(cif)
    rec = dict(name=name, cif=os.path.abspath(cif), N=S["N"], BETA=BETA, CONTACT_A=CONTACT_A)
    rec.update(coordination(S))
    rec.update(force_network(S))
    rec.update(grip(S))
    rec.update(interface_density(S))
    rec.update(soft_mode(S))
    rec.update(read_gate(name))
    if wt_S is not None:
        rec["M5_pocket_rmsd_vs_wt"] = pocket_rmsd_vs_wt(S, wt_S)
    else:
        rec["M5_pocket_rmsd_vs_wt"] = 0.0
    return rec, S


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cif")
    ap.add_argument("--name", default=None)
    ap.add_argument("--wt", default=WT_CIF)
    args = ap.parse_args()
    name = args.name or os.path.basename(args.cif).replace(".cif", "")
    wt_S = setup(args.wt) if (args.wt and os.path.abspath(args.wt) != os.path.abspath(args.cif)) else None
    rec, _ = score_cif(name, args.cif, wt_S=wt_S)
    print(json.dumps(rec, indent=2))


if __name__ == "__main__":
    main()
