#!/usr/bin/env python3
"""Task 2 (ClpX) -- Force-transmission network, single-chain ClpX design vs native 6PP5.

ClpX port of gp16_design/outputs/coordination/force_network.py, GNM-correlation variant
(NO MD available for ClpX -> we use the elastic-network cross-correlation as the edge
weight, i.e. the gp16 study's "gnm cross-check" path, not its MD-DCCM path).

Dynamic residue network:
  nodes  = Ca (native residues in the reduced set, 6 subunits)
  edges  = residue pairs in contact (Ca-Ca < CONTACT_CUT in the model)
  weight = -log(|C_ij|), C from GNM cross-correlation -> co-moving contacts are "short"
           = good force conduits.

Two comparisons:
  (A) NODE-MATCHED (fair, identical nodes for both, like gp16): joint common residue set.
      Report betweenness centrality of the trans arginine finger R307 (the ClpX coupling
      residue, gp16 R146 analog), Walker-A, sensor-II R370, and pore-grip residues
      155/157/202. This is the rigorous design-vs-native force-hub comparison.
  (B) SUPPLEMENTARY conduit ATP-site(Walker-A) -> pore-1 aromatic Y153 (substrate grip,
      gp16 Y129 analog) -> pore-contact, traced per subunit on each structure's FULL
      modeled Ca set (design: 6 complete copies; native: the ordered subunits where Y153
      is modeled -- the disengaged ADP seam F has a disordered pore-1 loop). Y153 lives in
      the disordered native seam so it cannot enter the node-matched set (A); (B) reports
      its betweenness with the node-count caveat noted.
"""
import sys, os, json
import numpy as np
import networkx as nx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import clpx_coord_common as cc
from prody import GNM, calcCrossCorr, confProDy
confProDy(verbosity='none')

CONTACT_CUT = 9.0     # Angstrom, Ca-Ca contact
GNM_CUTOFF = 10.0
EPS = 1e-6
NSUB = cc.NSUB


def build_graph(xyz, gnm_kind):
    N = len(xyz)
    gnm = GNM(gnm_kind)
    gnm.buildKirchhoff(xyz, cutoff=GNM_CUTOFF, gamma=1.0)
    gnm.calcModes(n_modes='all')
    C = np.asarray(calcCrossCorr(gnm))
    diff = xyz[:, None, :] - xyz[None, :, :]
    D = np.sqrt((diff * diff).sum(-1))
    contact = (D < CONTACT_CUT) & (D > 0)
    absC = np.clip(np.abs(C), EPS, 1.0)
    W = -np.log(absC)
    G = nx.Graph()
    G.add_nodes_from(range(N))
    ii, jj = np.where(np.triu(contact, 1))
    for i, j in zip(ii.tolist(), jj.tolist()):
        G.add_edge(i, j, weight=float(W[i, j]))
    return G, C


# ---------------- (A) node-matched betweenness on joint common set ----------------
def node_matched(kind, common):
    D = cc.load_common(kind, common_override=common)
    G, C = build_graph(D['xyz'], kind)
    lm = cc.landmark_rows(D)
    natres = D['natres']

    btw = nx.betweenness_centrality(G, weight='weight', normalized=True,
                                    k=min(500, G.number_of_nodes()), seed=1)

    def res_btw(resnum):
        rows = [lm['lut'][(k, resnum)] for k in range(NSUB) if (k, resnum) in lm['lut']]
        return float(np.mean([btw[r] for r in rows])) if rows else None

    arg_btw = res_btw(cc.ARGFINGER)
    s2_btw = res_btw(cc.SENSOR2)
    wa_btw = float(np.mean([btw[r] for k in range(NSUB) for r in lm['walkerA'][k]]))
    pore_common = [r for r in cc.PORE_CONTACT if r in set(D['common'])]
    pore_btw = float(np.mean([btw[lm['lut'][(k, r)]]
                              for k in range(NSUB) for r in pore_common
                              if (k, r) in lm['lut']]))

    # intra-subunit conduit ATP(Walker-A) -> nearest pore-grip residue (155/157/202)
    lens = []
    for k in range(NSUB):
        wA = lm['walkerA'][k]
        pore_rows = [lm['lut'][(k, r)] for r in pore_common if (k, r) in lm['lut']]
        try:
            dl, _ = nx.multi_source_dijkstra(G, set(wA), weight='weight')
            best = min(dl[t] for t in pore_rows if t in dl)
            lens.append(best)
        except (nx.NetworkXNoPath, ValueError):
            pass
    len_atp_pore = float(np.mean(lens)) if lens else None

    # rigidify ranking: top edge-betweenness residues by native identity
    ebtw = nx.edge_betweenness_centrality(G, weight='weight', normalized=True,
                                          k=min(400, G.number_of_nodes()), seed=1)
    node_eb = np.zeros(G.number_of_nodes())
    for (a, b), v in ebtw.items():
        node_eb[a] += v
        node_eb[b] += v
    agg = {}
    for n in range(G.number_of_nodes()):
        agg.setdefault(int(natres[n]), []).append(node_eb[n])
    agg_mean = {r: float(np.mean(v)) for r, v in agg.items()}
    top_rig = sorted(agg_mean.items(), key=lambda x: -x[1])[:12]

    return dict(
        kind=kind, N=G.number_of_nodes(), n_edges=G.number_of_edges(),
        R307_argfinger_betweenness=round(arg_btw, 6) if arg_btw is not None else None,
        R370_sensorII_betweenness=round(s2_btw, 6) if s2_btw is not None else None,
        WalkerA_betweenness=round(wa_btw, 6),
        pore_grip_betweenness=round(pore_btw, 6),
        pore_grip_residues=pore_common,
        len_ATP_to_pore=round(len_atp_pore, 4) if len_atp_pore is not None else None,
        top_rigidify_residues=[[r, round(v, 5)] for r, v in top_rig],
    )


# ---------------- (B) supplementary Y153 conduit on full modeled sets ----------------
def full_maps(kind):
    """Per-subunit CA map on each structure's OWN modeled residues (not intersected)."""
    sub_maps, _ = cc._ca_maps(kind)
    rows, xyz = [], []
    r = 0
    for k in range(NSUB):
        ca, _ = sub_maps[k]
        for nr in sorted(ca.keys()):
            rows.append((r, k, nr))
            xyz.append(ca[nr])
            r += 1
    return np.asarray(xyz, float), rows


def y153_conduit(kind):
    xyz, rows = full_maps(kind)
    G, C = build_graph(xyz, kind + '_full')
    lut = {(k, nr): r for (r, k, nr) in rows}
    btw = nx.betweenness_centrality(G, weight='weight', normalized=True,
                                    k=min(500, G.number_of_nodes()), seed=1)
    y_rows, conduit_lens, subs_used = [], [], []
    for k in range(NSUB):
        y = lut.get((k, cc.YFORCE))
        wA = [lut[(k, n)] for n in cc.WALKERA if (k, n) in lut]
        pore = [lut[(k, n)] for n in cc.PORE_CONTACT if (k, n) in lut and (k, n) != (k, cc.YFORCE)]
        if y is None or not wA or not pore:
            continue
        try:
            dl_ay, _ = nx.multi_source_dijkstra(G, set(wA), target=y, weight='weight')
            dl_yp, _ = nx.multi_source_dijkstra(G, {y}, weight='weight')
            best_yp = min(dl_yp[t] for t in pore if t in dl_yp)
        except (nx.NetworkXNoPath, ValueError):
            continue
        y_rows.append(btw[y])
        conduit_lens.append(dl_ay + best_yp)
        subs_used.append(k)
    return dict(
        kind=kind, N=G.number_of_nodes(),
        Y153_betweenness=round(float(np.mean(y_rows)), 6) if y_rows else None,
        len_ATP_Y153_pore=round(float(np.mean(conduit_lens)), 4) if conduit_lens else None,
        subunits_with_intact_conduit=subs_used,
        n_subunits_used=len(subs_used),
    )


def main():
    dcom = set(cc.load_common('design')['common'])
    ncom = set(cc.load_common('native')['common'])
    common = sorted(dcom & ncom)
    out = {'joint_common_n': len(common), 'contact_cut': CONTACT_CUT,
           'edge_weight': '-log|C_GNM|  (NO MD -> GNM cross-correlation, gp16 no-MD path)'}

    print("== (A) node-matched betweenness (joint common set) ==", flush=True)
    for kind in ('design', 'native'):
        out['A_' + kind] = node_matched(kind, common)
        a = out['A_' + kind]
        print(f"[{kind}] R307 btw={a['R307_argfinger_betweenness']}  "
              f"WalkerA btw={a['WalkerA_betweenness']}  pore btw={a['pore_grip_betweenness']}  "
              f"len ATP->pore={a['len_ATP_to_pore']}", flush=True)

    print("\n== (B) supplementary Y153 conduit (full modeled sets) ==", flush=True)
    for kind in ('design', 'native'):
        out['B_' + kind] = y153_conduit(kind)
        b = out['B_' + kind]
        print(f"[{kind}] Y153 btw={b['Y153_betweenness']}  "
              f"len ATP->Y153->pore={b['len_ATP_Y153_pore']}  "
              f"subunits used={b['n_subunits_used']} ({b['subunits_with_intact_conduit']})", flush=True)

    da, na = out['A_design'], out['A_native']
    db, nb = out['B_design'], out['B_native']
    out['comparison'] = {
        'R307_argfinger_betweenness_design_vs_native':
            [da['R307_argfinger_betweenness'], na['R307_argfinger_betweenness']],
        'R307_design_over_native':
            round(da['R307_argfinger_betweenness'] / na['R307_argfinger_betweenness'], 3)
            if na['R307_argfinger_betweenness'] else None,
        'pore_grip_betweenness_design_vs_native':
            [da['pore_grip_betweenness'], na['pore_grip_betweenness']],
        'Y153_betweenness_design_vs_native':
            [db['Y153_betweenness'], nb['Y153_betweenness']],
        'Y153_design_over_native':
            round(db['Y153_betweenness'] / nb['Y153_betweenness'], 3)
            if nb['Y153_betweenness'] else None,
        'gp16_reference_Y129_betweenness_design_over_native_GNM': 1.65,
        'note': ('lower path length = stiffer transmission; higher betweenness = more '
                 'central force-routing node. (B) native uses only seam-intact subunits '
                 '(fewer nodes) -> Y153 comparison is directional, not node-matched.'),
    }
    path = 'gp16_design/outputs/clpx_coupling/clpx_force_network_results.json'
    with open(path, 'w') as fh:
        json.dump(out, fh, indent=2)
    print("\n==== ClpX force-transmission network summary ====")
    print(json.dumps(out['comparison'], indent=2))
    print(f"\nwrote {path}")


if __name__ == '__main__':
    main()
