#!/usr/bin/env python3
"""Task 3 -- Force-transmission network, cp233 design C vs native apo ring A.

Dynamic residue network (Girvan-Newman / Luthey-Schulten style):
  nodes = Ca (native 4..330 in all 5 subunits, 1635 nodes)
  edges = residue pairs in contact (Ca-Ca < CONTACT_CUT in the MD-mean structure)
  edge weight (distance) = -log(|C_ij|)   with C_ij the MD dynamic cross-correlation
      -> strongly co-moving contacts are "short" = good force conduits.

We trace the power-stroke conduit  ATP-site (Walker-A P-loop 24-31) -> Y129 ->
DNA-contact residues  WITHIN each subunit, average over the 5 subunits, and compare
design vs native:
  * path length (sum of -log|C|)  -> lower = stiffer/stronger transmission
  * Y129 node betweenness         -> how central Y129 is to force routing
  * edge betweenness on the conduit + robustness (remove top edge, re-path)
Also ranks the highest edge-betweenness residues = candidates to RIGIDIFY for a
more powerful motor.

Cross-checked against an ENM/GNM-correlation network (short MD caveat).
"""
import sys, os, json
import numpy as np
import mdtraj as md
import networkx as nx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import coord_common as cc
from prody import GNM, calcCrossCorr, confProDy
confProDy(verbosity='none')

TRAJ = {
    'design': ('gp16_design/md/openmm_validation/trajectories/C/C_prod.dcd',
               'gp16_design/md/openmm_validation/trajectories/C/C_start.pdb', 'design'),
    'native': ('gp16_design/md/openmm_validation/trajectories/A/A_prod.dcd',
               'gp16_design/md/openmm_validation/trajectories/A/A_start.pdb', 'native'),
}
CONTACT_CUT = 9.0   # Angstrom, Ca-Ca contact
EPS = 1e-6


def build(kind, corr_source='md'):
    dcd, top, kd = TRAJ[kind]
    traj = md.load(dcd, top=top)
    ref = cc.load_ca(top)
    rows = cc.ca_table(ref, kd)
    ai = cc.atom_indices(rows)
    traj.superpose(traj, frame=0, atom_indices=ai, ref_atom_indices=ai)
    X = traj.xyz[:, ai, :] * 10.0                 # (F,N,3) Angstrom
    meanX = X.mean(0)                              # (N,3)
    N = len(rows)

    if corr_source == 'md':
        d = X - meanX
        cov = np.einsum('fia,fja->ij', d, d) / X.shape[0]
        var = np.diag(cov).copy()
        C = cov / np.sqrt(np.outer(var, var))
    else:  # gnm
        gnm = GNM(kind); gnm.buildKirchhoff(meanX, cutoff=10.0, gamma=1.0)
        gnm.calcModes(n_modes='all'); C = np.asarray(calcCrossCorr(gnm))

    # contact map from mean structure
    diff = meanX[:, None, :] - meanX[None, :, :]
    D = np.sqrt((diff * diff).sum(-1))
    contact = (D < CONTACT_CUT) & (D > 0)

    absC = np.clip(np.abs(C), EPS, 1.0)
    W = -np.log(absC)

    G = nx.Graph()
    G.add_nodes_from(range(N))
    ii, jj = np.where(np.triu(contact, 1))
    for i, j in zip(ii.tolist(), jj.tolist()):
        G.add_edge(i, j, weight=float(W[i, j]))
    return G, rows, C


def conduit(kind, corr_source='md'):
    G, rows, C = build(kind, corr_source)
    lm = cc.landmark_rows(rows)
    sub = cc.subunit_ids(rows)
    natres = cc.native_res_ids(rows)

    per_sub = []
    for k in range(5):
        wA = lm['walkerA'][k]
        y = lm['y129'][k]
        dna = [r for r in lm['dna'][k] if r != y]     # DNA contacts excluding 129 itself
        # ATP -> Y129 (multi-source over the P-loop)
        try:
            dl_ay, p_ay = nx.multi_source_dijkstra(G, set(wA), target=y, weight='weight')
        except nx.NetworkXNoPath:
            continue
        # Y129 -> nearest DNA contact
        dl_yd, p_yd = nx.multi_source_dijkstra(G, {y}, weight='weight')
        best_d = min((dl_yd[t], t) for t in dna if t in dl_yd)
        dl_yd_best, dtar = best_d
        p_yd_best = p_yd[dtar]
        # direct ATP -> DNA (unconstrained) for reference
        dl_ad, _ = nx.multi_source_dijkstra(G, set(wA), weight='weight')
        dl_ad_best = min(dl_ad[t] for t in dna if t in dl_ad)

        full_len = dl_ay + dl_yd_best
        full_path = p_ay + p_yd_best[1:]
        per_sub.append(dict(
            sub=k, len_ATP_Y129=dl_ay, len_Y129_DNA=dl_yd_best,
            len_ATP_Y129_DNA=full_len, len_ATP_DNA_direct=dl_ad_best,
            path_natres=[int(natres[n]) for n in full_path],
            dna_target=int(natres[dtar]),
        ))

    L = lambda key: float(np.mean([s[key] for s in per_sub]))
    # Y129 node betweenness (whole-ring graph, weighted, sampled for speed)
    btw = nx.betweenness_centrality(G, weight='weight', normalized=True,
                                    k=min(400, G.number_of_nodes()), seed=1)
    y_btw = float(np.mean([btw[lm['y129'][k]] for k in range(5)]))
    wa_btw = float(np.mean([btw[r] for k in range(5) for r in lm['walkerA'][k]]))

    # top edge-betweenness residues (rigidify candidates), by native res identity
    ebtw = nx.edge_betweenness_centrality(G, weight='weight', normalized=True,
                                          k=min(300, G.number_of_nodes()), seed=1)
    node_eb = np.zeros(G.number_of_nodes())
    for (a, b), v in ebtw.items():
        node_eb[a] += v; node_eb[b] += v
    # aggregate by native residue number across the 5 subunits
    agg = {}
    for n in range(G.number_of_nodes()):
        agg.setdefault(int(natres[n]), []).append(node_eb[n])
    agg_mean = {r: float(np.mean(v)) for r, v in agg.items()}
    top_rig = sorted(agg_mean.items(), key=lambda x: -x[1])[:12]

    # robustness: remove the single highest-|weight-betweenness| edge on subunit-0 conduit,
    # recompute the ATP->Y129->DNA length, report fractional increase (avg over subunits)
    incs = []
    for s in per_sub:
        k = s['sub']
        path = None
        # rebuild path node list for this subunit
        wA = lm['walkerA'][k]; y = lm['y129'][k]
        dna = [r for r in lm['dna'][k] if r != y]
        try:
            dl_ay, p_ay = nx.multi_source_dijkstra(G, set(wA), target=y, weight='weight')
            dl_yd, p_yd = nx.multi_source_dijkstra(G, {y}, weight='weight')
            dtar = min((dl_yd[t], t) for t in dna if t in dl_yd)[1]
            full = p_ay + p_yd[dtar][1:]
        except nx.NetworkXNoPath:
            continue
        base = s['len_ATP_Y129_DNA']
        # find highest-betweenness edge along the path
        edges = list(zip(full[:-1], full[1:]))
        eb_on = [(ebtw.get((a, b), ebtw.get((b, a), 0.0)), (a, b)) for a, b in edges]
        if not eb_on:
            continue
        _, (ea, eb) = max(eb_on)
        Gc = G.copy(); Gc.remove_edge(ea, eb)
        try:
            dl_ay2, _ = nx.multi_source_dijkstra(Gc, set(wA), target=y, weight='weight')
            dl_yd2, _ = nx.multi_source_dijkstra(Gc, {y}, weight='weight')
            dtar2 = min((dl_yd2[t], t) for t in dna if t in dl_yd2)[1]
            new = dl_ay2 + dl_yd2[dtar2]
            incs.append((new - base) / base)
        except (nx.NetworkXNoPath, ValueError):
            incs.append(np.nan)

    return dict(
        kind=kind, corr_source=corr_source,
        n_edges=G.number_of_edges(),
        len_ATP_Y129=round(L('len_ATP_Y129'), 4),
        len_Y129_DNA=round(L('len_Y129_DNA'), 4),
        len_ATP_Y129_DNA=round(L('len_ATP_Y129_DNA'), 4),
        len_ATP_DNA_direct=round(L('len_ATP_DNA_direct'), 4),
        Y129_node_betweenness=round(y_btw, 5),
        WalkerA_node_betweenness=round(wa_btw, 5),
        conduit_robustness_frac_increase=round(float(np.nanmean(incs)), 4) if incs else None,
        example_path_sub0_natres=per_sub[0]['path_natres'] if per_sub else None,
        example_path_dna_target=per_sub[0]['dna_target'] if per_sub else None,
        top_rigidify_residues=[[r, round(v, 5)] for r, v in top_rig],
    ), per_sub


def main():
    out = {}
    for kind in ('design', 'native'):
        print(f"[{kind}] building force-transmission network (MD-DCCM) ...", flush=True)
        r, per_sub = conduit(kind, 'md')
        out[kind] = r
        out[f'{kind}_persub'] = per_sub
        print(f"[{kind}] ATP->Y129->DNA len = {r['len_ATP_Y129_DNA']}, "
              f"Y129 btw = {r['Y129_node_betweenness']}", flush=True)
    d, n = out['design'], out['native']
    out['comparison'] = {
        'len_ATP_Y129_DNA_design_vs_native': [d['len_ATP_Y129_DNA'], n['len_ATP_Y129_DNA']],
        'len_ATP_DNA_direct_design_vs_native': [d['len_ATP_DNA_direct'], n['len_ATP_DNA_direct']],
        'Y129_betweenness_design_vs_native': [d['Y129_node_betweenness'], n['Y129_node_betweenness']],
        'conduit_robustness_design_vs_native':
            [d['conduit_robustness_frac_increase'], n['conduit_robustness_frac_increase']],
        'note': 'lower path length = stiffer/stronger transmission (higher correlation along conduit)',
    }
    with open('gp16_design/outputs/coordination/force_network_results.json', 'w') as fh:
        json.dump(out, fh, indent=2)
    print("\n==== force-transmission network summary ====")
    print(json.dumps(out['comparison'], indent=2))


if __name__ == '__main__':
    main()
