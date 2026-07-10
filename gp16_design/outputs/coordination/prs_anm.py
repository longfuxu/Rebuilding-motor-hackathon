#!/usr/bin/env python3
"""Task 1 -- Coordination via elastic-network perturbation-response scanning (PRS)
and GNM cross-correlation, cp233 design C vs native apo ring A.

Both models are built on the SAME residue set (native 4..330 in all 5 subunits,
1635 CA) with identical ANM/GNM parameters, so the numbers are directly comparable.

Coordination proxy = how strongly a perturbation at one subunit's ATP site
(Walker-A P-loop) propagates to its ring-NEIGHBOUR subunit's ATP site, averaged
around the whole ring, normalised by the intra-subunit self-response.

Outputs JSON + prints a compact summary.
"""
import sys, os, json, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import coord_common as cc

from prody import ANM, GNM, calcPerturbResponse, calcCrossCorr, confProDy
confProDy(verbosity='none')

TOP = {
    'design': 'gp16_design/md/openmm_validation/trajectories/C/C_start.pdb',
    'native': 'gp16_design/md/openmm_validation/trajectories/A/A_start.pdb',
}
ANM_CUTOFF = 15.0
GNM_CUTOFF = 10.0


def functional_neighbor(top, rows, xyz):
    """For each subunit k, the neighbour subunit m whose Walker-A centroid is
    closest to k's R146 CA (arginine-finger donor->acceptor direction)."""
    lm = cc.landmark_rows(rows)
    ai = cc.atom_indices(rows)
    wcent = {k: xyz[ai[lm['walkerA'][k]]].mean(0) for k in range(5)}
    r146xyz = {k: xyz[ai[lm['r146'][k]]] for k in range(5)}
    nb = {}
    for k in range(5):
        best, bm = 1e9, None
        for m in range(5):
            if m == k:
                continue
            d = np.linalg.norm(r146xyz[k] - wcent[m])
            if d < best:
                best, bm = d, m
        nb[k] = (bm, best)
    return nb


def analyse(kind):
    t0 = time.time()
    top = cc.load_ca(TOP[kind])
    rows = cc.ca_table(top, kind)
    ai = cc.atom_indices(rows)
    xyz = top.xyz[0][ai] * 10.0   # nm -> Angstrom, in reduced CA order
    N = len(rows)
    lm = cc.landmark_rows(rows)
    order, nxt, prv, cents = cc.ring_neighbor_order(top, rows, xyz=top.xyz[0])
    fnb = functional_neighbor(top, rows, top.xyz[0] * 10.0)

    # ---- ANM + PRS ----
    anm = ANM(kind)
    anm.buildHessian(xyz, cutoff=ANM_CUTOFF, gamma=1.0)
    anm.calcModes(n_modes='all')     # full covariance for PRS
    prs = calcPerturbResponse(anm)   # N x N ; row=source residue, col=response
    if isinstance(prs, tuple):
        prs = prs[0]
    prs = np.asarray(prs)

    def block(mat, src_rows, dst_rows):
        return float(mat[np.ix_(src_rows, dst_rows)].mean())

    wA = lm['walkerA']
    # intra self-response (source ATP site -> same ATP site, incl diagonal)
    intra = np.mean([block(prs, wA[k], wA[k]) for k in range(5)])

    # spatial nearest-neighbour ATP<->ATP coupling around the ring (symmetric avg)
    nn_pairs = [(k, nxt[k]) for k in range(5)]
    nn_couple = np.mean([0.5 * (block(prs, wA[a], wA[b]) + block(prs, wA[b], wA[a]))
                         for a, b in nn_pairs])

    # functional (arginine-finger direction) ATP->ATP coupling k -> neighbour it feeds
    fn_couple = np.mean([block(prs, wA[k], wA[fnb[k][0]]) for k in range(5)])

    # decay around the ring: perturb subunit `order[0]`s ATP site, response at ring dist 1,2
    # average over all subunits as source, group destinations by ring topological distance
    ring_pos = {order[i]: i for i in range(5)}
    by_dist = {1: [], 2: []}
    for k in range(5):
        for m in range(5):
            if m == k:
                continue
            d = min((ring_pos[m] - ring_pos[k]) % 5, (ring_pos[k] - ring_pos[m]) % 5)
            by_dist.setdefault(d, []).append(block(prs, wA[k], wA[m]))
    decay = {d: float(np.mean(v)) for d, v in by_dist.items()}

    # global inter- vs intra-subunit PRS (all residues) -> ring-wide coupling
    sub = cc.subunit_ids(rows)
    inter_mask = sub[:, None] != sub[None, :]
    intra_mask = ~inter_mask
    np.fill_diagonal(intra_mask, False)
    prs_inter_all = float(prs[inter_mask].mean())
    prs_intra_all = float(prs[intra_mask].mean())

    # ---- GNM cross-correlation ----
    gnm = GNM(kind)
    gnm.buildKirchhoff(xyz, cutoff=GNM_CUTOFF, gamma=1.0)
    gnm.calcModes(n_modes='all')
    ccm = calcCrossCorr(gnm)         # normalized cross-correlation, N x N
    ccm = np.asarray(ccm)
    # nearest-neighbour ATP-site cross-correlation
    gnm_nn = np.mean([block(ccm, wA[a], wA[b]) for a, b in nn_pairs])
    gnm_inter_all = float(ccm[inter_mask].mean())
    gnm_intra_all = float(np.abs(ccm[intra_mask]).mean())  # magnitude, intra has self-corr ~1

    res = dict(
        kind=kind, N=N, anm_cutoff=ANM_CUTOFF, gnm_cutoff=GNM_CUTOFF,
        ring_order=[int(x) for x in order],
        functional_neighbor={int(k): [int(fnb[k][0]), round(float(fnb[k][1]), 2)] for k in fnb},
        prs_ATP_intra_self=round(float(intra), 5),
        prs_ATP_NN_coupling=round(float(nn_couple), 5),
        prs_ATP_NN_over_intra=round(float(nn_couple / intra), 4),
        prs_ATP_funcfinger_coupling=round(float(fn_couple), 5),
        prs_decay_by_ringdist={int(k): round(v, 6) for k, v in decay.items()},
        prs_inter_subunit_all=round(prs_inter_all, 6),
        prs_intra_subunit_all=round(prs_intra_all, 6),
        prs_inter_over_intra_all=round(prs_inter_all / prs_intra_all, 4),
        gnm_ATP_NN_crosscorr=round(float(gnm_nn), 5),
        gnm_inter_subunit_all=round(gnm_inter_all, 5),
        seconds=round(time.time() - t0, 1),
    )
    return res


def main():
    out = {}
    for kind in ('design', 'native'):
        print(f"[{kind}] building ANM/GNM + PRS ...", flush=True)
        out[kind] = analyse(kind)
        print(f"[{kind}] done in {out[kind]['seconds']}s", flush=True)
    # ratios design/native
    d, n = out['design'], out['native']
    out['comparison'] = {
        'prs_ATP_NN_coupling_design_over_native':
            round(d['prs_ATP_NN_coupling'] / n['prs_ATP_NN_coupling'], 3),
        'prs_ATP_NN_over_intra_design_vs_native':
            [d['prs_ATP_NN_over_intra'], n['prs_ATP_NN_over_intra']],
        'prs_inter_over_intra_all_design_vs_native':
            [d['prs_inter_over_intra_all'], n['prs_inter_over_intra_all']],
        'gnm_ATP_NN_crosscorr_design_vs_native':
            [d['gnm_ATP_NN_crosscorr'], n['gnm_ATP_NN_crosscorr']],
    }
    path = 'gp16_design/outputs/coordination/prs_anm_results.json'
    with open(path, 'w') as fh:
        json.dump(out, fh, indent=2)
    print("\n==== PRS / ENM coordination summary ====")
    print(json.dumps(out['comparison'], indent=2))
    print(f"\nwrote {path}")


if __name__ == '__main__':
    main()
