#!/usr/bin/env python3
"""Task 1 (ClpX) -- Inter-subunit coordination via elastic-network perturbation-response
scanning (PRS) + GNM cross-correlation, single-chain ClpX design vs native 6PP5 hexamer.

Direct ClpX port of gp16_design/outputs/coordination/prs_anm.py. Both rings are built on
the SAME residue set (native residues present in all 6 native chains AND all 6 design
copies) with identical ANM/GNM parameters, so the numbers are directly comparable.

Coordination proxy = how strongly a force perturbation at one subunit's ATP site
(Walker-A P-loop 119-128) propagates to its ring-NEIGHBOUR subunit's ATP site, averaged
around the whole 6-membered ring, normalised by the intra-subunit self-response.
Higher = more inter-subunit coordinated.

NO MD for ClpX -> ENM/PRS + GNM only (see clpx_coord_common.py header for caveats).
"""
import sys, os, json, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import clpx_coord_common as cc

from prody import ANM, GNM, calcPerturbResponse, calcCrossCorr, confProDy
confProDy(verbosity='none')

ANM_CUTOFF = 15.0
GNM_CUTOFF = 10.0
NSUB = cc.NSUB


def block(mat, src_rows, dst_rows):
    return float(mat[np.ix_(src_rows, dst_rows)].mean())


def analyse(kind, common):
    t0 = time.time()
    D = cc.load_common(kind, common_override=common)
    xyz = D['xyz']
    N = len(D['rows'])
    lm = cc.landmark_rows(D)
    order, nxt, prv, cents = cc.ring_neighbor_order(D)
    fnb = cc.functional_neighbor(D)
    wA = lm['walkerA']

    # ---- ANM + PRS ----
    anm = ANM(kind)
    anm.buildHessian(xyz, cutoff=ANM_CUTOFF, gamma=1.0)
    anm.calcModes(n_modes='all')
    prs = calcPerturbResponse(anm)
    if isinstance(prs, tuple):
        prs = prs[0]
    prs = np.asarray(prs)

    # intra self-response: ATP site -> same ATP site
    intra = np.mean([block(prs, wA[k], wA[k]) for k in range(NSUB)])

    # spatial nearest-neighbour ATP<->ATP coupling around the ring (symmetric avg)
    nn_pairs = [(k, nxt[k]) for k in range(NSUB)]
    nn_couple = np.mean([0.5 * (block(prs, wA[a], wA[b]) + block(prs, wA[b], wA[a]))
                         for a, b in nn_pairs])

    # functional (arg-finger R307 donor->acceptor) ATP->ATP coupling: k feeds fnb[k]
    fn_couple = np.mean([block(prs, wA[k], wA[fnb[k][0]]) for k in range(NSUB)])

    # decay around the ring by topological distance (1,2,3 for hexamer)
    ring_pos = {order[i]: i for i in range(NSUB)}
    by_dist = {}
    for k in range(NSUB):
        for m in range(NSUB):
            if m == k:
                continue
            d = min((ring_pos[m] - ring_pos[k]) % NSUB, (ring_pos[k] - ring_pos[m]) % NSUB)
            by_dist.setdefault(d, []).append(block(prs, wA[k], wA[m]))
    decay = {int(d): float(np.mean(v)) for d, v in by_dist.items()}

    # global inter- vs intra-subunit PRS (all residues)
    sub = D['sub']
    inter_mask = sub[:, None] != sub[None, :]
    intra_mask = ~inter_mask
    np.fill_diagonal(intra_mask, False)
    prs_inter_all = float(prs[inter_mask].mean())
    prs_intra_all = float(prs[intra_mask].mean())

    # ---- GNM cross-correlation ----
    gnm = GNM(kind)
    gnm.buildKirchhoff(xyz, cutoff=GNM_CUTOFF, gamma=1.0)
    gnm.calcModes(n_modes='all')
    ccm = np.asarray(calcCrossCorr(gnm))
    gnm_nn = np.mean([block(ccm, wA[a], wA[b]) for a, b in nn_pairs])
    gnm_inter_all = float(ccm[inter_mask].mean())

    res = dict(
        kind=kind, N=N, n_residues=len(D['common']),
        anm_cutoff=ANM_CUTOFF, gnm_cutoff=GNM_CUTOFF,
        ring_order=[int(x) for x in order],
        functional_neighbor={int(k): [int(fnb[k][0]), fnb[k][1]] for k in fnb},
        prs_ATP_intra_self=round(float(intra), 6),
        prs_ATP_NN_coupling=round(float(nn_couple), 6),
        prs_ATP_NN_over_intra=round(float(nn_couple / intra), 5),
        prs_ATP_argfinger_coupling=round(float(fn_couple), 6),
        prs_decay_by_ringdist={k: round(v, 6) for k, v in decay.items()},
        prs_ringdist_retention_2_over_1=round(decay[2] / decay[1], 4) if 1 in decay and 2 in decay else None,
        prs_inter_subunit_all=round(prs_inter_all, 6),
        prs_intra_subunit_all=round(prs_intra_all, 6),
        prs_inter_over_intra_all=round(prs_inter_all / prs_intra_all, 5),
        gnm_ATP_NN_crosscorr=round(float(gnm_nn), 5),
        gnm_inter_subunit_all=round(gnm_inter_all, 5),
        seconds=round(time.time() - t0, 1),
    )
    return res


def main():
    # joint common residue set = intersection of design & native modeled residues
    dcom = set(cc.load_common('design')['common'])
    ncom = set(cc.load_common('native')['common'])
    common = sorted(dcom & ncom)
    print(f"joint common residue set: n={len(common)} "
          f"({common[0]}..{common[-1]}); design-only {len(dcom-ncom)}, native-only {len(ncom-dcom)}",
          flush=True)

    out = {'joint_common_n': len(common),
           'joint_common_range': [common[0], common[-1]]}
    for kind in ('design', 'native'):
        print(f"[{kind}] building ANM/GNM + PRS on N={NSUB*len(common)} ...", flush=True)
        out[kind] = analyse(kind, common)
        print(f"[{kind}] done in {out[kind]['seconds']}s  "
              f"PRS ATP->NN = {out[kind]['prs_ATP_NN_coupling']}", flush=True)

    d, n = out['design'], out['native']
    out['comparison'] = {
        'prs_ATP_NN_coupling_design_over_native':
            round(d['prs_ATP_NN_coupling'] / n['prs_ATP_NN_coupling'], 3),
        'prs_ATP_argfinger_coupling_design_over_native':
            round(d['prs_ATP_argfinger_coupling'] / n['prs_ATP_argfinger_coupling'], 3),
        'prs_ATP_NN_over_intra_design_vs_native':
            [d['prs_ATP_NN_over_intra'], n['prs_ATP_NN_over_intra']],
        'prs_ringdist2_coupling_design_over_native':
            round(d['prs_decay_by_ringdist'][2] / n['prs_decay_by_ringdist'][2], 3),
        'prs_inter_over_intra_all_design_vs_native':
            [d['prs_inter_over_intra_all'], n['prs_inter_over_intra_all']],
        'gnm_ATP_NN_crosscorr_design_vs_native':
            [d['gnm_ATP_NN_crosscorr'], n['gnm_ATP_NN_crosscorr']],
        'gp16_reference_prs_ATP_NN_design_over_native': 1.788,
    }
    path = 'gp16_design/outputs/clpx_coupling/clpx_prs_anm_results.json'
    with open(path, 'w') as fh:
        json.dump(out, fh, indent=2)
    print("\n==== ClpX PRS / ENM coordination summary ====")
    print(json.dumps(out['comparison'], indent=2))
    print(f"\nwrote {path}")


if __name__ == '__main__':
    main()
