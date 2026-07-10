#!/usr/bin/env python3
"""Task 2 -- MD dynamic cross-correlation (DCCM) of Ca fluctuations, cp233 design
C vs native apo ring A, from the existing OpenMM production trajectories.

Restricted to the SAME residue set (native 4..330 in all 5 subunits, 1635 CA) so
inter-subunit off-diagonal blocks are directly comparable.

CAVEAT (stated in the report too): the trajectories are SHORT (~70-80 frames,
sub-ns of sampling); DCCM off-diagonal terms are noisy. Treat as a qualitative
cross-check on the ENM/PRS result, not a converged equilibrium measurement.
"""
import sys, os, json
import numpy as np
import mdtraj as md

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import coord_common as cc

TRAJ = {
    'design': ('gp16_design/md/openmm_validation/trajectories/C/C_prod.dcd',
               'gp16_design/md/openmm_validation/trajectories/C/C_start.pdb', 'design'),
    'native': ('gp16_design/md/openmm_validation/trajectories/A/A_prod.dcd',
               'gp16_design/md/openmm_validation/trajectories/A/A_start.pdb', 'native'),
}


def dccm(kind):
    dcd, top, kd = TRAJ[kind]
    traj = md.load(dcd, top=top)
    ref = cc.load_ca(top)
    rows = cc.ca_table(ref, kd)
    ai = cc.atom_indices(rows)
    # superpose all frames on the reduced CA set (remove rigid-body motion)
    traj.superpose(traj, frame=0, atom_indices=ai, ref_atom_indices=ai)
    X = traj.xyz[:, ai, :] * 10.0          # (F, N, 3) Angstrom
    F, N, _ = X.shape
    mean = X.mean(axis=0)                    # (N,3)
    d = X - mean                             # fluctuations
    # covariance C[i,j] = <d_i . d_j> over frames
    # flatten to (F, N*3), compute per-pair dot
    cov = np.einsum('fia,fja->ij', d, d) / F      # (N,N)  <dri.drj>
    var = np.diag(cov).copy()
    denom = np.sqrt(np.outer(var, var))
    dc = cov / denom                          # normalized DCCM in [-1,1]

    lm = cc.landmark_rows(rows)
    sub = cc.subunit_ids(rows)
    order, nxt, prv, cents = cc.ring_neighbor_order(ref, rows)
    wA = lm['walkerA']

    def block(mat, a, b):
        return float(mat[np.ix_(a, b)].mean())
    def absblock(mat, a, b):
        return float(np.abs(mat[np.ix_(a, b)]).mean())

    inter_mask = sub[:, None] != sub[None, :]
    intra_mask = ~inter_mask
    np.fill_diagonal(intra_mask, False)

    # nearest-neighbour subunit block (whole-subunit)
    nn_pairs = [(k, nxt[k]) for k in range(5)]
    nn_sub = np.mean([block(dc, np.where(sub == a)[0], np.where(sub == b)[0]) for a, b in nn_pairs])
    nn_sub_abs = np.mean([absblock(dc, np.where(sub == a)[0], np.where(sub == b)[0]) for a, b in nn_pairs])

    # ATP-site nearest-neighbour block
    nn_atp = np.mean([block(dc, wA[a], wA[b]) for a, b in nn_pairs])
    nn_atp_abs = np.mean([absblock(dc, wA[a], wA[b]) for a, b in nn_pairs])

    res = dict(
        kind=kind, n_frames=int(F), N=int(N),
        dccm_inter_all_mean=round(float(dc[inter_mask].mean()), 5),
        dccm_inter_all_absmean=round(float(np.abs(dc[inter_mask]).mean()), 5),
        dccm_intra_all_absmean=round(float(np.abs(dc[intra_mask]).mean()), 5),
        dccm_NN_subunit_mean=round(float(nn_sub), 5),
        dccm_NN_subunit_absmean=round(float(nn_sub_abs), 5),
        dccm_NN_ATPsite_mean=round(float(nn_atp), 5),
        dccm_NN_ATPsite_absmean=round(float(nn_atp_abs), 5),
    )
    res['dccm_inter_over_intra_abs'] = round(res['dccm_inter_all_absmean'] /
                                             res['dccm_intra_all_absmean'], 4)
    return res, dc, rows


def main():
    out = {}
    mats = {}
    for kind in ('design', 'native'):
        print(f"[{kind}] computing DCCM ...", flush=True)
        r, dc, rows = analyse = dccm(kind)
        out[kind] = r
        mats[kind] = dc
        print(f"[{kind}] frames={r['n_frames']}", flush=True)
    d, n = out['design'], out['native']
    out['comparison'] = {
        'dccm_NN_ATPsite_absmean_design_vs_native':
            [d['dccm_NN_ATPsite_absmean'], n['dccm_NN_ATPsite_absmean']],
        'dccm_NN_subunit_absmean_design_vs_native':
            [d['dccm_NN_subunit_absmean'], n['dccm_NN_subunit_absmean']],
        'dccm_inter_all_absmean_design_vs_native':
            [d['dccm_inter_all_absmean'], n['dccm_inter_all_absmean']],
        'dccm_inter_over_intra_abs_design_vs_native':
            [d['dccm_inter_over_intra_abs'], n['dccm_inter_over_intra_abs']],
    }
    np.save('gp16_design/outputs/coordination/dccm_design.npy', mats['design'])
    np.save('gp16_design/outputs/coordination/dccm_native.npy', mats['native'])
    with open('gp16_design/outputs/coordination/md_crosscorr_results.json', 'w') as fh:
        json.dump(out, fh, indent=2)
    print("\n==== MD DCCM inter-subunit summary ====")
    print(json.dumps(out['comparison'], indent=2))


if __name__ == '__main__':
    main()
