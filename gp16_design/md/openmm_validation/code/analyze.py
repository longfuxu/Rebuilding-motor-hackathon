#!/usr/bin/env python3
"""Trajectory readouts for the gp16 ring MD, consistent with reproduce/score_m2.py.

Per frame: Ca-RMSD (vs prod frame 0), ring radius / radius_CV / planarity_rms,
per-sequential-interface R146(guanidinium)->neighbour Walker-A(24-31) min-dist,
n_engaged(<8A), and per-interface heavy-atom residue-contact count (pairs <4.5A,
over the t0 interface residue set). Plus per-residue RMSF.

Subunit + adjacency (LOCKED from t=0 so each interface tracks the same seam):
  A (apo, chains A..E):      donor i -> acceptor (i+1)%5   [A->B->C->D->E->A]
  B (7JQQ, chains A..E):     donor i -> acceptor (i-1)%5   [A->E,B->A,C->B,D->C,E->D]
  C (design, 1 chain,copies):donor i -> acceptor (i+1)%5   [copy1->2->..->5->1]
"""
import argparse, os, json
import numpy as np
import mdtraj as md

GUAN = ('NE','CZ','NH1','NH2')
WALKER = list(range(24,32))       # 24..31 (gp16 numbering)
R146 = 146
ENGAGED_A = 8.0

# design copy layout (absolute resSeq): copy start lo -> R146=lo+254, Walker=lo+132..lo+139
COPY_LOS = [1,353,705,1057,1409]

def subunits_and_adj(system, top):
    """Return (subs, adj) where subs[k] = dict(r146_guan_atomidx, walker_atomidx, ca_atomidx),
    adj[k] = acceptor subunit index."""
    subs = []
    if system in ('A','B'):
        chains = list(top.chains)
        assert len(chains) >= 5, f"expected >=5 chains, got {len(chains)}"
        for ch in chains[:5]:
            rs = {r.resSeq: r for r in ch.residues}
            r146 = rs.get(R146)
            guan = [a.index for a in r146.atoms if a.name in GUAN] if r146 else []
            walk = [a.index for rn in WALKER if rn in rs for a in rs[rn].atoms
                    if a.element is not None and a.element.symbol != 'H']
            ca = [a.index for r in ch.residues for a in r.atoms if a.name == 'CA']
            subs.append(dict(guan=guan, walk=walk, ca=ca,
                             resset={rn: [a.index for a in rs[rn].atoms
                                          if a.element.symbol!='H'] for rn in rs}))
        adj = [(i+1) % 5 for i in range(5)] if system == 'A' else [(i-1) % 5 for i in range(5)]
    else:  # C
        ch = list(top.chains)[0]
        rs = {r.resSeq: r for r in ch.residues}
        for lo in COPY_LOS:
            r146abs = lo + 254
            walkabs = [lo + o for o in range(132,140)]   # lo+132..lo+139
            hi = lo + 341                                  # copy span lo..lo+341
            r = rs.get(r146abs)
            guan = [a.index for a in r.atoms if a.name in GUAN] if r else []
            walk = [a.index for rn in walkabs if rn in rs for a in rs[rn].atoms
                    if a.element.symbol != 'H']
            ca = [a.index for rn in range(lo, hi+1) if rn in rs
                  for a in rs[rn].atoms if a.name == 'CA']
            resset = {rn: [a.index for a in rs[rn].atoms if a.element.symbol!='H']
                      for rn in range(lo, hi+1) if rn in rs}
            subs.append(dict(guan=guan, walk=walk, ca=ca, resset=resset))
        adj = [(i+1) % 5 for i in range(5)]
    return subs, adj

def ring_geometry(cent):
    C = np.asarray(cent); k = len(C)
    ctr = C.mean(0); X = C - ctr
    _,_,Vt = np.linalg.svd(X); normal = Vt[2]
    oop = X @ normal
    planarity = float(np.sqrt((oop**2).mean()))
    inplane = X - np.outer(oop, normal)
    radii = np.linalg.norm(inplane, axis=1)
    return float(radii.mean()), float(radii.std()/radii.mean()), planarity

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--system', required=True, choices=['A','B','C'])
    ap.add_argument('--top', required=True)      # <sys>_start.pdb
    ap.add_argument('--traj', required=True)     # <sys>_prod.dcd
    ap.add_argument('--out', required=True)
    ap.add_argument('--report_ps', type=float, default=20.0)
    ap.add_argument('--contact_stride', type=int, default=5)
    args = ap.parse_args()

    top = md.load_topology(args.top)
    traj = md.load_dcd(args.traj, top=args.top)
    n = traj.n_frames
    print(f"[analyze] {args.system}: {n} frames, {traj.n_atoms} atoms", flush=True)
    subs, adj = subunits_and_adj(args.system, top)

    # --- interface R146->Walker atom-pair index lists (donor guan x acceptor walk) ---
    all_ca = np.array(sorted(a.index for ch in top.chains for r in ch.residues
                             for a in r.atoms if a.name == 'CA'))

    # --- per-frame CA-RMSD vs prod frame 0 ---
    traj.superpose(traj, frame=0, atom_indices=all_ca)
    ca_rmsd = md.rmsd(traj, traj, frame=0, atom_indices=all_ca) * 10.0  # nm->A

    # --- RMSF (per CA, after superpose) ---
    rmsf = md.rmsf(traj, traj, frame=0, atom_indices=all_ca) * 10.0

    # --- interface distances per frame ---
    idist = np.full((n,5), np.nan)
    for k in range(5):
        g = subs[k]['guan']; w = subs[adj[k]]['walk']
        if not g or not w:
            continue
        pairs = np.array([(gi,wi) for gi in g for wi in w])
        d = md.compute_distances(traj, pairs) * 10.0  # A
        idist[:,k] = d.min(axis=1)
    n_eng = (idist < ENGAGED_A).sum(axis=1)

    # --- ring geometry per frame (subunit CA centroids) ---
    radius = np.zeros(n); radcv = np.zeros(n); plan = np.zeros(n)
    xyz = traj.xyz * 10.0  # A
    ca_lists = [np.array(subs[k]['ca']) for k in range(5)]
    for f in range(n):
        cent = [xyz[f, ca_lists[k]].mean(0) for k in range(5)]
        radius[f], radcv[f], plan[f] = ring_geometry(cent)

    # --- interface contacts (residue closest-heavy <4.5A) over t0 interface residue set ---
    contact_stride = max(1, args.contact_stride)
    fsel = np.arange(0, n, contact_stride)
    ncontact = np.full((len(fsel),5), np.nan)
    # build per-interface residue pair list: donor res x acceptor res with min-heavy<10A at frame0
    for k in range(5):
        dres = subs[k]['resset']; ares = subs[adj[k]]['resset']
        dkeys = list(dres); akeys = list(ares)
        # coarse prefilter by CA distance at frame0 to keep pair list small
        dca = {rn:[i for i in dres[rn]] for rn in dkeys}
        # map residue -> its CA index
        def ca_of(resset_map, rn, top=top):
            for a in top.atom(resset_map[rn][0]).residue.atoms:
                if a.name=='CA': return a.index
            return resset_map[rn][0]
        dca_idx = {rn: ca_of(dres, rn) for rn in dkeys}
        aca_idx = {rn: ca_of(ares, rn) for rn in akeys}
        f0 = xyz[0]
        pairlist = []
        for rd in dkeys:
            pd = f0[dca_idx[rd]]
            for ra in akeys:
                if np.linalg.norm(pd - f0[aca_idx[ra]]) < 14.0:
                    pairlist.append((rd, ra))
        if not pairlist:
            ncontact[:,k] = 0; continue
        respairs = np.array([[top.atom(dres[rd][0]).residue.index,
                              top.atom(ares[ra][0]).residue.index] for rd,ra in pairlist])
        sub = traj[fsel]
        dmat = md.compute_contacts(sub, contacts=respairs, scheme='closest-heavy')[0]*10.0
        ncontact[:,k] = (dmat < 4.5).sum(axis=1)

    # --- write timeseries ---
    t_ns = np.arange(n) * args.report_ps / 1000.0
    os.makedirs(args.out, exist_ok=True)
    ts_path = os.path.join(args.out, f'{args.system}_timeseries.csv')
    with open(ts_path,'w') as fh:
        cols = (['frame','time_ns','ca_rmsd_A','radius_A','radius_CV','planarity_A','n_engaged']
                + [f'd_if{k+1}_A' for k in range(5)])
        fh.write(','.join(cols)+'\n')
        for f in range(n):
            row = [f, f"{t_ns[f]:.4f}", f"{ca_rmsd[f]:.3f}", f"{radius[f]:.3f}",
                   f"{radcv[f]:.4f}", f"{plan[f]:.3f}", int(n_eng[f])] + \
                  [f"{idist[f,k]:.3f}" for k in range(5)]
            fh.write(','.join(map(str,row))+'\n')

    ct_path = os.path.join(args.out, f'{args.system}_contacts.csv')
    with open(ct_path,'w') as fh:
        fh.write('frame,time_ns,'+','.join(f'nc_if{k+1}' for k in range(5))+',nc_total\n')
        for i,f in enumerate(fsel):
            vals = [int(ncontact[i,k]) if not np.isnan(ncontact[i,k]) else '' for k in range(5)]
            tot = int(np.nansum(ncontact[i]))
            fh.write(f"{f},{t_ns[f]:.4f},"+','.join(map(str,vals))+f",{tot}\n")

    rmsf_path = os.path.join(args.out, f'{args.system}_rmsf.csv')
    ca_atoms = [top.atom(i) for i in all_ca]
    with open(rmsf_path,'w') as fh:
        fh.write('idx,chain,resSeq,resname,rmsf_A\n')
        for i,a in enumerate(ca_atoms):
            fh.write(f"{i},{a.residue.chain.index},{a.residue.resSeq},{a.residue.name},{rmsf[i]:.3f}\n")

    summary = dict(
        system=args.system, n_frames=int(n), t_total_ns=float(t_ns[-1]) if n else 0.0,
        ca_rmsd_final=float(ca_rmsd[-1]), ca_rmsd_mean_last20pct=float(ca_rmsd[int(0.8*n):].mean()),
        iface_dist_t0=[float(idist[0,k]) for k in range(5)],
        iface_dist_final=[float(idist[-1,k]) for k in range(5)],
        iface_dist_mean_last20pct=[float(np.nanmean(idist[int(0.8*n):,k])) for k in range(5)],
        n_engaged_t0=int(n_eng[0]), n_engaged_final=int(n_eng[-1]),
        n_engaged_mean_last20pct=float(n_eng[int(0.8*n):].mean()),
        radius_CV_t0=float(radcv[0]), radius_CV_final=float(radcv[-1]),
        radius_CV_mean_last20pct=float(radcv[int(0.8*n):].mean()),
        planarity_t0=float(plan[0]), planarity_final=float(plan[-1]),
        rmsf_mean=float(rmsf.mean()), rmsf_max=float(rmsf.max()),
    )
    with open(os.path.join(args.out, f'{args.system}_summary.json'),'w') as fh:
        json.dump(summary, fh, indent=2)
    print(f"[analyze] t0 iface (should match score_m2): "
          f"{['%.2f'%idist[0,k] for k in range(5)]}", flush=True)
    print(f"[analyze] final iface: {['%.2f'%idist[-1,k] for k in range(5)]}  "
          f"n_eng {n_eng[0]}->{n_eng[-1]}  radiusCV {radcv[0]:.3f}->{radcv[-1]:.3f}  "
          f"RMSD_final {ca_rmsd[-1]:.2f}A", flush=True)
    print("[analyze] wrote", ts_path, ct_path, rmsf_path, flush=True)

if __name__ == '__main__':
    main()
