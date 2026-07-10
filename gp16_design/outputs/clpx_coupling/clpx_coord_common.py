#!/usr/bin/env python3
"""Shared residue mapping + landmark selection for the ClpX coordination /
force-transmission analysis -- the ClpX generalization of the gp16 cp233 study
(gp16_design/outputs/coordination/). Question: is the covalent single-chain ClpX
more inter-subunit coordinated than the native non-covalent hexamer, as cp233 was
vs native gp16 (PRS 1.79x)?

Structures
----------
  DESIGN  : single-chain ClpX(dN) "good" pseudohexamer, OpenFold3 prediction, apo.
            gp16_design/cs_sessions/session2_denovo_clpxvalidation/clpx_out/clpx_good__of3.cif
            One chain A, resSeq 1..2378 = 6x ClpXdN(native 62..424, 363 res) joined by
            5x (GGGGS)x8 = 40-res tethers. Copy spans (construct 1-based) below; tether
            residues are DROPPED (no native counterpart). Mapping verified by aligning the
            Walker-A P-loop motif GPTGSGKT: native_res = construct_pos - span_lo + 62.
  NATIVE  : E. coli ClpXP cryo-EM hexamer 6PP5 (substrate-engaged spiral, 5x ATPgS + 1x ADP
            seam subunit F). gp16_design/cs_sessions/session3_clpx_vs_gp16_topology/6PP5.cif
            6 protein chains A-F (UniProt P0A6H1 numbering), substrate peptide chain G,
            nucleotides AGS/ADP, Mg.

Both rings are reduced to the SAME residue set (native residues present in ALL 6 native
chains AND all 6 design copies, within the 64..413 window) so the ENM numbers are
directly comparable, exactly as in the gp16 study.

Landmarks (native ClpX / UniProt P0A6H1 numbering)
--------------------------------------------------
  Walker-A / ATP P-loop : 119..128  (GPTGSGKT..; catalytic K125)  -> the ATP site
  trans arginine finger : R307      (supplied in trans into the neighbour ATP site;
                                      direct ClpX analog of gp16's R146)
  sensor-II             : R370
  pore-1 loop aromatic  : Y153      (GYVG substrate-gripping paddle; ClpX analog of
                                      gp16's force residue Y129)
  substrate-contact set : residues <8 A from the threaded substrate (chain G) in >=3 of
                          the 6 native subunits -> pore loops 1 (149-157) and 2 (199-202)

CAVEATS (honest, prominent)
  * NO MD available for ClpX -> ENM/PRS + GNM only (the gp16 numbers had a short-MD
    cross-check; here we rely on the GNM-correlation network, the gp16 no-MD path).
  * State mismatch: DESIGN is an apo OpenFold3 prediction that folds to a symmetric
    compact ring (radius_cv 0.003); NATIVE 6PP5 is a substrate-engaged SPIRAL with a
    disengaged ADP seam (F). Part of any single-chain "coordination" advantage may be
    the native spiral seam being an intrinsic discontinuity. Flagged in the report.
  * Native seam subunit F is heavily disordered (pore-1 loop 148-154 and 192-200
    unmodeled); those loops are dropped from the strict common set for BOTH structures.
"""
import numpy as np
import prody
prody.confProDy(verbosity='none')

# ---- landmark residue numbers (native ClpX numbering) ----
WALKERA = list(range(119, 129))         # 119..128 ATP P-loop
ARGFINGER = 307                         # trans arginine finger (gp16 R146 analog)
SENSOR2 = 370
YFORCE = 153                            # pore-1 aromatic, substrate grip (gp16 Y129 analog)
# substrate-contact residues from 6PP5 (chain G proximity, >=3 subunits); the pore "grip"
PORE_CONTACT = [149, 151, 152, 153, 154, 155, 157, 199, 200, 202]
GUAN = ('NE', 'CZ', 'NH1', 'NH2')       # arginine guanidinium atoms
RES_LO, RES_HI = 64, 413                # window modeled across all native chains
NSUB = 6

# ---- design copy spans (construct resSeq, 1-based) ----
DESIGN_SPANS = [(1, 363), (404, 766), (807, 1169),
                (1210, 1572), (1613, 1975), (2016, 2378)]
NATIVE_COPY_LEN = 363                   # native 62..424

DESIGN_CIF = ('gp16_design/cs_sessions/session2_denovo_clpxvalidation/'
              'clpx_out/clpx_good__of3.cif')
NATIVE_CIF = ('gp16_design/cs_sessions/session3_clpx_vs_gp16_topology/6PP5.cif')
NATIVE_CHAINS = list('ABCDEF')


def design_construct_pos(subunit, native_res):
    """design chain-A resSeq for copy `subunit` (0..5) and native residue number."""
    lo = DESIGN_SPANS[subunit][0]
    return lo + (native_res - 62)


def _ca_maps(kind):
    """Return per-subunit {native_res: (ca_xyz, {atomname: xyz})} for landmark side-chain
    atoms, plus the set of native residues modeled in every subunit."""
    if kind == 'design':
        st = prody.parseMMCIF(DESIGN_CIF)
        sub_maps = []
        for k in range(NSUB):
            ca = {}
            heavy = {}
            for nr in range(RES_LO, RES_HI + 1):
                cp = design_construct_pos(k, nr)
                sel = st.select(f'chain A and resnum {cp}')
                if sel is None:
                    continue
                caat = sel.select('name CA')
                if caat is None:
                    continue
                ca[nr] = caat.getCoords()[0]
                heavy[nr] = {a.getName(): a.getCoords()
                             for a in sel.select('not hydrogen').iterAtoms()}
            sub_maps.append((ca, heavy))
    else:
        st = prody.parseMMCIF(NATIVE_CIF)
        sub_maps = []
        for chid in NATIVE_CHAINS:
            ca = {}
            heavy = {}
            for nr in range(RES_LO, RES_HI + 1):
                sel = st.select(f'chain {chid} and resnum {nr} and protein')
                if sel is None:
                    continue
                caat = sel.select('name CA')
                if caat is None:
                    continue
                ca[nr] = caat.getCoords()[0]
                heavy[nr] = {a.getName(): a.getCoords()
                             for a in sel.select('not hydrogen').iterAtoms()}
            sub_maps.append((ca, heavy))
    # residues modeled in every subunit
    common = set(range(RES_LO, RES_HI + 1))
    for ca, _ in sub_maps:
        common &= set(ca.keys())
    return sub_maps, common


def load_common(kind, common_override=None):
    """Load a structure into reduced CA arrays on the common residue set.

    Returns dict with:
      xyz     : (N,3) Angstrom CA coords, ordered (subunit ascending, native_res ascending)
      rows    : list of (row, subunit, native_res)
      sub     : (N,) subunit id per row
      natres  : (N,) native residue number per row
      heavy   : per-subunit {native_res: {atomname: xyz}}   (for guanidinium distances)
      common  : sorted list of native residues used
    """
    sub_maps, common = _ca_maps(kind)
    if common_override is not None:
        common = set(common_override) & common
    common_sorted = sorted(common)
    rows = []
    xyz = []
    r = 0
    for k in range(NSUB):
        ca, _ = sub_maps[k]
        for nr in common_sorted:
            rows.append((r, k, nr))
            xyz.append(ca[nr])
            r += 1
    return dict(
        xyz=np.asarray(xyz, float),
        rows=rows,
        sub=np.array([k for (_, k, _) in rows], int),
        natres=np.array([nr for (_, _, nr) in rows], int),
        heavy=[h for _, h in sub_maps],
        common=common_sorted,
    )


def landmark_rows(D):
    """Row indices for each landmark set, per subunit, given a load_common() dict."""
    lut = {(k, nr): r for (r, k, nr) in D['rows']}
    common = set(D['common'])
    wa = {k: [lut[(k, n)] for n in WALKERA if n in common] for k in range(NSUB)}
    arg = {k: lut[(k, ARGFINGER)] for k in range(NSUB)} if ARGFINGER in common else {}
    yf = {k: lut[(k, YFORCE)] for k in range(NSUB)} if YFORCE in common else {}
    pore = {k: [lut[(k, n)] for n in PORE_CONTACT if n in common] for k in range(NSUB)}
    return dict(walkerA=wa, arg=arg, yforce=yf, pore=pore, lut=lut)


def ring_neighbor_order(D):
    """Spatial ring order of the 6 subunits from Walker-A centroids (projection onto the
    best-fit ring plane, sorted by angle). Returns order, next-map, prev-map, centroids."""
    lm = landmark_rows(D)
    xyz = D['xyz']
    cents = np.array([xyz[lm['walkerA'][k]].mean(0) for k in range(NSUB)])
    v = cents - cents.mean(0)
    _, _, vt = np.linalg.svd(v - v.mean(0))
    basis = vt[:2]
    ang = np.arctan2(v @ basis[1], v @ basis[0])
    order = list(np.argsort(ang))
    nxt = {order[i]: order[(i + 1) % NSUB] for i in range(NSUB)}
    prv = {order[i]: order[(i - 1) % NSUB] for i in range(NSUB)}
    return order, nxt, prv, cents


def functional_neighbor(D):
    """For each subunit k, the neighbour subunit m whose Walker-A centroid is closest to
    k's arginine-finger (R307) guanidinium -- the trans arg-finger donor->acceptor
    direction (k feeds m's ATP site). Falls back to R307 CA if guanidinium absent."""
    lm = landmark_rows(D)
    xyz = D['xyz']
    wcent = {k: xyz[lm['walkerA'][k]].mean(0) for k in range(NSUB)}
    nb = {}
    for k in range(NSUB):
        hv = D['heavy'][k].get(ARGFINGER, {})
        gpts = [hv[a] for a in GUAN if a in hv]
        if gpts:
            src = np.array(gpts)
        else:
            src = xyz[[lm['arg'][k]]] if lm['arg'] else wcent[k][None]
        best, bm = 1e9, None
        for m in range(NSUB):
            if m == k:
                continue
            d = np.linalg.norm(src[:, None, :] - wcent[m][None, None, :], axis=-1).min()
            if d < best:
                best, bm = d, m
        nb[k] = (bm, round(float(best), 2))
    return nb


if __name__ == '__main__':
    for kind in ('design', 'native'):
        D = load_common(kind)
        print(f"\n{kind}: {len(D['rows'])} CA in reduced set "
              f"({len(D['common'])} residues x {NSUB} subunits)")
        print(f"  common residue window: {D['common'][0]}..{D['common'][-1]} "
              f"(n={len(D['common'])}); Walker-A present: "
              f"{[n for n in WALKERA if n in set(D['common'])]}")
        order, nxt, prv, cents = ring_neighbor_order(D)
        print(f"  spatial ring order (Walker-A centroid angle): {order}")
        fnb = functional_neighbor(D)
        print(f"  arg-finger(R307) functional neighbour map k->m (dist A): "
              f"{ {k: fnb[k] for k in fnb} }")
