#!/usr/bin/env python3
"""
Aim C2 — ENM/ANM soft-mode overlap analysis for the gp16 P->H->P cycle.

INDEPENDENT implementation using ProDy (parsePDB / ANM / calcANM), on the
canonical validation inputs:
  P (planar)  = md/openmm_validation/inputs/A_apo.pdb            (native apo ring, 5 chains A-E, res 1-332)
  H (helical) = md/openmm_validation/inputs/B_7jqq_helical.pdb   (7JQQ staircase, 5 chains A-E, res 4-330)
  design      = md/openmm_validation/trajectories/C/C_start.pdb  (cp233 single chain, res 1-1750)

Question: is the P<->H transition a low-frequency (soft) collective mode of the
ring, and did the covalent single-chain DESIGN keep that soft mode?

Outputs (this dir):
  C2_enm.npz            Delta vector, native+design ANM eigvals/eigvecs, metadata
  C2_overlap_figure.png overlap bars + cumulative curve + eigenvalue spectrum
Console: full overlap table + native-vs-design mode-character comparison.
"""
import os, numpy as np
import prody
prody.confProDy(verbosity='none')

REPO = "/Users/longfu/Developer/claude-science-hackthon/.claude/worktrees/export-gp16-results/gp16_design"
OUT  = f"{REPO}/outputs/php_cycle/C2_enm"
A_APO = f"{REPO}/md/openmm_validation/inputs/A_apo.pdb"
B_JQQ = f"{REPO}/md/openmm_validation/inputs/B_7jqq_helical.pdb"
C_DES = f"{REPO}/md/openmm_validation/trajectories/C/C_start.pdb"
CUT   = 13.0          # match repo tools (reproduce/anm.py, overlap.py, enm_openclose.py)
NMODE = 20            # non-trivial modes to keep
CHAINS = list("ABCDE")
DES_COPIES = [(1,342),(353,694),(705,1046),(1057,1398),(1409,1750)]  # 5 subunits in single chain

def kabsch(P, Q):
    """Superpose P onto Q; return aligned P and RMSD."""
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ Wt)); D = np.diag([1, 1, d])
    R = V @ D @ Wt
    Pal = Pc @ R + Q.mean(0)
    return Pal, float(np.sqrt(((Pal - Q) ** 2).sum(1).mean()))

def ca_by_chain_res(path, resnums=None):
    """{chain: {resnum: xyz}} for CA atoms (parsed independently with ProDy)."""
    ag = prody.parsePDB(path, subset='ca')
    d = {}
    for at in ag:
        c = at.getChid(); r = int(at.getResnum())
        if resnums is not None and r not in resnums: continue
        d.setdefault(c, {})[r] = at.getCoords()
    return d

def ring_frame(cent):
    """Best-fit plane normal from 5 subunit centroids; radial unit vecs."""
    ctr = cent.mean(0)
    _, _, Vt = np.linalg.svd(cent - ctr); normal = Vt[2]
    rad = []
    for p in cent:
        v = (p - ctr) - ((p - ctr) @ normal) * normal
        rad.append(v / (np.linalg.norm(v) + 1e-12))
    return ctr, normal, np.array(rad)

# ---------------------------------------------------------------------------
# PART 1 : native P -> H difference vector
# ---------------------------------------------------------------------------
print("="*74)
print("PART 1 : native P(apo) -> H(7JQQ) difference vector")
print("="*74)
apo = ca_by_chain_res(A_APO)
jqq = ca_by_chain_res(B_JQQ)
common = sorted(set.intersection(*[set(apo[c]) for c in CHAINS],
                                 *[set(jqq[c]) for c in CHAINS]))
print(f"common residues per subunit: {len(common)} (res {common[0]}-{common[-1]}); "
      f"ring = {5*len(common)} CA")

# P ring in fixed order chain A..E, res=common
Xp = np.array([apo[c][r] for c in CHAINS for r in common])

# find best chain assignment of H (cyclic rotations + reflection) by min RMSD
best = None
cand = [CHAINS[i:]+CHAINS[:i] for i in range(5)] + \
       [(CHAINS[::-1])[i:]+(CHAINS[::-1])[:i] for i in range(5)]
for order in cand:
    Xh = np.array([jqq[order[k]][r] for k in range(5) for r in common])
    Xh_al, rmsd = kabsch(Xh, Xp)
    if best is None or rmsd < best[0]:
        best = (rmsd, order, Xh_al)
rmsd, order, Xh_al = best
print(f"best chain assignment apo(ABCDE) -> 7JQQ{order}")
print(f"whole-ring apo<->7JQQ Ca-RMSD (after Kabsch) = {rmsd:.2f} A")

Delta = (Xh_al - Xp)                      # H - P per-CA (matched, aligned)
Delta_flat = Delta.ravel()
Delta_hat = Delta_flat / np.linalg.norm(Delta_flat)

# --- axial sanity: intrinsic subunit-centroid axial spread in each state's own frame
def axial_spread(coords_by_subunit):
    cent = np.array([c.mean(0) for c in coords_by_subunit])
    _, n, _ = ring_frame(cent)
    ax = (cent - cent.mean(0)) @ n
    return float(ax.std()), float(ax.max()-ax.min()), ax
Xp_sub = [Xp[k*len(common):(k+1)*len(common)] for k in range(5)]
Xh_sub = [Xh_al[k*len(common):(k+1)*len(common)] for k in range(5)]
p_std, p_pp, _ = axial_spread(Xp_sub)
# helical spread from the *unaligned* H (its intrinsic staircase)
Xh_raw = np.array([jqq[order[k]][r] for k in range(5) for r in common])
Xh_raw_sub = [Xh_raw[k*len(common):(k+1)*len(common)] for k in range(5)]
h_std, h_pp, _ = axial_spread(Xh_raw_sub)
print(f"\nAXIAL SANITY (subunit-centroid spread in each state's own best-fit plane):")
print(f"  planar apo : std {p_std:.2f} A   peak-to-peak {p_pp:.2f} A")
print(f"  7JQQ (H)   : std {h_std:.2f} A   peak-to-peak {h_pp:.2f} A   (calib: staircase ~4.8 A vs planar ~0.1 A)")

# per-subunit axial component of Delta (apo frame)
ctr_p, normal_p, rad_p = ring_frame(np.array([s.mean(0) for s in Xp_sub]))
dsub = np.array([Delta[k*len(common):(k+1)*len(common)].mean(0) for k in range(5)])
dsub_ax = dsub @ normal_p
print(f"\nper-subunit MEAN axial component of Delta (H-P) along apo ring normal (A):")
for k, c in enumerate(CHAINS):
    print(f"    subunit {c}: {dsub_ax[k]:+.2f}")
print(f"    -> axial spread of the difference vector (std) = {dsub_ax.std():.2f} A "
      f"(staircase-forming motion present: {'YES' if dsub_ax.std()>0.5 else 'no'})")

# ---------------------------------------------------------------------------
# PART 2/3 : ANM on planar apo (common CA)  + overlap with Delta
# ---------------------------------------------------------------------------
print("\n" + "="*74)
print(f"PART 2/3 : ANM on planar apo (common {5*len(common)} CA, cutoff {CUT} A) + overlap")
print("="*74)
# build a coordinate-only AtomGroup in the SAME order as Xp / Delta
apo_ag = prody.AtomGroup('apo_common')
apo_ag.setCoords(Xp)
anm = prody.ANM('apo_common')
anm.buildHessian(Xp, cutoff=CUT, gamma=1.0)
anm.calcModes(n_modes=NMODE, zeros=False)     # excludes 6 trivial modes
eig_n = anm.getEigvals()
V_n = anm.getEigvecs()                          # (3N, NMODE), columns are unit modes

# manual overlaps (mode_i . Delta_hat); cross-check with ProDy calcOverlap
ov = np.abs(V_n.T @ Delta_hat)
ov_prody = np.abs(prody.calcOverlap(anm, prody.Vector(Delta_flat, 'PtoH')))
assert np.allclose(ov, ov_prody, atol=1e-6), "ProDy vs manual overlap mismatch"
cum = np.sqrt(np.cumsum(ov**2))                 # cumulative overlap (sqrt sum sq)

print(f"{'ProDy#':>6} {'repo#':>6} {'eigval':>9} {'1/eig(soft)':>11} {'|overlap|':>10}  {'cumOv':>6} {'cum%':>6}")
best_mode = int(np.argmax(ov))
for i in range(NMODE):
    tag = "  <== best single mode (P->H)" if i == best_mode else ""
    print(f"{i+1:>6} {i+7:>6} {eig_n[i]:>9.4f} {1.0/eig_n[i]:>11.3f} "
          f"{ov[i]:>10.3f}  {cum[i]:>6.3f} {100*cum[i]**2:>5.0f}%{tag}")

print(f"\nbest single ANM mode: ProDy#{best_mode+1} (repo#{best_mode+7}), "
      f"overlap {ov[best_mode]:.3f}")
for K in (1,3,5,10,20):
    print(f"  cumulative overlap, {K:>2} softest modes: {cum[K-1]:.3f}  "
          f"({100*cum[K-1]**2:.0f}% of the P->H transition captured)")

# cutoff-15 sensitivity (does the story survive a different cutoff?)
anm15 = prody.ANM('apo15'); anm15.buildHessian(Xp, cutoff=15.0, gamma=1.0)
anm15.calcModes(n_modes=NMODE, zeros=False)
ov15 = np.abs(anm15.getEigvecs().T @ Delta_hat); cum15 = np.sqrt(np.cumsum(ov15**2))
print(f"\ncutoff=15 A sensitivity: best single {ov15.max():.3f} (mode {int(np.argmax(ov15))+1}); "
      f"cum 5={cum15[4]:.3f}({100*cum15[4]**2:.0f}%) 10={cum15[9]:.3f}({100*cum15[9]**2:.0f}%) "
      f"20={cum15[19]:.3f}({100*cum15[19]**2:.0f}%)")

# ---------------------------------------------------------------------------
# mode character on planar apo (helical/out-of-plane vs radial-breathing)
# ---------------------------------------------------------------------------
def mode_character(V, eig, X, subunit_slices, normal, rad, label):
    print(f"\n{label}: per-mode subunit-centroid character (softest {V.shape[1]} modes)")
    print(f"{'ProDy#':>6} {'1/eig':>8} {'helical/oop':>11} {'radial':>8} {'local':>6}  character")
    soft0 = 1.0/eig[0]
    chars = []
    for m in range(V.shape[1]):
        v = V[:, m].reshape(-1, 3)
        dc = np.array([v[sl].mean(0) for sl in subunit_slices])
        mag = np.linalg.norm(dc, axis=1)
        oop = np.abs(dc @ normal)
        radc = np.array([(dc[i] - (dc[i] @ normal) * normal) @ rad[i] for i in range(len(dc))])
        helical = np.sqrt((oop**2).mean()) / (np.sqrt((mag**2).mean()) + 1e-9)
        breathing = abs(radc.mean()) / (np.sqrt((mag**2).mean()) + 1e-9)
        local = mag.max() / (mag.mean() + 1e-9)
        lab = ("HELICAL/oop" if helical > 0.5 else
               ("breathing" if breathing > 0.5 else "in-plane/other"))
        chars.append(helical)
        if m < 14:
            print(f"{m+1:>6} {(1.0/eig[m])/soft0:>7.2f}x {helical:>11.2f} "
                  f"{breathing:>8.2f} {local:>6.2f}  {lab}")
    return np.array(chars)
slices_p = [slice(k*len(common),(k+1)*len(common)) for k in range(5)]
helical_native = mode_character(V_n, eig_n, Xp, slices_p, normal_p, rad_p,
                                "PLANAR apo (native ring)")

# ---------------------------------------------------------------------------
# PART 4 : ANM on the DESIGN (single chain) + mode-character comparison
# ---------------------------------------------------------------------------
print("\n" + "="*74)
print(f"PART 4 : ANM on DESIGN cp233 single chain (cutoff {CUT} A)")
print("="*74)
des_ag = prody.parsePDB(C_DES, subset='ca')
Xd = des_ag.getCoords()
rd = np.array([int(a.getResnum()) for a in des_ag])
# subunit label per design CA from copy ranges (linker -> -1)
lab_d = -np.ones(len(rd), dtype=int)
for k,(lo,hi) in enumerate(DES_COPIES):
    lab_d[(rd>=lo)&(rd<=hi)] = k
print(f"design CA: {len(Xd)} ; per-subunit core counts: "
      + ", ".join(str(int((lab_d==k).sum())) for k in range(5))
      + f" ; linker CA excluded: {(lab_d<0).sum()}")

anm_d = prody.ANM('design'); anm_d.buildHessian(Xd, cutoff=CUT, gamma=1.0)
anm_d.calcModes(n_modes=NMODE, zeros=False)
eig_d = anm_d.getEigvals(); V_d = anm_d.getEigvecs()

# design ring frame from the 5 subunit cores
des_cent = np.array([Xd[lab_d==k].mean(0) for k in range(5)])
ctr_d, normal_d, rad_d = ring_frame(des_cent)
slices_d = [np.where(lab_d==k)[0] for k in range(5)]
des_ax_std, des_ax_pp, _ = axial_spread([Xd[lab_d==k] for k in range(5)])
print(f"design subunit-centroid axial spread: std {des_ax_std:.2f} A pp {des_ax_pp:.2f} A "
      f"({'planar' if des_ax_std<1 else 'non-planar'})")
helical_design = mode_character(V_d, eig_d, Xd, slices_d, normal_d, rad_d,
                                "DESIGN cp233 (single chain)")

# softest HELICAL mode rank in native vs design
def first_helical(chars, thr=0.5):
    idx = np.where(chars > thr)[0]
    return (int(idx[0])+1, float(chars[idx[0]])) if len(idx) else (None, None)
nh_i, nh_v = first_helical(helical_native)
dh_i, dh_v = first_helical(helical_design)
print(f"\nsoftest HELICAL/out-of-plane mode:")
print(f"  native : ProDy#{nh_i} (helical frac {nh_v:.2f})")
print(f"  design : ProDy#{dh_i} (helical frac {dh_v:.2f})")
n_help_top5 = int((helical_native[:5] > 0.5).sum())
d_help_top5 = int((helical_design[:5] > 0.5).sum())
print(f"  # helical modes among softest 5:  native {n_help_top5}/5   design {d_help_top5}/5")

# ---------------------------------------------------------------------------
# breathing content (enm_openclose metric) on the SAME two inputs
# ---------------------------------------------------------------------------
def breathing_vector(X, subunit_slices):
    cent = np.array([X[sl].mean(0) for sl in subunit_slices])
    ctr, normal, _ = ring_frame(cent)
    B = np.zeros(3*len(X))
    for i, p in enumerate(X):
        v = (p - ctr) - ((p - ctr) @ normal) * normal
        n = np.linalg.norm(v)
        if n > 1e-6: B[3*i:3*i+3] = v/n
    return B/np.linalg.norm(B)
def cum_breathing(V, B, ks=(1,3,5,10,20)):
    ov = V.T @ B
    c = np.cumsum(ov**2)
    return {k: float(c[min(k,len(c))-1]) for k in ks}
# native breathing on FULL apo ring (1660 CA) for parity with repo enm_openclose native
apo_full = prody.parsePDB(A_APO, subset='ca'); Xpf = apo_full.getCoords()
rpf = np.array([a.getChid() for a in apo_full])
slices_pf = [np.where(rpf==c)[0] for c in CHAINS]
anm_pf = prody.ANM('apo_full'); anm_pf.buildHessian(Xpf, cutoff=CUT, gamma=1.0)
anm_pf.calcModes(n_modes=NMODE, zeros=False)
B_n = breathing_vector(Xpf, slices_pf)
B_d = breathing_vector(Xd, slices_d)
cb_n = cum_breathing(anm_pf.getEigvecs(), B_n)
cb_d = cum_breathing(V_d, B_d)
print(f"\nbreathing (radial ring-opening) content captured by softest modes (indep. ProDy):")
print(f"  native apo (1660 CA): " + " ".join(f"top{k}={cb_n[k]:.2f}" for k in (1,3,5,10,20)))
print(f"  design C_start (1750): " + " ".join(f"top{k}={cb_d[k]:.2f}" for k in (1,3,5,10,20)))

# ---------------------------------------------------------------------------
# save npz + figure
# ---------------------------------------------------------------------------
np.savez_compressed(f"{OUT}/C2_enm.npz",
    Delta=Delta, Delta_hat=Delta_hat, common_res=np.array(common),
    chain_order=np.array(order), ring_rmsd=rmsd,
    native_eigvals=eig_n, native_eigvecs=V_n, native_overlap=ov, native_cumoverlap=cum,
    native_eigvals_c15=anm15.getEigvals(), native_overlap_c15=ov15,
    design_eigvals=eig_d, design_eigvecs=V_d,
    helical_native=helical_native, helical_design=helical_design,
    breathing_native=np.array([cb_n[k] for k in (1,3,5,10,20)]),
    breathing_design=np.array([cb_d[k] for k in (1,3,5,10,20)]),
    cutoff=CUT, best_mode_prody=best_mode+1, Xp=Xp, Xh_aligned=Xh_al)
print(f"\nsaved {OUT}/C2_enm.npz")

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
x = np.arange(1, NMODE+1)
# panel 1: per-mode overlap + cumulative
bars = ax[0].bar(x, ov, color='#4C78A8', label='per-mode |overlap|')
bars[best_mode].set_color('#E45756')
ax0b = ax[0].twinx()
ax0b.plot(x, 100*cum**2, '-o', color='#222', ms=4, lw=1.6, label='cumulative %')
ax[0].set_xlabel('ANM mode (ProDy#, non-trivial)'); ax[0].set_ylabel('|overlap| with P->H')
ax0b.set_ylabel('cumulative % of P->H captured'); ax0b.set_ylim(0, 60)
ax[0].set_title(f'P->H overlap  (best single = ProDy#{best_mode+1}, {ov[best_mode]:.2f})')
ax[0].set_xticks(x[::2])
# panel 2: eigenvalue spectrum native vs design (normalized to softest)
ax[1].plot(x, eig_n/eig_n[0], '-o', color='#4C78A8', ms=4, label='native apo')
ax[1].plot(x, eig_d/eig_d[0], '-s', color='#F58518', ms=4, label='design cp233')
ax[1].set_xlabel('ANM mode (ProDy#)'); ax[1].set_ylabel('eigenvalue / softest eigenvalue')
ax[1].set_title('low-mode stiffness spectrum'); ax[1].legend(); ax[1].set_xticks(x[::2])
# panel 3: helical/out-of-plane character native vs design
w = 0.4
ax[2].bar(x-w/2, helical_native, w, color='#4C78A8', label='native apo')
ax[2].bar(x+w/2, helical_design, w, color='#F58518', label='design cp233')
ax[2].axhline(0.5, color='#888', ls='--', lw=1)
ax[2].set_xlabel('ANM mode (ProDy#)'); ax[2].set_ylabel('helical / out-of-plane fraction')
ax[2].set_title('planar<->helical character of soft modes'); ax[2].legend(); ax[2].set_xticks(x[::2])
plt.tight_layout()
plt.savefig(f"{OUT}/C2_overlap_figure.png", dpi=150)
print(f"saved {OUT}/C2_overlap_figure.png")
print("\nDONE.")
