# Aim C2 — ENM / ANM soft-mode overlap analysis (gp16 P→H→P cycle)

**Question (necessary-condition test).** The P→H→P model says the φ29 gp16 ring cycles
between a **planar** ring (P) and a **helical staircase** (H). Because both AF3 and Boltz
predict the ring planar for apo/ATP/ADP (a predictor-basin bias), H-state accessibility
has to be probed by *physics*, not prediction. The cheapest such probe is an elastic
network model: **is the P↔H transition a low-frequency (soft) collective mode of the ring
architecture, and did the covalent single-chain *design* keep that soft mode, or did the
linkers stiffen it out?**

**Method (this analysis is an INDEPENDENT re-derivation).** ProDy 2.6.1 (`parsePDB`,
`ANM.buildHessian`/`calcModes`, `calcOverlap`) — a different implementation and different
input structures from the repo tools it is cross-checked against. Cα-only ANM, cutoff
13 Å, γ = 1, lowest 20 non-trivial modes (the 6 trivial rigid-body modes are excluded).
Script: `c2_enm_analysis.py` (this dir). All numbers below are the ProDy values.

Inputs:
| role | file | Cα |
|---|---|---|
| P (planar) | `md/openmm_validation/inputs/A_apo.pdb` | 5 chains, 1660 (res 1–332) |
| H (helical) | `md/openmm_validation/inputs/B_7jqq_helical.pdb` | 5 chains, 1635 (res 4–330) |
| design | `md/openmm_validation/trajectories/C/C_start.pdb` | 1 chain, 1750 (5 copies) |

> **Mode indexing.** ProDy returns *non-trivial* modes, so **ProDy# = repo#(1-indexed) − 6**.
> ProDy#4 = the repo's "mode 10". Both columns are shown in the table.

---

## 1. P→H difference vector (native) — the target motion is a real staircase

Mapped A_apo ↔ B_7jqq on the common (chain, resSeq) Cα (**327 residues/subunit, res 4–330,
1635-Cα ring**), best cyclic+reflected chain assignment apo(ABCDE) → 7JQQ **(B,A,E,D,C)**,
Kabsch-superposed. Δ = H − P.

- **Whole-ring apo↔7JQQ Cα-RMSD = 6.65 Å** (repo cross-check 6.6 Å ✓).
- **Axial sanity (subunit-centroid spread, each state's own best-fit plane):**
  planar apo **std 0.05 Å / peak-to-peak 0.14 Å**; 7JQQ **std 1.80 Å / peak-to-peak 4.76 Å**
  → reproduces the calibration (staircase ≈ 4.8 Å vs planar ≈ 0.1 Å ✓). The H reference is
  a genuine axial staircase, not a distorted planar ring.
- **Per-subunit mean axial component of Δ** (along the apo ring normal), Å:
  A −0.26, B −1.48, C −3.87, **D +4.48**, E +1.13 → axial std **2.78 Å**. The difference
  vector *is* a staircase-forming (out-of-plane, rank-ordered) displacement — exactly the P→H
  motion the CV in C3 must drive.

## 2/3. ANM on the planar ring + overlap with P→H

ANM on planar apo (1635 common Cα, cutoff 13 Å). |overlap| = |mode·Δ̂|; cumulative overlap =
√(Σ overlapᵢ²); "cum %" = cumulative² = fraction of the P→H *variance* captured.

| ProDy# | repo# | eigval | 1/eig (soft) | \|overlap\| | cumOv | cum % |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 7 | 0.0897 | 11.15 | 0.100 | 0.100 | 1% |
| 2 | 8 | 0.0918 | 10.90 | 0.046 | 0.110 | 1% |
| 3 | 9 | 0.1178 | 8.49 | 0.244 | 0.267 | 7% |
| **4** | **10** | **0.1286** | **7.78** | **0.329** | **0.424** | **18%** | ← best single mode |
| 5 | 11 | 0.1301 | 7.69 | 0.122 | 0.441 | 19% |
| 6 | 12 | 0.1583 | 6.32 | 0.029 | 0.442 | 20% |
| 7 | 13 | 0.1596 | 6.27 | 0.036 | 0.444 | 20% |
| 8 | 14 | 0.1629 | 6.14 | 0.290 | 0.530 | 28% |
| 9 | 15 | 0.2357 | 4.24 | 0.215 | 0.572 | 33% |
| 10 | 16 | 0.2545 | 3.93 | 0.010 | 0.572 | 33% |
| … | | | | | | |
| 17 | 23 | 0.3328 | 3.01 | 0.169 | 0.641 | 41% |
| 20 | 26 | 0.3723 | 2.69 | 0.002 | 0.642 | 41% |

**Best single mode:** ProDy#4 (= repo#10), overlap **0.329**.
**Cumulative overlap:** 1 mode 0.100 (1%), 3 modes 0.267 (7%), **5 modes 0.441 (19%)**,
**10 modes 0.572 (33%)**, **20 modes 0.642 (41%)**.

**Cutoff robustness (15 Å):** best single 0.323 (same mode #4); cumulative 19% / 30% / 37%
at 5/10/20 modes. The story is not a cutoff artifact.

**Mode character (planar apo, softest 14 modes).** Helical/out-of-plane fraction of the
subunit-centroid motion: modes 1,2 = 0.77–0.79; modes 4,5 = **1.00**; modes 10–14 = 0.76–0.88;
radial "breathing" concentrates in modes 8,9 (0.69–0.91). **The very softest modes of the
planar ring are helical/out-of-plane** — i.e. the planar↔helical *direction* is intrinsically
soft, independent of any predictor.

## 4. Design vs native — did the covalent single chain keep the soft mode?

ANM on the cp233 single chain (1750 Cα, 5 copies of 342 core residues at resSeq starts
1/353/705/1057/1409; 40 linker Cα present but excluded from subunit centroids). C_start is
planar (centroid axial spread 0.05 Å).

- **Stiffness spectrum (fig. panel 2):** native and design low-mode eigenvalue spectra are
  nearly superimposable — the design was **not globally rigidified**.
- **Helical/out-of-plane soft modes are retained.** Design helical fractions: modes 1,2 =
  0.55–0.56; modes 7,8,10 = 0.84–0.88; modes 11,12 = 0.74–0.75. The softest HELICAL mode is
  ProDy#1 in *both* (native 0.77, design 0.55).
- **One honest nuance:** within the *very softest 5* modes, native = 5/5 helical whereas
  design = 2/5 (design modes 3,4 are in-plane, mode 5 is radial breathing). So the covalent
  linkers reshuffle the softest band slightly — the pure planar↔helical character first
  appears at design mode 1–2 and again strongly at 6–8, rather than saturating modes 1–5.
  The mode is **kept and still soft**, just not as exclusively dominant in the lowest band.
- **Ring-opening/breathing content** (radial collective coordinate captured by softest modes,
  independent ProDy): native apo **top20 = 0.24**, design C_start **top20 = 0.14** (top5:
  native 0.00, design 0.10). Comparable order — the design retains ring-opening character in
  its low modes; it is **not over-rigidified**.

## Cross-check reconciliation (adversarial verification)

Independent ProDy values vs the repo-tool numbers (`reproduce/overlap.py`,
`reproduce/anm.py`, `reproduce/enm_openclose.py`, re-run and confirmed this session). Note the
repo tools use *different* structures (Boltz native ring + `7JQQ.cif`; AF3 cp233), so agreement
is a genuine independent confirmation, not a tautology.

| quantity | repo tools | ProDy (this work) | verdict |
|---|---|---|---|
| ring apo↔7JQQ RMSD | 6.6 Å | 6.65 Å | ✓ agree |
| best chain assignment | (B,A,E,D,C) | (B,A,E,D,C) | ✓ agree |
| best single-mode overlap | 0.33 (mode 10) | 0.329 (ProDy#4 = repo#10) | ✓ agree |
| cumulative, 5 modes | 19% | 19% | ✓ agree |
| cumulative, 10 modes | 33% | 33% | ✓ agree |
| cumulative, 20 modes | 41% | 41% | ✓ agree |
| softest modes helical | 0.88–1.0 (modes 7,8,10,11) | 0.77–1.0 (ProDy#1,2,4,5) | ✓ agree |
| breathing, native top20 | 0.24 | 0.24 | ✓ agree |
| breathing, design top20 | 0.13 (AF3 cp233) | 0.14 (C_start cp233) | ✓ agree |

**No material disagreement.** Two independent codebases, two independent input sets, same answer.

---

## Honest read

- **Is P↔H a soft mode? Partly — and in the right way.** The planar ring's *softest* modes
  point in the planar↔helical (out-of-plane, staircase) direction (helical fraction 0.77–1.0
  in modes 1,2,4,5), so the transition direction is **geometrically accessible and low-energy**,
  built into the fold. But **no single soft mode dominates**: the best single ANM mode captures
  only ~11% of the P→H displacement (overlap 0.33), and it takes ~10–20 modes to reach 33–41%.
  The transition is **collective/multi-mode**, not one clean hinge.
- **Did the design keep it? Yes.** The cp233 single chain retains low-frequency helical/out-of-
  plane and ring-breathing modes with a stiffness spectrum matching native. The covalent linkers
  did **not** stiffen the P↔H motion out. (Minor caveat: the softest 5-mode band is slightly less
  purely-helical than native, 2/5 vs 5/5, but helical modes remain soft.)
- **What ENM does NOT show.** ENM reports which motions are *low-energy and available*, not their
  *driven energetics*, not the barrier, and not the mechanism (concerted vs hand-over-hand). This
  is a **necessary-condition** pass, not proof. A low overlap would have been the falsifying result;
  it came out accessible-but-collective.

### Implication for C3 (enhanced sampling)
Because the P→H transition is **collective (≈33% at 10 modes, no dominant single mode)** rather
than a single soft hinge, the driven-MD collective variable in C3 **should not be a single ANM
mode or a scalar planarity metric** — either would under-drive the true path. Use a
**difference-RMSD to the H (7JQQ) endpoint** (the full Δ vector, saved here) or a **string / path
CV** spanning P→H. The saved Δ and the 20 native ANM modes (`C2_enm.npz`) can seed that string or
serve as a low-dimensional TMD/SMD basis.

## Files
- `C2_REPORT.md` — this report
- `C2_overlap_figure.png` — per-mode overlap + cumulative curve; native/design stiffness spectrum; helical character
- `C2_enm.npz` — `Delta`, `Delta_hat`, native `eigvals`/`eigvecs`, cutoff-15 arrays, design `eigvals`/`eigvecs`, helical-fraction arrays, breathing arrays, chain order, `Xp`/`Xh_aligned`
- `c2_enm_analysis.py` — the analysis script (reproducible; ProDy 2.6.1, `/Users/longfu/miniforge3/bin/python3`)
