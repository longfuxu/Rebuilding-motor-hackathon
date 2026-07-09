# OpenMM MD physics-validation of the gp16 single-chain ring — methods

**Goal.** A predictor-independent, physics-based check of the single-chain design
`cp233_int15_inter10`: does the designed ring stay closed/stable under molecular
dynamics like a genuine closed apo ring, or does it open / collapse / dissociate?
This is orthogonal evidence to the Boltz/OF3/AF3 fold-in-context M1–M4 gate — MD
knows nothing about the structure predictors.

## Three systems (identical protocol, so results are directly comparable)

| id | system | source | ligands | starting state |
|----|--------|--------|---------|----------------|
| **A** | apo closed ring | Boltz native-ring `outputs/structures/cycle3/native_ring__boltz__s1__native_ring_model_0.pdb` | none | closed, C5-symmetric, planar |
| **B** | 7JQQ ATP-bound helical | cryo-EM `data/raw/7JQQ.cif`, protein chains A–E only | **stripped** (AGS/ATPγS, Mg, pRNA, DNA removed) | asymmetric 3-ATP intermediate (tight ATP-side fingers + open apo seam) |
| **C** | design single chain | `outputs/structures/af3_sweep/cp233/inter10/fold_2026_07_08_cp233_int15_inter10_model_0.cif` | none | closed, sequential-consistent (1750-aa, 5 covalent copies) |

**Why A is the Boltz apo ring, not 7JQQ-stripped.** 7JQQ's protein is frozen in the
ATP-driven helical conformation; stripping ligands does not relax it to the apo basin
on an MD timescale. The genuine *closed, symmetric* reference the design must resemble
is the apo predicted ring (A). B (7JQQ, ligands stripped) is retained as the
asymmetric ATP-bound-conformation reference — its ligand-free MD tests whether the
helical state is intrinsically metastable, and gives a "what an open/asymmetric ring
looks like" contrast.

**ATPγS/Mg note.** For tractability and robustness on a headless GPU VM, B is run
protein-only (no ATPγS parameterization). B therefore probes the metastability of the
ATP-bound *backbone conformation*, not ATP-driven dynamics. Explicit ATPγS+Mg
(GAFF/AM1-BCC) is a documented Phase-2 escalation, not run here.

## t=0 baselines (from `reproduce/score_m2.py`, reproduced exactly by the MD analyzer)

| system | M2 engaged | R146→WalkerA (Å) | radius (Å) | radius_CV | planarity_rms (Å) |
|--------|-----------|------------------|-----------|-----------|-------------------|
| A apo  | 5/5 | 6.5–6.7 uniform | 26.9 | 0.00 | 0.1 (planar) |
| B 7JQQ | 3/5 | 2.8/3.2 fingers + 15.5/25.5 open seam | 29.1 | 0.02 | 1.8 (non-planar) |
| C design | 5/5 | 7.3–7.6 uniform | 25.8 | 0.01 | 0.1 (planar) |

## Pipeline (`md_driver.py`)

1. **Prep** — PDBFixer: fill missing sidechain atoms (missing *loops* NOT modelled:
   `missingResidues={}`), remove heterogens, add hydrogens at pH 7.0. Author residue
   numbering preserved through prep and PDB writes (`keepIds=True`) so gp16 numbering
   (R146=146, Walker-A 24–31) survives — required for the analyzer.
2. **Force field** — `amber14-all.xml`. Implicit solvent **GBSA-OBC2**
   (`implicit/obc2.xml`), soluteDielectric 1.0 / solventDielectric 78.5,
   **CutoffNonPeriodic 2.0 nm** (see benchmark below).
3. **Integrator** — LangevinMiddle, 300 K, 1 ps⁻¹ friction; HMR (H=1.5 amu) + HBonds
   constraints → 4 fs timestep.
4. **Protocol** — energy minimize (5000 it) → 100–200 ps NVT equilibration (300 K) →
   NVT production per system. Trajectory (DCD, 20 ps/frame) + energy (StateDataReporter
   CSV) written incrementally; checkpoint every 500 ps (restartable). **Production
   length: A = 3 ns, B/C = 2 ns.** Colab sessions on this account self-terminated at
   ~20–30 min (a documented Colab flakiness — see below), so B/C were run at 2 ns
   (~22 min/session) to complete inside one session; A's own curves show the RMSD /
   radius_CV plateau is reached by ~1–1.5 ns, so the shorter B/C runs lose no signal.
   All three use identical force field, solvent, integrator, and readouts.

   *Session robustness:* runs are driven by a disconnect-tolerant orchestrator
   (`run_tolerant.py`) that owns the session, retries transient CLI websocket drops
   without relaunching, downloads each system's results the instant it finishes, and
   only recreates the session (reinstall+reupload+relaunch) after `colab sessions`
   confirms a true teardown. The A100 self-terminated ~twice during the campaign; the
   orchestrator recovered each time and no completed result was lost.
5. **Platform** — OpenMM 8.5.2 on Colab-Pro **A100-40GB**. (PyPI wheel exposed no CUDA
   plugin on this VM; **OpenCL** platform used — GPU-accelerated on the A100.)

## Compute budget (Phase 0 benchmark, measured on the 27,770-atom design, A100/OpenCL)

| build | ns/day | note |
|-------|--------|------|
| implicit GBSA **NoCutoff** (O(N²)), 4 fs HMR | 19.8 | too slow for a screen |
| implicit GBSA NoCutoff, 2 fs | 9.9 | — |
| implicit GBSA **2.0 nm cutoff**, 4 fs HMR | **145** | **7.3× faster — used for the screen** |
| explicit TIP3P+PME (~250k atoms), 4 fs | *failed* | NaN blow-up (shallow bench minimize); would need staged equilibration and, at ~10–20 ns/day on OpenCL, ~1 day/system — **deferred**, not run |

Screen cost: 3 systems × (5 ns prod + 200 ps equil) ≈ 15.6 ns ≈ **~2.6 h A100** wall.
Explicit all-atom (the gold standard) is left as a documented Phase-2 escalation: it is
~10× the per-ns cost and needs CUDA-enabled OpenMM (conda-forge) to be practical.

## Readouts over time (`analyze.py`, consistent with `score_m2.py`)

- **Cα-RMSD(t)** vs the production start (backbone drift).
- **Per-residue RMSF** (flexibility).
- **Ring closure** — for each *sequential* interface, R146 guanidinium (NE/CZ/NH1/NH2)
  → neighbour Walker-A (res 24–31) minimum heavy-atom distance. Adjacency is **locked
  at t=0** so each trace follows the same physical seam. engaged < 8 Å.
- **Radial symmetry** — `radius_CV` (std/mean of subunit-centroid radii; 0 = perfect ring).
- **Interface contacts** — inter-subunit residue pairs (closest-heavy < 4.5 Å) per seam.

The analyzer's t=0 interface distances reproduce `score_m2.py` exactly for all three
systems (A 6.57/6.62/6.53/6.64/6.69; B 15.51/3.17/2.81/5.79/25.46; C 7.32/7.31/7.44/7.59/7.32).

## Decision criterion

Design C is physically plausible if, under the **same** protocol, its stability and
ring-closure track the genuine apo ring A (RMSD plateaus, radius_CV stays low, ≥4/5
interfaces stay engaged). If C rapidly opens / collapses / dissociates while A stays
closed → the design is not physically reasonable. B contextualizes what an asymmetric
ATP-conformation ring does under the same force field.
