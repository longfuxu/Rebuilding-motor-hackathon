# Goal-4 in-silico P→H→P session log (2026-07-12) — attempts, methods, params, results, honest reliability

Owner: longfuxu@berkeley.edu. Worktree: `.../export-gp16-results/gp16_design`. Compute: GCP A100-40GB Spot
(OpenCL), local Mac (analysis/NIM). **This log is written to be reproducible and HONEST about what each result
can and cannot support.** Reliability tags: 🟢 solid · 🟡 suggestive/needs-more · 🔴 illustrative-model-only.

---

## 0. The goal
Computationally probe the φ29 gp16 **planar↔helical (P→H→P)** packaging cycle on the **single-chain** design
(cp233), and turn it into falsifiable predictions for the lab's single-molecule experiments (3D-MINFLUX, optical
tweezers, smFRET). We do NOT try to get absolute rates or ΔG (out of reach / de-scoped); we ask **mechanism /
pathway / feasibility** questions.

## 1. What we tried, in order (each: why · method+params · result · reliability)

### 1a. Umbrella sampling → P→H free-energy landscape G(ξ)  🟡→🔴 (absolute ΔG unreliable)
- **Why:** get the shape of the P→H landscape + coupling along the path.
- **Method:** `md/driven/tmd_umbrella.py`, staircase axial CV, implicit GBSA-OBC2, 27,770-atom single chain
  (`C_start.pdb`). `--nwin 24 --equil_ps 300 --sample_ps 1200 --kcv 30000 --passes both`. WHAM.
- **Result:** forward 24 windows done; coupling held (n_eng 5→4.79); graded unequal per-subunit motion.
  **Absolute ΔG UNRELIABLE:** WHAM gave +81 kcal/mol but that exceeds the 35–53 fast-pull work bound; root
  cause = `kcv=30000` too stiff → window overlap ratio 1.07 (>1 = insufficient) → WHAM under-determined.
  **Reverse STOPPED at the gate** (same kcv → same overlap → can't fix). `state_fwd_end.xml` saved (resumable).
  Output: `outputs/php_cycle/C3_umbrella/`. **Fix for a trustworthy landscape: `kcv≈8000` or 2× windows.**

### 1b. Driven-MD P→H→P cycle campaign (3 fwd + 3 rev seeds)  🟡
- **Why:** characterise both half-strokes with replicates (ascent P→H + descent H→P), ring+dsDNA.
- **Method:** `run_gcp_cycle_campaign.sh` → `tmd_staircase.py` on `C_plus_dna_relaxed.pdb` (31,578 atoms).
  fwd = `--lam0 0 --lam1 1`; rev = `--pos_from helical_endpoint.pdb --lam0 1 --lam1 0`. `--pull_ps 500
  --equil_ps 50 --kcv 30000`. Analysis `analyze_cycle.py` (fixed planar-ring axis).
- **Result:** reversible (descent returns ring to planar, planarity→0.09); **graded unequal treads highly
  reproducible** −2.23/−1.24/+0.03/+1.28/+2.15 Å ± ~0.08 (sub3 = hinge); coupling held both halves; work asc
  28±6 vs desc 2.9±1.8 kcal/mol. **DNA translocation NULL:** asc +0.36±1.43 vs desc +0.35±0.79 Å (error>mean).
  (A single seed rev_s1 gave +1.45 Å and looked positive → erased by the 3-seed average → run replicates.)
  Output: `outputs/php_cycle/cycle_campaign/` (`PHP_CYCLE_REPORT.md`, `cycle.mp4`, `cycle_stats.png`).

### 1c. Ring+DNA animation (single P→H drive, 200 ps)  🟢 (visual, honest)
- **Why:** the "motor gripping DNA" visual + mechanistic-observable dashboard.
- **Method:** `tmd_staircase.py` on `C_plus_dna_relaxed.pdb`, `--pull_ps 200 --report_ps 2`; render
  `render_dna_anim.py` (side-view 3D + λ vs planarity/span/n_engaged panels).
- **Result:** span 0.2→4.27 Å, n_eng 5→5 with DNA threaded, DNA nudged ~0.5 Å (not processive).
  Output: `outputs/php_cycle/dna_translocation/`.

### 1d. Continuous slow stroke (8 ns triangle λ 0→1→0)  🟡 (n=1)
- **Why:** a CONTINUOUS cycle (no fwd/rev seed discontinuity) + test "is the DNA-null just too-fast time?"
- **Method:** `tmd_staircase.py --triangle 1 --pull_ps 8000 --report_ps 8` (16× slower than the campaign).
- **Result:** DNA net **+3.54 Å** (vs ~0 fast) — moved on the descent; but single noisy trajectory (range
  −2.5..+4.4 Å). ⇒ **time/rate matters partly, but not enough alone.** Output: `dna_translocation/cycle_slow/`.

### 1e. Per-position dead-seat scan (A)  🟡 (static predictor)
- **Why:** does a single dead seat's effect depend on WHICH position (is a seat "special")?
- **Method:** `pipelines/php_cycle/deadseat_per_position.py` — WT + R146A at each of the 5 seats (chain res
  255/607/959/1311/1663), Boltz-2 NIM, 3 reps, handedness-robust M2 (`reproduce/score_m2.py`).
- **Result:** WT 5/5; R146A at EVERY position = 4/5, 3/3 reps, **spread = 0 → position-INDEPENDENT**.
  ⇒ no static-structure "special" seat; the special subunit (if real) is **dynamic/timing**.
  Output: `outputs/php_cycle/special_subunit/PER_POSITION_RESULTS.md`.

### 1f. Mechanochemical ratchet (C): grip + ATP-clock, concerted vs sequential  🔴 (imposed model)
- **Why:** add a **grip** (the missing translocation ingredient) + an ATP-firing **clock**; see if DNA is driven,
  and how concerted vs sequential differ.
- **Method:** `md/driven/tmd_ratchet.py` — per-subunit power stroke on `--mode concerted|sequential` + a
  catch-and-carry grip (spring couples dsDNA to a subunit during its down-stroke, releases on reset).
  `--ncycles 4 --cycle_ps 300 --step_A 3.4 --kgrip 12000`. Compare `compare_ratchet.py`.
- **Result:** **sequential translocates DNA −16 Å (~1.2 bp/cycle); concerted −1.6 Å**; both keep coupling 5/5.
  Sequential per-subunit z(t) = travelling wave; concerted = in-phase. Output: `outputs/php_cycle/ratchet/`
  (`ratchet_compare.png`, `RATCHET_RESULTS.md`). **NB: grip + clock are IMPOSED — this is an illustrative model
  of the LOGIC (continuous grip ⇒ translocation), NOT a measurement or proof of φ29's mechanism.**

## 2. ⚠️ Model-naming honesty (important)
- **Helical↔planar (spiral/lock-washer) model** = the field's φ29 hypothesis = the WHOLE ring's conformational
  cycle. That is what 1a/1b/1c/1d modelled (a collective P↔H drive). Finding: the conformational cycle **alone**
  (no grip) does **not** translocate DNA.
- **Ratchet (1f) is a HAND-OVER-HAND FIRING model**, NOT the helical↔planar model: subunits rest planar and bob
  on a per-subunit clock; it does not form a helical staircase. Its "sequential" is AAA+/rotary firing, a
  *different axis* from helical↔planar. **We have NOT yet run the integrated model** (helical↔planar WITH a grip
  active during the H→P descent) — that is the honest next step and the closest to the field's actual mechanism.

## 3. What ties together (the story, as hypotheses)
- Conformational cycling is **necessary but not sufficient** for translocation → the **grip/ratchet is the
  load-bearing element** (1b/1d null vs 1f positive).
- Efficient translocation needs **continuous grip** → favours sequential-like timing (1f).
- The special subunit is **dynamic/timing**, not static-structure (1e) → candidate **phase-setter** for the
  sequential wave (1e+1f).
- Single-molecule discriminator: per-subunit z(t) **travelling wave (sequential)** vs **in-phase (concerted)**,
  and DNA **substeps vs one-step** — MINFLUX + OT test it.

## 4. Reproducibility — exact commands
```bash
# umbrella (needs A100): bash md/driven/run_gcp_umbrella.sh ; analyze: python md/driven/analyze_umbrella.py
# cycle campaign:        bash md/driven/run_gcp_cycle_campaign.sh ; python md/driven/analyze_cycle.py
# ring+DNA anim:         TRIANGLE=0 OUTNAME=dna_anim bash md/driven/run_gcp_dna_anim.sh ; python md/driven/render_dna_anim.py --dcd ... --top ...
# continuous slow:       TRIANGLE=1 PULL_PS=8000 OUTNAME=cycle_slow bash md/driven/run_gcp_dna_anim.sh
# dead-seat (local NIM):  python pipelines/php_cycle/deadseat_per_position.py
# ratchet both modes:    NCYCLES=4 CYCLE_PS=300 KGRIP=12000 bash md/driven/run_gcp_ratchet.sh ; python md/driven/compare_ratchet.py
# every GCP run: poll_*.sh downloads + tears down the VM (cap=1 guard). Report credit, verify 0 instances.
```
Trajectories (`*.dcd`, ~1 GB) are gitignored — regenerate with the commands above; the small analysis outputs
(png/json/csv/md + fold CIFs) are committed.

## 5. Honest reliability / publishability (read §6 of the reports too)
- 🟢 **Solid / publishable now:** the single-chain design's predictor-robust ring closure + designed sequential
  coupling (prior work); the *feasibility/plausibility* framing of the driven-MD; the reproducible unequal-tread
  GEOMETRY.
- 🟡 **Suggestive, needs more before a mechanism claim:** the DNA time-dependence (1d, n=1), the dead-seat
  position-independence (1e, static predictor + geometric proxy).
- 🔴 **Illustrative model only (NOT a result to publish as mechanism):** the ratchet concerted-vs-sequential
  numbers (1f) — grip + clock are imposed; the −16 vs −1.6 is parameter-dependent; it demonstrates a known
  PRINCIPLE + a testable SIGNATURE, not a discovery about φ29.
- **To be community-accepted as MECHANISM:** need converged sampling (softer-kcv umbrella / metadynamics),
  explicit solvent + ATP·Mg, a DNA-coupled model with REAL grip residues (AWSEM+3SPN or all-atom), replicates,
  and consistency with the experimental observables. The present work is **hypothesis-generation + a framework +
  predictions**, with MINFLUX/OT as the arbiter — frame it that way.

## 6. Backlog / next
- Integrated **helical↔planar WITH grip** model (the field's actual mechanism; not yet run).
- Softer-kcv umbrella for a trustworthy G(ξ); replicates for 1d/1f; concerted grip variant that doesn't
  self-conflict (fairer concerted test); grip on real DNA-binding residues + explicit solvent.
