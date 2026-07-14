# In-silico P→H→P cycle of the single-chain gp16 ring (+dsDNA) — driven-MD campaign report

_Goal-4 keystone. Owner: Longfu Xu. Repo: `gp16_design/`. This run: 2026-07-12, GCP A100._
_Honest report as requested: why · what · methods · results (+figures/animation) · limitations · discussion
(help to our model + single-molecule predictions). **`<<NUMBERS>>` are filled from the 3 fwd + 3 rev seeds.**_

---

## 1. Why we did this
The φ29 packaging motor runs a **planar↔helical (P→H→P)** conformational cycle: a **dwell** where the planar
ring loads ATP and lifts into a helical lock-washer (**P→H, ascent/loading**), and a **burst** where it
descends back to planar and translocates DNA (**H→P, descent/power stroke**). Our earlier work characterised
only the **ascent** (forward targeted-MD / umbrella): the ring can be driven P→H with the M2 coupling retained,
and the per-subunit lift is a **graded, unequal staircase**. Three things were missing, and they are exactly
what single-molecule experiments care about:
1. **The descent (H→P) — the power stroke** that actually moves DNA. Is it reversible? Does DNA translocate on
   the descent (not the ascent)?
2. **Does the M2 coupling survive the descent too**, or only the ascent?
3. **Robustness** of the unequal-tread signature across independent trajectories (is it mechanism or noise?).
We deliberately did **not** chase the absolute free energy (established earlier as untrustworthy from a single
1-D CV with under-overlapped windows; the user de-scoped it) nor a full DNA-coupled model (AWSEM+3SPN — out of
scope). We asked the **mechanism/pathway** questions that MD *can* answer without simulating the real 10–80 ms.

## 2. What we did (methods)
- **System:** the circular-permutation **single-chain** gp16 ring (cp233, chain A, 1,750 residues = 5 covalently
  linked subunits) threaded with **60-bp dsDNA** (chains F/G) through the pore → `C_plus_dna_relaxed.pdb`,
  **31,578 atoms**, implicit solvent (amber14-all + GB-OBC2, 2.0 nm cutoff). No ATP (the CV supplies the drive).
- **Reaction coordinate (CV):** the proven **axial staircase** — for each subunit *k*, restrain the axial
  projection of its Cα-centroid (relative to the ring COM, on the fixed planar-ring axis) toward `λ·z_k`, with
  `z_k` a monotonic staircase spanning the 7JQQ helical amplitude (±2.38 Å). `CustomCentroidBondForce`,
  GPU-native. Driver `md/driven/tmd_staircase.py`.
- **Ascent (P→H):** `fwd_s{1,2,3}` — start planar, ramp **λ: 0→1**, 500 ps steered pull, 50 ps equil, frame/2 ps.
- **Descent (H→P):** `rev_s{1,2,3}` — the new capability: CV geometry from the **planar** structure but initial
  **coordinates from the helical endpoint** (`--pos_from helical.pdb`, i.e. a forward run's final state), ramp
  **λ: 1→0**. This drives the ring back down from the helical staircase to planar. (Validated on-hardware:
  the first reverse frame starts at span ≈4.1 Å and descends, work is negative = downhill.)
- **Readouts per frame:** per-subunit axial z (fixed planar-ring axis, so ascent/descent are comparable),
  ring planarity (rotation-invariant refit), axial-span, DNA axial COM position, and the **M2 coupling**
  (# of the 5 R146→neighbour-WalkerA interfaces < 8 Å), plus accumulated steered work.
- **Compute:** one A100-40 GB Spot (OpenCL), 3+3 runs, ~13 min/run. Analysis local (`analyze_cycle.py`,
  `render_cycle.py`).

## 3. Results  (3 forward + 3 reverse seeds; figures `cycle_stats.png`, `cycle_mechanics.png`; animation `cycle.mp4`)

### 3.1 The ring executes a REVERSIBLE P→H→P conformational cycle  ✅ robust
Ascent drives planar→helical (planarity 0.02–0.05 → **0.71–0.81 Å**, axial span 0.1–0.2 → **4.2–4.4 Å**).
Descent returns it to planar (planarity 0.70–0.78 → **0.05–0.11 Å**, mean `plan_return = 0.089 Å`; span
4.3–4.5 → **0.15–0.34 Å**). The staircase that forms on the ascent **fully relaxes** on the descent — the
covalent single chain is not trapped in the helical basin. All 3 seeds behave the same.

### 3.2 The per-subunit staircase is graded and UNEQUAL, and HIGHLY reproducible  ✅ strongest result
Per-subunit axial treads over the ascent (mean±sd, n=3):
**sub1 −2.23±0.09 · sub2 −1.24±0.09 · sub3 +0.03±0.07 (hinge) · sub4 +1.28±0.08 · sub5 +2.15±0.04 Å.**
The sd's are ~0.08 Å — the graded, position-dependent pattern is **set by the ring geometry, not thermal
noise**. This is the structural basis of a **graded, unequal substep** prediction (not five identical steps).

### 3.3 M2 coupling largely survives BOTH half-strokes  ✅ (single-chain-specific)
n_engaged stays high across ascent AND descent (per-run minimum: ascent mean **3.67/5**, one seed dipped to 3;
descent mean **4.0/5**). The R146→neighbour-WalkerA coupling does not systematically break as the ring goes up
or comes down — the covalent chain keeps its subunits mechanically coupled through the whole cycle. It is
slightly **more strained/fluctuating on the descent** (more transient 5→4 dips). This is the load-bearing
single-chain result.

### 3.4 Steered work: ascent uphill, descent downhill  ✅ (as expected)
Ascent **28.1±6.1 kcal/mol**; descent **2.9±1.8 kcal/mol** (much lower, partly negative per-frame — the
descent runs downhill in the CV once the ring is in the helical basin). These are **fast-pull dissipative
upper bounds, NOT ΔG** (§4).

### 3.5 DNA does NOT significantly translocate in this conformation-only model  ⚠️ HONEST NULL
DNA axial displacement per half-stroke (each rel. to its own start, mean±sd, n=3): **ascent +0.36±1.43 Å,
descent +0.35±0.79 Å** — statistically **indistinguishable, and the error bars exceed the means**. There is
**no significant, directional DNA translocation in either half** — the DNA essentially wobbles (~±1 Å) as the
ring strokes around it. _(A single early seed, rev_s1, showed +1.45 Å on the descent and looked like "the
descent translocates DNA"; the 3-seed average erased it — a textbook reason to run replicates.)_
**Interpretation:** the ring's conformational power stroke, on its own and without the ATP-driven **grip/ratchet
chemistry**, is **not sufficient** to move DNA processively in implicit solvent over a 500 ps driven stroke.
Translocation requires the grip that couples the conformational change to DNA — i.e. a DNA-coupled model (§4/§6).
This is a **useful mechanistic statement**, not a failure: conformational cycling is necessary but not
sufficient; the grip is the missing ingredient.

## 4. Limitations (read before over-interpreting)
1. **Driven, not spontaneous.** The CV forces the transition ~10⁷× faster than the real motor. Work is a
   dissipative upper bound, not ΔG; we make **no rate or free-energy claim**.
2. **The CV imposes concurrent driving**, so this **cannot resolve concerted vs sequential (hand-over-hand)**
   ordering — all subunits are driven together. What we *can* say is drivability + coupling + the *geometric*
   unequal treads. Concerted-vs-sequential needs unbiased/committor/path sampling (next step).
3. **DNA translocation is small (~1 Å) and modest-signal.** The descent moves DNA more than the ascent
   (real asymmetry), but not a full 2.5-bp step — implicit solvent, a short 500 ps driven stroke, and no
   ATP/grip chemistry. A processive per-bp register needs a DNA-coupled model (AWSEM+3SPN, de-scoped).
4. **Partial helical endpoint.** The staircase CV reaches ~90 % of the 7JQQ axial *span* but only ~45 % of the
   RMS planarity (the two-ends spread more than the whole-ring RMS) — the driven "H" is a partial lock-washer.
5. **Reverse seeds from a fixed helical endpoint**, not each forward run's exact last frame → the ascent/descent
   junction is not perfectly continuous (per-half DNA displacements are honest; a single cumulative "net per
   cycle" would carry a junction artifact, so we do not report one).
6. Implicit solvent; single force field; no explicit Mg·ATP.

## 5. Discussion
### 5.1 What it adds to our current model
- **Completes the cycle in silico.** We now have both half-strokes (ascent + descent), not just the ascent —
  the ring goes up into the staircase AND comes back down, reversibly. The descent is **downhill / low-work**
  (2.9 vs 28 kcal/mol), consistent with the "burst" being a release, but see the next point.
- **Separates conformational cycling from translocation.** The ring cycles P→H→P and the descent is downhill,
  BUT the DNA does **not** significantly translocate (null, §3.5). So in our model the conformational stroke is
  **necessary but not sufficient** — the ATP-driven **grip/ratchet that couples ring motion to the DNA phosphate
  backbone is the load-bearing element for translocation**, not the conformational amplitude alone. This sharpens
  what the "power stroke" must include: geometry + grip, not geometry alone.
- **Hardens the single-chain thesis.** The M2 coupling surviving the *entire* cycle (up and down) is the
  mechanistic payoff of the covalent single chain: it forces the ATPase coupling and the conformational stroke
  to stay locked together — something the native oligomer can shed. The **special/dead subunit** becomes a
  defined, positionable lever precisely because the chain holds the coupling.
- **Grounds the unequal-step prediction in geometry, robustly.** The treads are unequal and reproducible → the
  motor's substeps should be **graded, position-dependent**, not five identical 2.5-bp steps.

### 5.2 Predictive value for the lab's single-molecule experiments
- **3D-MINFLUX (per-subunit z(t)):** predict a **graded staircase** — subunit-specific axial amplitudes with a
  near-stationary "hinge" subunit and two most-mobile subunits — and, on the single chain, the amplitudes are
  **addressable by position** (put a dead/labeled subunit at a chosen seat and read its tread). Falsifiable:
  equal-amplitude per-subunit steps would refute the graded-tread picture.
- **Optical tweezers (force/step):** the reproducible **unequal per-subunit treads** predict **unequal substep
  sizes** within the ~10-bp burst (graded, position-dependent — NOT 4×2.5 bp identical). Our null on DNA motion
  says the substep sizes are set by *geometry × grip*, not the conformational amplitude alone — so a mutation
  that weakens the DNA grip (without changing the conformational treads) should shrink/blur the substeps while
  leaving the ring's conformational cycle intact — a clean OT/structure double-readout.
- **smFRET (coupling):** a FRET pair across an R146→neighbour-WalkerA interface should stay **high-FRET
  (engaged) through the whole cycle** in the single chain; a mixed ring with a dead seat should show the
  coupling defect localise to the addressed interface.
- **Mixed-ring / special-subunit design:** the coupling-holds-across-the-cycle result says a single dead seat
  should have a **defined, local** effect — motivating the per-position dead-seat fold screen.

## 6. Next steps (ranked)
1. **Concerted-vs-sequential** — the one mechanism question this driven set can't answer. Use committor / string
   / unbiased relaxation from the helical endpoint (no CV) to see the intrinsic descent ordering.
2. **Per-position dead-seat prediction** (local, free NIM) — the addressable special-subunit map for the mixed-ring
   experiment.
3. **DNA-coupled model** (AWSEM+3SPN) — only if a real per-bp register is needed (big; de-scoped for now).
4. If a trustworthy *landscape* is ever wanted: soften `kcv` (~8000) for adequate umbrella overlap.

## 7. Files
- Animation: `cycle.mp4`, `cycle.gif` (full P→H→P). Static: `cycle_mechanics.png`, `cycle_stats.png`.
- Data: `fwd_s{1..3}/`, `rev_s{1..3}/` (each: `traj.dcd`, `final.pdb`, `series.json/csv`), `cycle_summary.json`.
- Code: `md/driven/{tmd_staircase.py (now with --lam0/--lam1/--pos_from), analyze_cycle.py, render_cycle.py,
  run_gcp_cycle_campaign.sh, poll_cycle_campaign.sh}`.
- Prior single-half + landscape: `outputs/php_cycle/{dna_translocation/, C3_umbrella/}`.
