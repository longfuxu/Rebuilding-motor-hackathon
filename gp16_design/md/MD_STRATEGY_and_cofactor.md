# Simulation strategy — the scientific question, the hackathon/paper angle, and Cofactor integration

## 1. What scientific question does this simulation line answer?
The *design* work builds a **position-addressable gp16 ring** (put one defined subunit under control). The *biology*
that construct exists to probe is the **motor mechanism**: the ring opens **progressively with ATP occupancy**
(each ATP opens it a little; ADP closes it), through a continuum of nucleotide-indexed states. The core mechanistic
question — for gp16 and for AAA+ rings generally (cf. the ClpX central coupler, target #2) — is:

> **How does per-subunit nucleotide state control the ring's conformation, and does the firing ORDER
> (sequential vs concerted) determine the mechanical output?**

Static structure predictors cannot touch this (they give one closed ring regardless of ATP). Our cheap
occupancy elastic-network model already produces a concrete, testable answer: **the opening is progressive,
cooperative, and pattern-dependent** (3 ATP adjacent vs spread give different openings) — i.e. *sequential vs
concerted firing produce different ring states*. That is exactly the readout the addressable construct is built for.

## 2. Why this makes the paper high-impact (the arc)
It turns "we designed a protein" into **"we built the instrument AND the theory for a fundamental mechanistic question."**
- (i) a **transferable design framework** for position-addressable rings (methods contribution);
- (ii) demonstrated on gp16 (+ ClpX retrospective);
- (iii) a **mechanistic model** (occupancy-dependent, cooperative, order-dependent opening) that the construct is
  built to test — connecting *design → mechanism → (future) single-molecule*;
- (iv) a clean **"right tool for each layer"** methodological spine (below) — itself a contribution.
The construct makes a *prediction the model can be tested against*: place a catalytically-dead (E119Q) subunit at
position k → the occupancy model predicts how the opening pattern shifts → single-molecule optical tweezers
(Bustamante lab) measures it. Design, theory, and experiment close the loop.

## 3. The tool ladder — right tool for each question (a methods figure/table)
| Question | Tool | Cost | Status |
|---|---|---|---|
| Does a designed construct fold / reform the interface? | structure predictors (Boltz/OF3/AF3) + M1–M4 | free NIM | done |
| What are the ring's intrinsic soft motions? | elastic-network **NMA** | local, seconds | done (`reproduce/anm.py`) |
| Is the real opening a soft mode? | **mode–transition overlap** vs 7JQQ | local, seconds | done (`reproduce/overlap.py`) — enriched but not one mode |
| Does each ATP open the ring; does the pattern matter? | **occupancy-dependent ENM** (linear response) | local, minutes | done (`md/occupancy_enm.py`) — progressive + cooperative + pattern-dependent |
| The transition path / barrier / dynamics | **Cα dual-basin Gō MD** (apo+7JQQ basins) | GPU ~1 h | script ready (`md/gomodel_dualbasin_modal.py`, Modal) |
| Less-biased CG dynamics of the opening | **MARTINI CG-MD** | GPU hours | documented, later |
| Local stability of a designed construct | **all-atom short MD** (OpenMM) | GPU | for CP/B1 constructs |
| The real per-subunit mechanism | **single-molecule optical tweezers** + cryo-EM at controlled occupancy | wet lab | the ground truth the whole stack points to |

Static predictors for static folds; NMA for soft motions; occupancy-ENM for which-subunit-matters; Gō/CG-MD for
transition dynamics; single-molecule for the truth. **Choosing the tool to the question IS part of the contribution.**

## 4. Cofactor integration (design — Cofactor stays the owner's private instrument, not a hackathon deliverable)
The cheap, local analyses are ideal Cofactor "instruments" (stdlib+numpy+scipy, no GPU), mirroring how `fold.py`
is a local-MSA + NIM-backed tool. Proposed Cofactor tools:
- `cofactor.nma(structure)` → soft modes + character (from `reproduce/anm.py`).
- `cofactor.mode_overlap(state_a, state_b)` → how much of a conformational change is soft-mode (from `overlap.py`).
- `cofactor.occupancy_enm(apo, liganded, patterns)` → opening vs occupancy number/pattern (from `md/occupancy_enm.py`).
- `cofactor.transition_path(state_a, state_b)` → aANM/NEB path + barrier estimate (local; refine with Gō).
- `cofactor.go_md(state_a, state_b, ...)` → **compute-backed** tool that dispatches the dual-basin Gō MD to Modal
  GPU (`md/gomodel_dualbasin_modal.py`) — same "local-cheap default + cloud for the heavy step" pattern as `fold.py`.
This extends Cofactor from **structure prediction → conformational dynamics** — a natural, high-value capability
jump, and it makes the whole "tool ladder" reproducible from one instrument.

## 5. What is done vs next (the local-first, escalate-as-needed plan)
- Done, local, free: NMA (`anm.py`), mode–transition overlap (`overlap.py`), occupancy-ENM (`occupancy_enm.py`),
  NEB-lite (single-basin, monotone — motivates the dual-basin need).
- Ready to run on Modal top-up: **dual-basin Gō MD** (`gomodel_dualbasin_modal.py`, syntax-checked; tune eps/T/length
  on the first run) → the transition path/barrier and the effect of occupancy on the dynamics.
- **(c) LAST rung, documented:** MARTINI (or a multi-well Gō) CG-MD on H100 for a less-biased, physics-based view of
  the occupancy-driven opening — only if the dual-basin Gō leaves questions.
