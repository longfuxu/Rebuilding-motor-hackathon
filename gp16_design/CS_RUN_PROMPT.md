# Claude Science autonomous run — gp16 E119Q + single-chain closure (≈4–5 h, continuous)

> Paste everything below the line into Claude Science. It is self-contained and designed to run
> unattended for 4–5 hours, working through a priority queue without stopping to ask.

---

## RUN MODE (read first)
Run **autonomously for ~4–5 hours**. Work the PRIORITY QUEUE in order; when a task finishes, immediately start the next — **do not stop to ask for confirmation or approval**. Save results **incrementally** after every fold (never hold everything to the end). On a single fold failure: retry once, then log it and move on. If you finish the queue early, **deepen** it (more seeds, more linker lengths) rather than idling. Produce the SUMMARY (last section) at the end. **All compute is FREE — NVIDIA NIM (Boltz-2 + OpenFold3) + ColabFold MSA. Do NOT use Modal.**

## CONTEXT (self-contained)
Project: design a **single-chain, position-addressable φ29 gp16 packaging-motor ring** (one genetically-addressable seat, like single-chain ClpX). Public inputs: PDB **7JQQ**, UniProt **P11014** (gp16, 332 aa). Workspace: this sandbox `~/.claude-science/orgs/34666d1b-…/workspaces/9c944ce5-…/`. Fold tool: the user's `cofactor.fold` (stdlib ColabFold→NIM wrapper), key via `NVKEY_DIRECT`. Cached tiled MSAs already exist: `handoff/tiled_msa/{dimer,trimer,tetramer,pentamer}.a3m` and `handoff/core.a3m` (327-aa core, 1168 records) — **reuse them; only re-tile when the linker length changes** (re-gap the linker columns; the per-copy core homologs are unchanged).

State so far (cross-checked, NOT validated; apo):
- **B1 (covalent dimer + 3 WT, "dimer-in-ring") is the robust position-addressable construct** — Boltz-2 15/15 and Chai-1 47/50 reform the catalytic trans-R146 interface like native.
- **The fully single-chain pentamer is NOT predictor-robust:** Boltz-2 closes 5/5 (4.9–5.5 Å) but OpenFold3 closes only 2/5 — OF3 builds a **compact but topologically scrambled** ring (each open R146 buried against the WRONG partner) despite *higher* global pTM (0.58 vs 0.48).
- **E119Q constructs are prepped but NOT folded** (`fold_inputs/e119q/`). E119 = Walker-B catalytic glutamate; E119Q = Glu→Gln (a conservative, catalytically-dead mutation).

## THE METRIC CONTRACT — score by THIS, never by global pTM
Steering by global pTM/pLDDT is the known trap (it rewarded OF3's scrambled ring). Judge every fold by **M2**, exactly as validated in `gp16_design/reproduce/score_m2.py` (it reproduces the cycle-3 and ladder reference numbers to 0.01 Å):

> **M2 = per-interface trans-R146 arginine-finger engagement.** For each subunit's R146 **guanidinium** (atoms NE, CZ, NH1, NH2), the minimum heavy-atom distance to a **neighbour's Walker-A (res 24–31)**. "engaged" = **< 8 Å** (native 7JQQ ≈ 3.2; predicted apo ring ≈ 6.6). Pass = **≥ 4/5** interfaces engaged.
> For covalent single-chain constructs, score each copy against its **DESIGNED sequential ring partner** (copy k → copy k+1, cyclic) — **not** the nearest subunit — so a compact-but-scrambled ring is correctly scored as failed. Also record, per open interface, the *nearest* subunit (if it's the wrong partner, flag "scrambled"). For a res4–330 core, R146 = copy-position 143, Walker-A = 21–28 (use `--copy_start_res 4`).

Also report per-interface pLDDT (M4) and ring-closure geometry (M1: are the 5 protomer centroids on a common circle / is the ring closed vs an open spiral). **Report cross-predictor divergence loudly — it is the confidence signal, not a failure to hide.**

## COMPUTE ROUTING
Per fold: ColabFold core MSA (cached) → block-diagonal tile per copy → Boltz-2 NIM (+ OpenFold3 NIM cross-check). Reuse `handoff/core.a3m`. Multi-chain constructs (native ring, B1): give each WT chain the core MSA; give the covalent chain its tiled MSA. `proteinmpnn-nim` is available for Task C. RFdiffusion is NOT on NIM — if Task C needs backbone design, note it as a follow-up, don't block.

---

## PRIORITY QUEUE

### Task A — E119Q assembly-tolerance on B1 (HIGHEST; the "why we built it" construct)
Fold, cross-checked by Boltz-2 + OpenFold3 (≥3 Boltz seeds each):
1. **B1_L40 dimer-in-ring with E119Q on the tethered module + 3 WT** (`fold_inputs/e119q/`).
2. **Untethered control:** 5 separate WT chains, one carrying E119Q (dead seat present but position not covalently defined).
3. **WT B1_L40** baseline (re-fold for a matched comparison).
Score all by M2. **Question answered at the fold level:** does the dead-seat mutation stay **assembly-compatible** — does the ring still close (M2 ≥ 4/5) with E119Q in place, and is any geometric perturbation **local to the known E119Q seat** vs delocalized across the ring? Compare the tethered (defined-position) construct to the untethered control. **Be honest:** E119Q is a conservative mutation, so the apo fold may be ~identical to WT — the fold-level result is *assembly tolerance + a valid defined-dead-seat construct*; the functional "defect stays local" claim is future wet-lab. Save structures + a per-interface M2 table.

### Task B — single-chain pentamer robustness sweep (attack the Boltz↔OF3 divergence)
Goal: find a construct where **both** Boltz-2 AND OpenFold3 close ≥ 4/5 by M2 (designed-sequential).
1. **Seed depth:** re-fold the current pentamer (linker (GGGGS)×8) with **≥5 Boltz-2 seeds** and **≥3 OF3 seeds**; is Boltz consistently 5/5? does OF3 ever un-scramble? Report the distribution.
2. **Linker-length sweep:** rebuild the pentamer with linkers **(GGGGS)×6, ×8, ×10, ×12** (30/40/50/60 aa). Re-tile the core MSA per length. Fold each with Boltz-2 (3 seeds) + OF3. Score M2. Look for a length where OF3 also closes.
3. **(optional) register:** if a linker length gets Boltz robust but OF3 still scrambles, try shifting the seam (which junction is the non-covalent one) by one position.
Save a table: construct × predictor × seed → M2 engaged/5 + scrambled flags. **If no linker/register makes it cross-predictor-robust → that is the documented trigger for generative escalation (Task C); state it explicitly, it's a result not a failure.**

### Task C — (conditional on B failing) ProteinMPNN self-repair
Only if Task B yields no cross-predictor-robust pentamer. Take the best **near-closed** pentamer backbone; use `proteinmpnn-nim` to redesign the **linker + interface-adjacent positions**, with **catalytic/interface residues FIXED** (Walker A res 24–31, Walker B res 115–119, R146 arginine finger, ATP pocket). Re-fold the redesigned sequence (Boltz-2 + OF3), score M2. This is the framework's self-repair path. If it needs backbone (RFdiffusion) rather than sequence design, log it as the next-session follow-up.

### Task D — housekeeping (interleave whenever a fold is running)
- Copy the MSA structures (`outputs/structures/ladder/*__boltz2nim_msa.cif`, `*__openfold3_msa.cif`) and the score CSVs so they are captured with provenance (they belong in the repo for data availability).
- Keep `folding_results.md` updated incrementally with each new result + job/request IDs.

---

## GUARDRAILS (hold throughout)
- "Cross-checked by ≥2 predictors," **never** "validated." Folds are apo (no ATP/pRNA/DNA). Constructs address a **geometric seat**, not a guaranteed fixed function.
- **Steer by M2 (designed-sequential), not global pTM.** Do not pick a "winner" on pTM.
- Report predictor divergence explicitly. Do not smooth it over.
- Do not over-claim E119Q locality from an apo fold (see Task A).
- Save incrementally; assume the session could end at any time.

## SUMMARY (produce at the end, write to `outputs/RUN_SUMMARY.md`)
1. Table: every construct folded × predictor × M2 engaged/5 (designed-sequential) + scrambled flags + M1 closure call.
2. **Task A:** does E119Q-B1 still close (assembly-tolerant)? tethered vs untethered.
3. **Task B:** best linker/register; is any pentamer now cross-predictor-robust? If not, state the generative-escalation trigger is met.
4. **Task C** (if run): did ProteinMPNN redesign close it under both predictors?
5. One-paragraph "what is robust / what is not / recommended next," in the honest cross-checked-not-validated voice.
6. List all saved structures + job/request IDs for provenance.
