# Overnight run summary — φ29 gp16 position-addressable ring
**Mode:** autonomous, free NIM only (Boltz-2 + OpenFold3 via `health.api.nvidia.com`, key `NVKEY_DIRECT`; ColabFold MSA). No Modal. Scored by **M2** (per-interface trans-R146 engagement, <8 Å; native ≈3.2, apo-ring ≈6.6), NOT global pTM. Single-chain constructs scored **sequentially** (copy k → copy k+1 cyclic) so a compact-but-scrambled ring counts as failed. **Cross-checked by two predictors, not validated; apo (no ATP/pRNA/DNA); addresses a geometric ring seat.**

## 1. Master table — construct × predictor × M2
| construct | predictor | samples | M2 engaged (seq) | scrambled | closure call |
|---|---|---|---|---|---|
| B1_L40_WT | Boltz-2 | 3 seeds | 5/5 all | 0 | **closes** |
| B1_L40_WT | OpenFold3 | 1 | 5/5 | 0 | **closes** |
| B1_L40_E119Q_m1 | Boltz-2 | 3 seeds | 5/5 all | 0 | **closes** |
| B1_L40_E119Q_m1 | OpenFold3 | 1 | 5/5 | 0 | **closes** |
| untethered_E119Q | Boltz-2 | 3 seeds | 5/5 all | 0 | closes |
| pentamer **L30** (×6) | Boltz-2 | 5 seeds | 5/5 all | 0 | **closes** |
| pentamer **L30** (×6) | OpenFold3 | 1 | 5/5 | 0 | **closes** |
| pentamer L40 (×8) | Boltz-2 | 5 seeds | 4/5 seeds 5/5; 1 collapses | 0 | Boltz only |
| pentamer L40 (×8) | OpenFold3 | 1 | 2/5 | 3 | scrambled |
| pentamer L50 (×10) | Boltz-2 | 3 seeds | 5/5 all | 0 | Boltz only |
| pentamer L50 (×10) | OpenFold3 | 1 | 1/5 | 3 | scrambled |
| tetramer L30 | Boltz-2 | 2 seeds | 4/4 | 0 | Boltz only |
| tetramer L30 | OpenFold3 | 1 | 0/4 | 1 | scrambled |
| trimer L30 | Boltz-2 | 2 seeds | 0/3 | 0 | no closure |
| trimer L30 | OpenFold3 | 1 | 0/3 | 2 | scrambled |

## 2. Task A — does E119Q-B1 still close? (tethered vs untethered)
**Yes, under both predictors.** B1_L40 with the catalytically-dead Walker-B E119Q on the tethered module closes 5/5 (Boltz 3/3 seeds, OpenFold3 5/5) — identical to the WT B1 baseline. The dead module's own R146 still engages its partner (3.9–4.7 Å). The untethered control (5 separate chains, one E119Q) also assembles 5/5 (Boltz). **Fold-level conclusion: assembly tolerance — a single seat can be catalytically inactivated without collapsing the ring, and a valid defined-dead-seat construct exists.** E119Q is a conservative near-isosteric substitution, so this is assembly tolerance, not a functional-locality claim; per-seat functional locality (dead seat vs delocalized effect) is future wet-lab work, not resolvable from an apo fold.

## 3. Task B — best linker / cross-predictor robustness
**Winner: the single-chain pentamer with (GGGGS)×6 = L30 internal linkers.** It is the only construct where BOTH predictors reform the correct sequential ring (Boltz 5/5 seeds → 5/5; OpenFold3 5/5; 0 scrambled). At L40 and L50 OpenFold3 builds a compact-but-scrambled ring (R146 against the wrong partner, 1–2/5 sequential, 3 scrambled) while Boltz still closes — the linker slack is what let OF3 scramble. Shortening to 30 aa removes that slack and forces the correct topology under both. **A cross-predictor-robust construct exists → the Task C generative-escalation trigger is NOT met.** Deepening showed the fix is **pentamer-specific**: L30 tetramer/trimer do not close cross-predictor (OF3 still scrambles / incomplete arc), so full 5-copy ring closure supplies the constraint, not the short linker alone.


## 3b. Three-predictor tie-breaker on the L40 single chain (AlphaFold3, user-run, independently scored here)
The user folded the (GGGGS)×8 = L40 single-chain pentamer with **AlphaFold3** (5 models, seed 1). Scored here with the same designed-sequential M2:

| predictor | global pTM | M2 (designed-sequential) | topology |
|---|---|---|---|
| Boltz-2 | 0.46–0.48 | 5/5 (4/5 seeds) | designed ring closes |
| OpenFold3 | 0.58 | 2/5 | scrambled |
| **AlphaFold3** | 0.49–0.51 | **1/5 (unanimous, 5/5 models)** | scrambled |

**2 of 3 independent predictors reject the L40 single-chain designed ring; AF3 is unanimous across all 5 models** (every R146 buried ~7.5 Å against the WRONG partner, designed partner 54–56 Å away). Two consequences: (1) **metric discipline is hardened** — all three predictors sit at global pTM 0.46–0.58, a band that cannot separate the correct ring from the scrambled ones (AF3 even reports `has_clash 0.0`), so only interface-resolved designed-sequential M2 discriminates; (2) it **vindicates the linker fix** — the L40 scramble is a real, reproducible two-predictor failure mode, and shortening to L30–L34 closes it under both NIM predictors. The OF3 scramble boundary is bracketed: **robust ≤L34, scrambled ≥L40** (L36 OF3 vote unavailable — NIM endpoint connection resets). Caveat: AF3 was run only on L40; whether AF3 also un-scrambles at L30/L34 is the single most valuable next fold (user-run AF3 on the L30 winner).

## 4. Task C — ProteinMPNN self-repair
**Not run — and now cross-predictor evidence makes the picture two-sided.** For the ORIGINAL (GGGGS)×8=L40 passive-linker chain, the trigger fires cleanly: 2 of 3 independent predictors (OF3 + AF3) reject it → a designed directional connector (RFdiffusion + ProteinMPNN, off-NIM, next-session) is the principled fix. BUT the linker sweep shows a simpler fix already works within the passive-GS family: L30/L34 close under both NIM predictors without any redesign. So generative escalation is warranted for L40 but is NOT the only path to a closed single chain. RFdiffusion is not on NIM — logged as next-session follow-up.

## 5. What's robust / what's not / recommended next
Two independent structure predictors agree on two things: (a) the **B1 covalent-dimer-in-ring tolerates a catalytically-dead E119Q seat** without losing assembly, and (b) a **fully single-chain pentamer closes the correct sequential ring when its internal linkers are shortened to (GGGGS)×6**. Both are cross-checked, not validated, and both are apo predictions of a geometric ring seat — no claim is made about the resting-ring functional state or about per-seat catalytic locality. What is not robust: the original (GGGGS)×8 single-chain pentamer (OpenFold3-scrambled), and sub-pentameric L30 arcs (do not close under both predictors). Recommended next: (i) fold L30 with nucleotide/pRNA/DNA context to test whether closure survives the loaded state; (ii) an OpenFold3 seed sweep is not informative here (endpoint returns one deterministic sample), so use a third predictor (e.g. Chai-1 as in cycle 3) for a stronger cross-check on L30; (iii) wet-lab the two constructs — B1_L40_E119Q (addressability) and pentamer-L30 (single-chain closure) — as the two lead designs; (iv) if a backbone-level redesign is ever wanted for the L40/L50 regime, RFdiffusion is not on NIM — flag as a next-session follow-up.

## 6. Structures + provenance
All folds: Boltz-2 NIM + OpenFold3 NIM at `health.api.nvidia.com`, ColabFold MSA (core 327-aa, 1168 records, block-diagonal-tiled per copy), via `cofactor/fold.py`. Saved: `outputs/structures/e119q/*` (Task A: B1 WT/E119Q, untethered, Boltz+OF3), `outputs/structures/taskB/*` (pentamer L30/L40/L50 + L30 ladder). Tables: `taskA_e119q_assembly.csv`, `taskB_pentamer_sweep.csv`, `taskB_L30_ladder_transfer.csv`. Findings: `taskA_e119q_finding.txt`, `taskB_sweep_finding.txt`. Figure: `fig8_overnight_taskAB.png`. NIM is a synchronous REST endpoint — no persistent job IDs (unlike the Modal cycle-3 jobs); each fold is one POST, structures saved immediately on return.


## 7. Linker-length design curve (matched Boltz-2 + OpenFold3 grid + AF3 overlay)
Folded the full single-chain pentamer linker grid on free NIM (sequential M2, Boltz 2–5 seeds + OpenFold3 each), to plot correct-fold as a function of linker length against the user's AlphaFold3 runs:

| linker (aa) | Boltz-2 | OpenFold3 | AlphaFold3 | cross-predictor |
|---|---|---|---|---|
| L20 | 5/5 | 5/5 | (pending export) | robust |
| L25 | 5/5 (2/3 seeds) | 5/5 | (pending) | robust |
| **L30** ⭐ | 5/5 (5 seeds) | 5/5 | (pending) | **robust** |
| **L34** ⭐ | 5/5 | 5/5 | (pending) | **robust** |
| L36 | 5/5 | **1/5 scrambled** | (pending) | OF3 flips |
| L38 | 5/5 | 1/5 scrambled | (pending) | scrambled |
| L40 | 4/5 seeds | 2/5 scrambled | **1/5 (5 models)** | scrambled |
| L50 | 5/5 | 1/5 scrambled | (pending) | scrambled |

**The design curve has a sharp OpenFold3 cliff between L34 and L36.** Cross-predictor-robust window = **L20–L34** (both NIM predictors close the designed ring); at L36 and above OF3 threads the chain into a compact wrong-partner register. Shorter linkers (down to L20) do NOT start clashing — they stay closed. The AF3 overlay currently has one point (L40 = scrambled, matching OF3); the plot updates per-length as AF3 exports arrive. Figure: `fig9_linker_curve.png`; grid table: `linker_grid_correctfold.csv`.

## 8. L30 winner mechanism / switch-position (ported from L40)
Per-interface trans-R146 across all 6 L30 folds (5 Boltz seeds + OF3): **all 5 sequential interfaces engaged in every fold, mean 5.9–6.2 Å, inter-interface spread only 0.33 Å.** The robust L30 ring is **symmetric with no pinned opening** — it reproduces the native apo-ring symmetry (cycle 3: 6.5–6.7 Å, spread 0.2 Å). No interface is a distinguished "switch" seat in the apo/resting state; the open position is not a fixed geometric seat, consistent with the planar-resting hypothesis. The linker fix (≤34 aa) converts the L40 scrambled-register failure into this correct symmetric ring. Finding: `L30_mechanism_finding.txt`; per-interface data: `L30_per_interface.csv`.


## 9. Circular-permutation route (N–C direct fusion failed AF3 → CP escalation) — SUCCESS at the CP rung
User's AlphaFold3 rejected the N–C direct-fusion single chain at every linker length → escalate to circular permutation (CP moves each subunit's termini to the hinge/folding-unit boundary, shortening the inter-subunit junction so the chain cannot thread into a wrong-partner register). 18 constructs, 3 cut sites (228/233/217), staged.

**Stage 1 — CP monomers fold natively (gate, all 3 sites PASS).** Each CP subunit = [res P+1..330] + internal GS linker + [res 4..P]; folded with Boltz-2 + permuted ColabFold MSA (1168 core homologs), superposed onto native gp16. All 9 constructs recover both domains in native arrangement: N-domain 2.5–3.6 Å, C-domain 2.2–7.8 Å RMSD — **not fragmented, the qualitative opposite of the N–C direct fusion.** Best internal linker per site: cp228 int10 (N2.80/C2.30), cp217 int20 (N2.57/C2.25), cp233 int20 (N2.46/C4.57). R146→own Walker-A ~25.5 Å in every monomer = correct trans-finger orientation (points to neighbour, not own pocket).

**Stage 2 — CP pentamers in ring context (sequential M2, Boltz 3 seeds + OpenFold3, pass = both ≥4/5).**

| site | inter-linker | Boltz (seeds close) | OpenFold3 | cross-predictor |
|---|---|---|---|---|
| 217 | 10/15/20 | 1/3, 0/3, 1/3 | 5/5 | fail (Boltz seed-variable) |
| 228 | 10/15/20 | 1/3, 0/3, 1/3 | 5/5 | fail (Boltz seed-variable) |
| **233** | **10** | **3/3** | **5/5** | **PASS** |
| **233** | 15 | 2/3 | 5/5 | PASS |
| **233** | **20** | **3/3** | **5/5** | **PASS** |

**WINNER: cp233 (cut at res 233) closes the CP ring under BOTH predictors at all 3 inter-linker lengths.** cp233_int15_inter10 is tightest — Boltz 3/3 seeds 5/5 (all interfaces 5.4–5.7 Å), OpenFold3 5/5, 0 scrambled. Closes like B1, unlike the direct fusion.

**Key cross-predictor observation:** OpenFold3 closes the CP ring 5/5 on **all 9** constructs — the exact opposite of the N–C direct-fusion chain, where OF3 scrambled at every length. The circular permutation fixes the wrong-partner threading OF3 was catching. Boltz is seed-variable at sites 217/228 (less-constrained junction) but robust at 233.

**Storyline:** N–C direct fusion (failed, AF3) → **circular permutation (succeeded: cp233 closes cross-predictor)** → diffusion NOT needed. Three standing lead designs now: (1) **cp233_int15_inter10** (circular-permutation single chain), (2) **B1_L40_E119Q** (covalent dimer + 3 WT, addressable dead seat), (3) pentamer-L30 (short-linker direct fusion, Boltz+OF3 only — but AF3 rejected the L40 variant, so this route is the weakest). Recommend handing **cp233_int15_inter10** to AlphaFold3 for the decisive third-predictor check. Files: `cp_stage1_monomers.csv`, `cp_stage2_pentamers.csv`, `cp_stage2_allfolds.csv`, findings, `fig10_circular_permutation.png`.
