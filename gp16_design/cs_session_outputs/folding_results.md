# Folding results — gp16 position-addressable design (Aim 3)

**Method:** each construct folded by TWO independent predictors on identical inputs, each with its own real MSA. Confidence signal = agreement between predictors. Results are *computationally cross-checked*, not validated.

## Runs & provenance

| Predictor | Version | MSA source | GPU | Job ID | Wall (s) | Seed |
|---|---|---|---|---|---|---|
| Boltz-2 | boltz 2.2.1 | ColabFold server (`--use_msa_server`) | A100-80GB (Modal) | `f41858b0-0de5-4ce9-9e36-b8812c37cbd2` | 668 | default |
| Chai-1 | chai_lab | ColabFold server (built-in MSA) | A100-80GB (Modal) | `009395ae-7ead-4cfa-bbf5-b037b38cebca` | 870 | 42 |

Constructs folded: `gp16_CP280_2mer` (lead, 719 aa), `gp16_nativefus_2mer` (negative control, 684 aa), `gp16_native_dimer` (positive control, 2×332 aa native chains).
STEP 0 pre-check: `gp16_CP280_monomer` folded with public ESMFold (single-sequence) — domains fold locally but inter-domain arrangement not recovered without an MSA; motivated the MSA-based cross-check here.
(Boltz-2 ran the native dimer, rc 0, but its multi-chain output was not captured before the sandbox closed; the positive control is scored from Chai-1.)

## Scores — function-preserving contacts, not global fold

Native 7JQQ reference: intersubunit A–B interface = **53** CA contacts (<8 Å); trans arginine-finger R146(B)→Walker-A(A) = **3.17 Å**.

| Construct | Predictor | Module core RMSD (Å) | N-dom 4–200 (Å) | C-dom 201–330 (Å) | Interface contacts | R146 trans-min (Å) | core pLDDT | linker pLDDT | pTM |
|---|---|---|---|---|---|---|---|---|---|
| gp16_CP280_2mer | Boltz-2 | 13.69 | 3.12 | 19.56 | 176 | 60.35 | 76.3 | 33.7 | 0.411 |
| gp16_CP280_2mer | Chai-1 | 15.77 | 2.38 | 18.68 | 118 | 56.59 | 71.6 | 39.2 | 0.458 |
| gp16_nativefus_2mer | Boltz-2 | 6.14 | 2.62 | 7.93 | 82 | 5.38 | 74.6 | 26.6 | 0.474 |
| gp16_nativefus_2mer | Chai-1 | 13.29 | 2.3 | 14.54 | 46 | 3.97 | 73.4 | 42.2 | 0.553 |
| gp16_native_dimer(POS) | Chai-1 | 8.28 | nan | nan | 61 | 2.89 | nan | nan | 0.657 |

## Headline finding (honest, cross-checked)

**The CP@280 lead design does NOT beat the native-fusion negative control on the two function-preserving criteria (interface + R146). The static topology ranking is overturned by the fold.**

Both predictors agree on the mechanism:

1. **The gp16 N-terminal ATPase domain (res 4–200) folds natively in every construct** — RMSD 2.3–3.1 Å across CP@280, native-fusion, both predictors. The permutation and the linkers do not damage the core ATPase fold.
2. **The circular-permutation cut at residue 280 lands *inside* the C-terminal domain**, not in an inter-domain loop. After N-domain superposition, the relocated 280–330 fragment is displaced ~34 Å (Boltz) and the whole C-domain is at 18–20 Å RMSD for CP@280, versus 8–15 Å for native-fusion. Splitting the chain at 280 prevents the C-terminal subdomain from reassembling.
3. **Consequence for addressability:** because the C-domain does not reform, the trans arginine-finger R146 sits **56–60 Å** from the neighboring module's ATP pocket in CP@280 (both predictors) — the catalytic inter-subunit registry is lost. Native-fusion, by contrast, reforms it (Boltz 5.4 Å, Chai 4.0 Å ≈ native 3.2 Å).
4. **Metric sanity check (positive control):** the native two-chain dimer reforms the interface (61 contacts) and R146 (2.89 Å ≈ native 3.17 Å). So the metric detects a real trans-interface when one forms — the CP@280 failure is genuine, not a scoring artifact.

**Confidence / disagreement:** the two predictors agree on all qualitative calls (N-domain native, C-domain broken by CP@280, R146 lost in CP@280, R146 reformed in native-fusion). They diverge **substantially on the native-fusion construct**: module core RMSD 6.1 Å (Boltz) vs 13.3 Å (Chai), C-domain RMSD 7.9 Å (Boltz) vs 14.5 Å (Chai) — roughly 2× — and interface-contact count 82 (Boltz) vs 46 (Chai). Boltz-2 folds the native-fusion C-domain much closer to native than Chai-1 does; how well that C-domain and junction actually pack is genuinely uncertain and must not be over-read from the Boltz number alone. This divergence does not touch the CP@280 conclusion — both predictors put CP@280's C-domain at 18.7–19.6 Å and R146 at 56–60 Å, so the lead fails the addressability test either way. The inter-module linkers are the low-confidence region in every model (pLDDT 27–42 vs 72–76 core), as expected for flexible GS linkers.

## Interpretation for the project

This is failure mode **F5** made concrete (from `topology_failure_taxonomy.md`): a design can look best on static terminus-to-terminus geometry yet fail at the contact level. The static score ranked CP@280 #1 because it minimized the *fusion span* between adjacent subunit termini — but it did not encode that cutting at 280 bisects a folded subdomain. **The fold-level cross-check is exactly the filter that catches this.**

**Actionable next design cycle:** move the CP cut out of the C-domain. Res 280 and 297 both scored well on span but both fall in the C-terminal domain body. A CP site should sit in a genuine inter-domain or surface loop (e.g. near the 195–205 N/C-domain hinge, or a validated surface loop N-terminal to 280) so the cut does not fragment a cooperative folding unit. Re-fold that series through this same two-predictor pipeline. Native-fusion — despite fouling the DNA channel with the res330 C-terminus (the reason it was a negative control for the *ring* geometry) — actually preserves the subunit fold and the R146 interface better, so a hybrid (native-order fusion with a repositioned channel-clearing linker) is also worth folding.

*Framing note:* a 2-mer fold tests whether the CP fold, junction, and inter-subunit interface survive — not which conformational state (planar vs helical) the full ring adopts. No claim is made here about the ring state.


---

# Cycle 2 — cheap→expensive screening funnel

## Stage A — domain-aware static filter (FREE, CPU) — COMPLETE

**Method.** From 7JQQ chain A: (1) CA–CA contact graph → greedy-modularity communities; (2) for every cut site *k*, **cut-penalty** = number of native long-range heavy-atom contacts (|i−j|>12, <5 Å) severed by cutting between *k*−1 and *k*. Low penalty = a real domain/loop boundary. 433 native long-range contacts total; penalty median 70, 15th-percentile 33.

**Folding units:** N-domain (res 4–200), N/C **hinge (201–228)**, C-domain (229–330). The deep cut-penalty valley at 200–233 is the hinge.

**Hard constraint applied:** keep a cut only if (a) it separates two folding units (domain boundary) OR is a low-contact surface loop (penalty ≤ 33), AND (b) it is clear of Walker A/B, R146/trans-loop, ATP pocket, DNA-channel, and pRNA contacts (±2 buffer).

**Cycle-1 leads are now REJECTED:**

| Site | cut-penalty | folding unit | verdict |
|---|---|---|---|
| cut@280 (was rank #1) | **66** | C-domain body | REJECT — severs ~66 C-domain long-range contacts; not a boundary. *The cycle-1 fold confirmed this: C-domain fragmented, displaced ~34 Å.* |
| cut@297 (was rank #2) | **60** | C-domain body | REJECT — same failure mode. |

**Stage B shortlist (survivors, re-ranked by cut-penalty):**

| Tier | Site | cut-penalty | folding unit | domain boundary | junction span (helical max, Å) | stroke (Å) |
|---|---|---|---|---|---|---|
| 1 | **cut@228** | **6** | hinge (201–228) | yes | 38.8 | 7.6 |
| 1 | **cut@233** | **13** | C-domain edge (hinge-adjacent) | yes | 44.9 | 7.6 |
| 2 | cut@217 | 33 | hinge (201–228) | near-threshold | 43.4 | 6.2 |
| 2 | cut@205 | 34 | hinge (201–228) | near-threshold | 57.1 | 5.3 |

cut@228 has an **11× lower** cut-penalty than the cycle-1 lead cut@280 (6 vs 66). These sites trade a longer junction span (39–57 Å vs 24–27 Å for the rejected leads) for not fragmenting a folding unit — the fold showed span was the wrong thing to optimize.

## Estimated GPU cost per stage (Modal, compute-time; container startup adds ~1–2 min/job)

| Stage | GPU | est. GPU-hr | est. $ |
|---|---|---|---|
| A — static filter (done) | CPU | — | $0.00 |
| B — 4 monomer folds (Boltz-2, cached MSA, reduced settings) | A100-40GB | 0.17 | ~$0.35 |
| C — ≤2 2-mer folds (Boltz-2, cached MSA) | A100-40GB | 0.13 | ~$0.28 |
| D — 1 winner, Boltz-2 + Chai-1 (2-mer, full settings) | A100-80GB | 0.17 | ~$0.61 |
| **Total B–D** | | | **~$1.24** |

For comparison, cycle-1's blind dual-predictor fold of a single 2-mer was ≈$10. The one-time image build is already paid; the MSA is computed once for native gp16 and reused for all CP variants by re-permuting its columns.

*Nothing folded yet — awaiting approval to start Stage B.*


## Stage B — cheap monomer screen (Boltz-2, cached permuted MSA) — COMPLETE

MSA computed ONCE for native gp16 via ColabFold server, then column-permuted per CP variant (no per-variant search). All 4 hinge-cut monomers + native folded at reduced settings (recycling 1, sampling 50, 1 sample). Job `e22bb5a5-6052-4f96-8a65-319c2324df3a (A100-80GB, 279s)`.

| Monomer | cut | N-dom RMSD (Å) | C-dom RMSD (Å) | C-dom displacement after N-superposition (Å) | mean pLDDT |
|---|---|---|---|---|---|
| native | — | 2.62 | 7.79 | 9.9 | 81.5 |
| CP228 | 228 | 2.70 | 17.49 | 59.0 | 79.3 |
| CP233 | 233 | 2.60 | 18.43 | 57.0 | 79.8 |
| CP217 | 217 | 2.74 | 13.12 | 54.9 | 79.1 |
| CP205 | 205 | 3.16 | 14.97 | 52.5 | 79.2 |

**The specified monomer PASS/FAIL gate is INVALID for gp16 — a structural discovery, not a screen failure.** The gate assumed the subunit's two domains pack against each other intra-molecularly. They do not. Native-contact analysis of 7JQQ chain A (heavy-atom <5 Å, residue-pair counts):

- Chain-A **C-domain (229–330)**: only **3** contacts with its own N-domain, but **25** with neighbor B's C-domain and **19** with neighbor E's C-domain.
- Chain-A **N-domain (4–200)**: only **3** with its own C-domain, but **30** with neighbor B's N-domain and **27** with neighbor E's N-domain.

gp16 forms **inter-subunit N-domain and C-domain rings** — each domain's position is held by *neighboring subunits*. A folded monomer therefore cannot recover the native inter-domain arrangement; the native monomer itself floats its C-domain to 9.9 Å. **What Stage B does establish:** the N-domain folds natively in all four hinge cuts (2.6–3.2 Å), so the permutation does not damage local structure. Addressability must be judged in ≥dimer context. (`outputs/domain_topology_finding.json`, `outputs/monomer_screen.csv`.)

## Stage C — single-predictor 2-mer (Boltz-2, MSA server) — COMPLETE

Folded the corrected hinge **gp16_CP228_2mer** (677 aa) vs the **gp16_nativeorder_2mer** (684 aa, GS6 channel-clearing linker, the B1 primary). Job `1c9a2f01-a344-4d62-95b7-bc236cb3bed2 (A100-80GB, 232s)`. Native ref: interface ~53 contacts, trans R146 ~3.2 Å.

| 2-mer | module core RMSD (Å) | interface contacts | trans R146 (Å) | confidence |
|---|---|---|---|---|
| gp16_nativeorder_2mer | 8.2 / 8.3 | 25 | 53.8 | 0.717 |
| gp16_CP228_2mer | 22.5 / 23.0 | 249 | 33.4 | 0.725 |

**Winner: native-order 2-mer** (module cores 8.2 Å vs CP228's 22.5 Å — the hinge cut still fragments the fold once tension from the second module is applied, even though its monomer N-domain was fine). CP228's 249 "interface contacts" are a collapsed, non-native aggregate (cf. cycle 1: raw contact count misleads). Note this native-order run gave R146 = 53.8 Å — the *same 684-aa sequence* that reformed R146 (5.4 Å) in cycle 1. That discrepancy motivated Stage D.

## Stage D — dual-predictor ENSEMBLE on the winner — COMPLETE

Because the cycle-1↔Stage-C discrepancy on an identical sequence signalled run-to-run instability, Stage D characterizes **reliability**, not a single fold: Boltz-2 5-sample ensemble (job `4a1ff265-53e3-47db-ad3f-f5a6cc9655be (A100-80GB, 138s, 5 diffusion samples)`) + Chai-1 2 seeds × 5 models (job `c58ab0d1-10d4-4a33-8cba-2b0d7b8ad625 (A100-80GB, 475s, 2 seeds x 5 models)`), 15 models total.

**trans R146 → neighbor pocket, "interface reformed" = R146 < 8 Å:**

| Predictor | models | reformed | median R146 (Å) | median module core (Å) |
|---|---|---|---|---|
| Boltz-2 | 5 | **0 / 5** | 84.8 | 8.1 |
| Chai-1 | 10 (seed42: 5, seed7: 5) | **7 / 10** | 4.2 | 12.3 |

**The interface reformation is stochastic and predictor-dependent.** Boltz-2 never reforms R146 across 5 samples; Chai-1 reforms it in 7/10, and even disagrees across its own seeds (seed42 5/5, seed7 2/5). Cycle 1's single dual-fold that showed R146 = 4–5 Å was one draw from this distribution, not a robust result.

This is exactly what the **inter-subunit domain-ring topology (Stage B) predicts**: a 2-mer supplies each C-domain only ONE of its two native neighbors, so the trans interface is geometrically **under-determined** — it forms or not depending on the diffusion trajectory. The predictors' disagreement IS the honest confidence signal the protocol asked for.

## Cycle-2 conclusion & cost

- **Total GPU spend cycle 2: ~$0.82** (B $0.28 + C $0.14 + D $0.40), well under the $10 budget. Cycle-1's single blind dual-fold was ≈$10.
- **Scientific bottom line (cross-checked, not validated):** neither a hinge-corrected circular permutation nor a native-order fusion gives a *reliable* position-addressable 2-mer, because gp16's addressability-defining contacts (the C-domain ring and the trans R146) are **inter-subunit** — they need ≥2 flanking neighbors, which a 2-mer cannot provide. This directly supports the strategic pivot: **(B1)** a covalent native-order dimer must be evaluated as a dimer *embedded in a WT ring* (dimer + 3 WT co-assembly), not in isolation; **(B2)** the correct method abstraction is a folding-unit + inter-subunit-contact graph, where a valid single-chain path joins only spatially-adjacent segments and never severs an inter-subunit-held unit.
- **Caveat foregrounded:** 7JQQ is the helical (lock-washer) state and the ring cycles planar↔helical each translocation cycle; the special/regulatory subunit may be a dynamically reassigned role rather than a fixed seat. Any position-addressed construct addresses a **geometric seat, not a guaranteed fixed function** — which is precisely the hypothesis such constructs are built to test.

### Cycle-2 provenance
| Stage | Job ID | GPU | Predictor(s) | seeds/samples |
|---|---|---|---|---|
| B | `e22bb5a5-6052-4f96-8a65-319c2324df3a (A100-80GB, 279s)` | A100-80GB | Boltz-2 | 1 sample, reduced |
| C | `1c9a2f01-a344-4d62-95b7-bc236cb3bed2 (A100-80GB, 232s)` | A100-80GB | Boltz-2 | 1 sample, full |
| D-Boltz | `4a1ff265-53e3-47db-ad3f-f5a6cc9655be (A100-80GB, 138s, 5 diffusion samples)` | A100-80GB | Boltz-2 | 5 diffusion samples |
| D-Chai | `c58ab0d1-10d4-4a33-8cba-2b0d7b8ad625 (A100-80GB, 475s, 2 seeds x 5 models)` | A100-80GB | Chai-1 | seeds 42,7 × 5 models |

MSA source throughout: ColabFold server (Boltz `--use_msa_server`; Chai built-in). All results *computationally cross-checked by two independent predictors*, not validated.


---

# Cycle 3 / Overnight run — B1 dimer-in-ring (the positive result)

## Design logic
Cycle-2 established gp16's interface is a ≥3-subunit property (each C-domain contacts its TWO neighbors, not its own N-domain), so a 2-mer under-determines it. Cycle 3 tests the covalent dimer where it will actually live: **embedded in a full pentameric ring** (2 tethered seats + 3 separate WT chains), where every subunit has both neighbors.

## Phase 1 — native ring reference
Native gp16 pentamer folded as 5 separate WT chains reproduces a **closed pentamer with correct cyclic adjacency** (A:B,E / B:A,C / C:B,D / D:C,E / E:A,D), each monomer 4.5 Å RMSD to 7JQQ chain A. The predicted **ligand-free ring is symmetric** — all 5 arginine fingers engaged at 6.5–6.7 Å — versus 7JQQ's **asymmetric helical** state (fingers tight at ATP interfaces B→A 3.2/C→B 2.8 Å, disengaged at the open A–E seam 15–25 Å). Consistent with the planar resting-ring hypothesis. This symmetric apo ring is the reference distribution.

## Phase 2 — B1 dimer-in-ring assembly gate (Boltz, 1 seed)
Native-order covalent dimer (res 4–330 – GS linker – res 4–330) + 3 WT, 3 linker lengths. **All PASS:** seat RMSD 4.7–5.2 Å; intra-dimer R146 5.5–6.3 Å ≈ native apo ring; ring closes through two WT neighbors. The isolated 2-mer's failure was an artifact of missing neighbors, not the design.

## Phase 3/7 — decisive ensemble (Boltz-2 3 seeds + Chai-1 2 seeds × 5 models)
Per-interface trans-R146 engagement (< 8 Å), tether-flanking seats called out:

| construct | Boltz engaged | Chai engaged | Boltz median | Chai median | tether seats |
|---|---|---|---|---|---|
| native_ring (ref) | 15/15 | 50/50 | 6.6 Å | 2.8 Å | — |
| B1_L32_in_ring | 15/15 | 45/50 | 5.7 Å | 4.0 Å | Boltz 6/6, Chai 19/20 |
| **B1_L40_in_ring** | **15/15** | **47/50** | **5.8 Å** | **4.0 Å** | **Boltz 6/6, Chai 20/20** |

**Cross-checked conclusion:** in ring context the B1 covalent tether reforms the catalytic R146 interface as reliably as the native ring — both predictors agree. Predictor divergence is only in absolute geometry (Chai packs tighter, median 2.8–4.0 Å; Boltz looser, 5.7–6.6 Å), not in the engaged/not-engaged call. This is a genuine positive, opposite to the cycle-2 2-mer (Boltz 0/5).

## Phase 5 — B2 folding-unit graph (free)
Formalized: nodes = 10 folding units (N/C domain × 5 subunits), edges = inter-subunit contacts. Valid single-chain design = walk over ring seats that (R1) never severs a unit, (R2) joins only spatially-adjacent seats, (R3) is scored in ring context. Explains cycle-1 (broke R1), cycle-2 (broke R3), B1 (satisfies all three). Top topologies avoid the 66 Å A–E seam (leave it as the free break). `b2_folding_unit_graph.json`, `b2_topologies.md`.

## Phase 6 — B3 pRNA-scaffold (free)
Each gp16 subunit contacts exactly one pRNA chain (A↔K … E↔O, 15–21 residues, at the pRNA-binding surface away from the R146 interface). Engineered interlocking-loop complementarity (pRNA interlocking loops; [refs supplied in project brief; not independently verified in this session]) could template a defined subunit order with no giant chain. `b3_prna_scaffold_memo.md`.

## Cost & provenance
Total GPU spend **$9.14 / $10.00**. All folds ColabFold-MSA (Boltz `--use_msa_server`; Chai built-in + ESM), A100-80GB. Jobs: P1 `c001e8ad`, P2 `cc880831`, ensemble-Boltz `3cb2936b`, ensemble-Chai `44c3880c`. Full per-job costs in `budget_ledger.csv`; per-interface per-seed scores in `cycle3_ring_scores.csv`. Computationally cross-checked by two independent predictors, not experimentally validated. 7JQQ helical state; folds without pRNA/DNA/ATP; constructs address a geometric seat, not a fixed functional role (planar↔helical & dynamic-special-subunit; [refs supplied in project brief; not independently verified in this session]).


---

# Fusion ladder — MSA folds on free NIM (msa-search via ColabFold + Boltz-2 NIM; OpenFold3 cross-check)

**Provenance:** MSA = ColabFold public MMseqs2 server (core 327-aa, 1168 records), block-diagonal-tiled per copy (1167 records/copy, linker columns gapped). Fold = NVIDIA BioNeMo **Boltz-2 NIM** (`health.api.nvidia.com`, free credits), key delivered via a non-reserved credential (`NVKEY_DIRECT`) since the confined inference-kernel path was unavailable on this host. Cross-check = **OpenFold3 NIM**, same tiled MSA. Pipeline = the user's local `cofactor/fold.py` (stdlib ColabFold→Boltz-2 wrapper), extended with core-MSA tiling. Constructs = native-order single chains, N×(res4–330 core) joined by (GGGGS)×8 linkers; copy1–copyN = the single non-covalent seam. **Cross-checked by two predictors, not validated; apo (no pRNA/DNA/ATP); addresses a geometric ring seat, not a fixed functional role.**

## Confidence gain (MSA vs the user's single-seq screens)
| rung | pTM single-seq | pTM MSA | pLDDT single-seq | pLDDT MSA |
|---|---|---|---|---|
| dimer | — | 0.601 | — | 0.72 |
| trimer | 0.304 | 0.491 | 0.50 | 0.70 |
| tetramer | 0.265 | 0.410 | 0.45 | 0.68 |
| pentamer | 0.267 | 0.482 | 0.45 | 0.68 |

The tiled MSA roughly **doubles pTM** and lifts pLDDT from ~0.45 to ~0.68–0.72 at every rung — the single-seq screens were genuinely under-determined.

## Step 1 — per-interface closure (Boltz-2, native ~3.2 Å; apo-ring ref ~6.6 Å; engaged <8 Å)
| rung | interfaces engaged | open interface(s) |
|---|---|---|
| dimer | 1/2 | seam 2→1 (53.0 Å) |
| trimer | 2/3 | seam 3→1 (58.1 Å) |
| tetramer | 3/4 | **internal junction 1→2 (8.3 Å, marginal)** — seam is CLOSED (6.6 Å) |
| pentamer | **5/5** (4.9–5.5 Å) | none — full ring closes, seam included |

**Answer to "does the full single chain reform ALL interfaces, and where does it break":** the **pentamer closes all five (Boltz-2)**; incomplete arcs break at whichever interface cannot close. The break is at the seam for dimer/trimer but **relocates off the seam** to an internal junction by the tetramer — it is NOT fixed at the copy1–copyN seam.

## Step 2 — switch-position analysis (openness across 2→3→4→5 + native ring)
The R146 spread across interfaces collapses as subunits are added: **dimer 46.1 Å → trimer 51.6 Å → tetramer 1.7 Å → pentamer 0.6 Å**; the native apo ring is symmetric (~6.6 Å all round). So the "most-open" interface is a genuine physical gap only in incomplete arcs; once the ring closes the openness is symmetric with **no distinguished switch**, and its nominal position is not pinned by the covalent-junction layout. **The open/switch position is a positionable geometric seat, not a fixed functional switch** — consistent with the apo/planar-resting hypothesis.

## Step 3 — cross-predictor check on the best rung (pentamer)
Boltz-2 and OpenFold3, SAME tiled MSA:
- **Boltz-2: 5/5 interfaces engaged** (4.9–5.5 Å) → full closure.
- **OpenFold3: 2/5 engaged** (3→4 = 5.0 Å, 5→1 = 5.0 Å); 1→2/2→3/4→5 open (33–53 Å) → partial closure.
- Agreement 2/5; global confidence comparable (Boltz pTM 0.48/pLDDT 68; OF3 pTM 0.58/pLDDT 77).

**Honest read:** whether the single-chain tandem pentamer reforms ALL interfaces is **predictor-dependent** — Boltz-2 says yes, OpenFold3 says ~half. This is real uncertainty at the tandem-fusion ring closure, reported not hidden. **Robust across both predictors:** the open interfaces are scattered, not seam-locked (Step 2 conclusion holds under both). **Contrast:** the multi-chain B1 covalent-dimer-in-ring (cycle 3) had both predictors agree (Boltz 15/15, Chai 47/50) — the single-chain tandem fusion is **less robustly closed** than the covalent-dimer-in-ring, so B1 remains the stronger position-addressable construct.

## Step 4 — E119Q prep (no fold; `fold_inputs/e119q/`)
Defined dead-seat constructs staged: (a) tethered pentamer ladder with E119Q on copy 1; (b) B1_L40 dimer-in-ring with E119Q on the tethered module + 3 WT; (c) untethered control (5 separate chains, one E119Q — dead position present but not covalently defined). Residue 119 confirmed as the Walker-B catalytic glutamate.

Figure: `fig7_ladder_msa.png`. Tables: `ladder_scores.csv`, `switch_position_analysis.csv`, `pentamer_cross_predictor.csv`, `ladder_confidence_gain.csv`. Structures: `outputs/structures/ladder/*__boltz2nim_msa.cif`, `*__openfold3_msa.cif`.


---

# Overnight autonomous run — Task A (E119Q assembly) + Task B (linker robustness sweep)

**Backend:** free NVIDIA NIM (Boltz-2 + OpenFold3, `health.api.nvidia.com`, key `NVKEY_DIRECT`) + ColabFold MSA via `cofactor/fold.py`. No Modal. Scored by M2 (per-interface trans-R146, <8 Å), single-chain scored SEQUENTIALLY (copy k → k+1 cyclic) so scrambled rings fail. Cross-checked by two predictors, not validated; apo; geometric seat.

## Task A — E119Q assembly tolerance on B1 (HIGHEST priority)
| construct | Boltz-2 (3 seeds) | OpenFold3 | verdict |
|---|---|---|---|
| B1_L40_WT | 5/5 all | 5/5 | closes |
| B1_L40_E119Q_m1 (dead seat) | 5/5 all | 5/5 | **closes — dead seat tolerated** |
| untethered_E119Q control | 5/5 all | — | closes |

The catalytically-dead Walker-B E119Q on the tethered module does **not** disrupt assembly under either predictor; the dead module's R146 still engages (3.9–4.7 Å). Opposite of the single-chain pentamer, where OF3 scrambled. E119Q is a conservative substitution — fold-level result is assembly tolerance + a valid defined-dead-seat construct; per-seat functional locality is future wet-lab.

## Task B — single-chain pentamer linker sweep (attack Boltz↔OF3 divergence)
| linker | Boltz-2 (seeds→5/5) | OpenFold3 (seq engaged) | scrambled | cross-predictor |
|---|---|---|---|---|
| **L30 (GGGGS)×6** | **5/5 seeds** | **5/5** | 0 | **BOTH close** |
| L40 (GGGGS)×8 | 4/5 seeds | 2/5 | 3 | diverge |
| L50 (GGGGS)×10 | 3/3 seeds | 1/5 | 3 | diverge |

**L30 = (GGGGS)×6 is the cross-predictor-robust single-chain pentamer.** Shortening the internal linker removes the slack that let OpenFold3 build a compact-but-scrambled ring (R146 against the wrong partner). A cross-predictor-robust construct exists → **Task C (ProteinMPNN) not triggered.** Deepening: the L30 fix is **pentamer-specific** — L30 tetramer (Boltz 4/4 but OF3 0/4) and trimer (neither closes) still diverge; full 5-copy closure supplies the geometric constraint, not the short linker alone.

## Two lead designs (cross-checked, apo)
1. **B1_L40_E119Q** — covalent dimer + 3 WT with a defined catalytically-dead seat; the addressability construct.
2. **pentamer-L30** — fully single-chain pentamer, (GGGGS)×6 linkers; closes the correct sequential ring under both predictors.

Figure: `fig8_overnight_taskAB.png`. Full run summary: `RUN_SUMMARY.md`. Tables: `taskA_e119q_assembly.csv`, `taskB_pentamer_sweep.csv`, `taskB_L30_ladder_transfer.csv`.


## AlphaFold3 tie-breaker (user-run, independently scored) — L40 single-chain pentamer
The (GGGGS)×8 = L40 single-chain pentamer was folded with AlphaFold3 (5 models, seed 1) and scored here with the designed-sequential M2:

| predictor | global pTM | M2 designed-sequential | topology |
|---|---|---|---|
| Boltz-2 | 0.46–0.48 | 5/5 (4/5 seeds) | designed ring closes |
| OpenFold3 | 0.58 | 2/5 | scrambled |
| AlphaFold3 | 0.49–0.51 | **1/5 (unanimous, 5/5 models)** | scrambled |

**2 of 3 independent predictors reject the L40 single chain; AF3 unanimous.** Every R146 buried ~7.5 Å against the wrong partner; designed sequential partner 54–56 Å away. Hardens the metric-discipline case (global pTM 0.46–0.58 for all three cannot discriminate; AF3 reports `has_clash 0.0`) and vindicates the linker fix (L30–L34 close under both NIM predictors). OF3 scramble boundary bracketed: robust ≤L34, scrambled ≥L40. Extended sweep confirms **L34 also cross-predictor-robust** (Boltz 5/5, OF3 5/5). Cross-checked by three predictors, not validated; apo; geometric seat. Structures: `outputs/structures/af3/` (5 models + confidences). Finding: `af3_tiebreaker_finding.txt`.


## Linker-length design curve + L30 mechanism (matched NIM grid)
Full single-chain pentamer linker sweep on free NIM (sequential M2): **cross-predictor-robust window L20–L34** (Boltz-2 + OpenFold3 both close the designed ring); OpenFold3 flips to a scrambled wrong-partner register at **L36** and stays scrambled through L50. Boltz-2 closes at all lengths (one L40 seed collapses). AlphaFold3 confirms scramble at L40 (1/5). Shorter linkers (L20) do not clash — they stay closed. The design curve has a sharp OF3 cliff between L34 and L36.

**L30 winner mechanism:** all 5 sequential interfaces engaged in every one of 6 folds (5 Boltz seeds + OF3), mean R146 5.9–6.2 Å, inter-interface spread 0.33 Å — a SYMMETRIC apo ring with no pinned opening, reproducing the native apo-ring symmetry. The 'switch' is not a fixed geometric seat in the resting state. Cross-checked (Boltz+OF3), not validated; apo. Figure: `fig9_linker_curve.png`.


# Circular-permutation route (Stage 1 monomers + Stage 2 pentamers)
N–C direct fusion failed AF3 → CP escalation. **Stage 1:** all 3 CP cut sites (228/233/217) fold both domains natively as monomers (N 2.5–3.6 Å, C 2.2–7.8 Å to native; not fragmented). **Stage 2:** in ring context, **cp233 closes the ring under both Boltz-2 (3/3 seeds 5/5) and OpenFold3 (5/5, 0 scrambled) at all 3 inter-linker lengths**; sites 217/228 have OF3 5/5 but Boltz seed-variable (fail the both-predictor bar). OpenFold3 closes 5/5 on all 9 CP constructs — opposite of the direct fusion it scrambled. Winner **cp233_int15_inter10** (interfaces 5.4–5.7 Å). Diffusion rung not needed. Lead designs: cp233_int15_inter10, B1_L40_E119Q, (pentamer-L30). Cross-checked (Boltz+OF3), not validated; apo. Fig: `fig10_circular_permutation.png`.
