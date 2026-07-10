# AF3 Scoring Report — gp16 single-chain ring designs (2026-07-09)

AlphaFold3 is the **3rd independent structure predictor** applied to these designs
(Boltz-1 and OpenFold3 came first). AF3 was run with its own MSA. Judged strictly by
**M1** (compact, in-sequence-order ring) and **M2** (trans arginine-finger R146 → next
subunit's Walker-A P-loop < 8 Å, in sequential register). **M2 pass = best model ≥ 4/5
engaged interfaces** for the pentamer. Global pTM / pLDDT are *not* used as pass criteria
(high pTM on a scrambled ring is a known failure mode — steer by interface metrics).

Every one of the 9 folds was re-scored by an **independent adversarial verifier** (biotite
geometry recomputed from the .cif, *not* re-running `score_m2.py`). All 9 verifier verdicts
returned **CONFIRMED** with `matches_reported = true`.

## 1. Scoring table

| Fold | Type | Best M2 | M1 ring | M1 seq | Verify | iface pLDDT |
|---|---|:---:|:---:|:---:|:---:|:---:|
| cp233_novel_d5 | cp233_novel (de-novo) | **5/5** ✅ | yes | yes | CONFIRMED | 72.3 |
| cp233_novel_d2 | cp233_novel (de-novo) | **5/5** ✅ | yes | yes | CONFIRMED | 71.9 |
| cp233_novel_d3 | cp233_novel (de-novo) | **5/5** ✅ | yes | yes | CONFIRMED | 70.4 |
| cp233_novel_d7 | cp233_novel (de-novo) | 2/5 ❌ | yes | no | CONFIRMED | 71.3 |
| cp233_novel_d8 | cp233_novel (de-novo) | 2/5 ❌ | yes | no | CONFIRMED | 72.3 |
| cp285_int15_inter10 | cp285_CP | **5/5** ✅ | yes | yes | CONFIRMED | 68.1 |
| cp297_int15_inter10 | cp297_CP | **5/5** ✅ | yes | yes | CONFIRMED | 67.7 |
| rfdiff_ring_sal_l50_0_d0_1835aa | rfdiff_ring (ModeB) | 0/5 ❌ | yes | no | CONFIRMED | 81.1 |
| rfdiff_ring_sal_l50_0_d2_1835aa | rfdiff_ring (ModeB) | 2/5 ❌ | yes | no | CONFIRMED | 73.0 |

Per-model M2 (5 AF3 models each):
- cp233_novel_d5 `[5,5,5,5,5]` · d2 `[3,5,4,5,5]` · d3 `[3,5,5,5,4]` · d7 `[2,2,2,2,2]` · d8 `[2,1,2,2,1]`
- cp285 `[5,5,5,5,5]` · cp297 `[5,5,5,5,5]`
- rfdiff d0 `[0,0,0,0,0]` · rfdiff d2 `[1,2,2,2,2]`

**Note on M1 ring vs M1 seq:** every fold forms a geometrically compact, planar ring
(M1 ring = yes for all 9). The discriminating column is **M1 seq** (is the subunit
*register* the designed sequential order?). The three failing classes (d7, d8, both
rfdiff) close a ring but in a **scrambled register**, which is exactly why their M2 is low.

## 2. Per-group verdicts

### RFdiffusion Mode-B rings — do they close on AF3? **No.**
AF3 with its own MSA does **not** flip the RFdiffusion rings that were 0/5 under
single-sequence Boltz. Both fail M2:
- **d0: 0/5, no model near passing.** The ProteinMPNN sequence deleted the native
  arginine finger in 4 of 5 copies (copy 1 keeps ARG@143; copies 2–5 are ALA/ASN/VAL/SER),
  so a trans arg-finger contact is only *measurable* at one interface — and even that one
  sits at ~19 Å. Verifier confirms 0.
- **d2: 2/5.** Compact planar pentamer, but only the designed A1–A2 / A2–A3 pairs engage
  (~6.9–7.0 Å); A3/A4/A5 R146 point at non-adjacent subunits (34–54 Å). Scrambled register.
- **iface pLDDT is a trap here:** d0 has the *highest* interface pLDDT of the whole set
  (81.1) while being the *worst* M2 (0/5). Confirms the "don't judge by pLDDT/pTM" rule —
  a confidently-predicted ring can be confidently in the wrong register.

### cp285 / cp297 secondary CP sites — do they close like cp233? **Yes.**
Both close cleanly on AF3:
- **cp285: 5/5 in all 5 models**, compact sequential ring, verifier CONFIRMED.
- **cp297: 5/5 in all 5 models**, compact sequential ring, verifier CONFIRMED.
- Caveat: their **interface pLDDT (~67–68) is lower** than the cp233 family (~70–72), and
  the cp297 M2 contacts cluster **marginally** just under the cutoff (all five in a narrow
  7.4–7.7 Å band). Engagement is real and directionally correct (reverse-neighbour control
  ~54 Å) but the interfaces are looser than cp233_novel_d5's. Net: the closable-CP-site set
  extends beyond cp233 to at least cp285 and cp297.

### cp233_novel de-novo re-sequences — which are the best closers? **d5 > d2 ≈ d3; d7/d8 fail.**
5 de-novo sequences (~53% identity) on the cp233 scaffold; **3 of 5 pass M2, 2 fail.**
- **d5 — cleanest closer.** The only design that is **5/5 in all five AF3 models**
  (`[5,5,5,5,5]`), highest iface pLDDT of the passers (72.3), radius CV ~0.01. Best pick.
- **d2 — pass, robust.** Best model 5/5, three models at 5/5 (`[3,5,4,5,5]`), iface 71.9.
- **d3 — pass, robust.** Best model 5/5, three models at 5/5 (`[3,5,5,5,4]`), iface 70.4;
  best model near-perfect (radius_CV 0.01, planarity 0.2 Å).
- **d7 — fail.** Uniform 2/5 across all models; compact ring, scrambled register.
- **d8 — fail.** Best 2/5 (`[2,1,2,2,1]`); compact ring, scrambled register.

## 3. Discrepancies flagged by the adversarial verifier

**No result-level discrepancies:** all 9 verifier verdicts are CONFIRMED and match the
reported M2 counts. Two methodological items are worth recording because they *could* have
produced false negatives and will bite anyone who re-scores:

1. **rfdiff d2 — copy-span registration (near-miss false negative).** The copy spans handed
   to the scorer were a **uniform 367-aa split** (1835/5). That is **wrong** for copies 2–5.
   The true repeat is **377 aa** (327-aa gp16 domain + 50-aa linker; the last copy is
   truncated to 327, i.e. 377×4 + 327 = 1835). Evidence: the canonical Walker-A P-loop
   `GARGIGKS` recurs at residues 21/398/775/1152/1529 (spacing 377), and the R146 anchor at
   136/513/890/1267/1644 (spacing 377). Under the wrong 367-spacing the R146 landmark lands
   on ALA/ASN/VAL/SER and Walker-A on non-P-loop residues for copies 2–5, giving a **spurious
   0/5**. With structurally-correct spans the verifier reproduces **2/5**. *Meaning:* the
   `l50` RFdiffusion constructs need span-aware scoring; a naive equal-split mis-registers
   landmarks and under-counts. (Does not change the conclusion — the fold still fails.)

2. **rfdiff d0 — landmark validity caveat.** Only copy 1 carries a genuine arginine finger +
   Walker-A motif; MPNN mutated the finger away in copies 2–5, so M2 is only physically
   defined for one interface. That interface is still ~19 Å (not engaged). M2 = 0 stands.

**Tooling note:** `gp16_design/reproduce/score_m2.py` crashes (None-format, ~line 196) when
an interface is *not* engaged, because it tries to format a missing "nearest" distance. It
was worked around cosmetically / via a temp-patched scorer for the failing rfdiff folds; the
M1/M2/M4 numbers are unaffected. Fix the None-guard before the next batch so non-engaged
interfaces print instead of crashing.

## 4. Bottom line for the paper's cross-predictor table

- **The cp233 scaffold family closes on a 3rd independent predictor (AF3, own MSA):** the
  cp233 design plus **3 of 5 de-novo re-sequences (d2/d3/d5, ~53% id)** give sequential,
  compact **5/5 M2 rings**; **d5 is unanimous 5/5 across all 5 AF3 models** — the cleanest
  de-novo closer.
- **Two additional CP sites close:** cp285 and cp297 are **5/5 in every AF3 model**,
  extending the closable-circular-permutation set beyond cp233 (at modestly lower interface
  pLDDT, ~67–68).
- **RFdiffusion Mode-B de-novo rings do not close** on AF3 (0/5 and 2/5), consistent with
  single-sequence Boltz — they form geometrically compact rings (highest pLDDT of the set is
  the *worst* M2), but the subunit register is scrambled.

---
## RECONCILIATION (2026-07-09, after tiled_msa_fold results)
Two corrections to the RFdiffusion verdict above, after cross-checking with the Boltz-2 tiled-MSA folds:
1. **Copy span was wrong** here (uniform 367; true repeat is 377 = 327 domain + 50 linker). Re-scored with correct 377 landmarks + handedness-robust M2.
2. **The RFdiffusion ring is PREDICTOR-SPLIT, not simply "does not close":** the identical 1835aa sequence folds to **5/5 clean ring under Boltz-2 + tiled MSA** but **2/5 scrambled (radius_CV 0.09) under AF3** — a genuine predictor disagreement (confirmed with matched correct scoring), NOT a scoring artifact. Contrast cp233 (all 3 predictors agree 5/5). So the generative rung yields a marginal, non-robust design.
3. **Handedness caveat:** covalent-ring M2 must be read handedness-robust (both k→k+1 and k→k-1); mirror-wound folds otherwise misread as 0/5. The cp233_novel d7/d8 and cp285/cp297 verdicts here used direction-specific scoring — cp285/297 are corroborated 5/5 by Boltz-tiled too (robust); d7/d8 should be re-checked handedness-robust (may or may not flip).

### d7/d8 handedness re-check (closed 2026-07-09)
Re-scored cp233_novel d7/d8 handedness-robust: **still 2/5 in BOTH directions (fwd 2, rev 0)** → GENUINELY scrambled, NOT a mirror-winding artifact (contrast E119Q: 0 fwd / 5 rev = mirror-wound-but-closed). So the de-novo verdict is final and correct: **d5 (unanimous 5/5) > d2 ≈ d3 pass; d7, d8 genuinely fail.** The handedness-robust scorer cleanly separates "closed ring wound the other way" from "no coherent ring."
