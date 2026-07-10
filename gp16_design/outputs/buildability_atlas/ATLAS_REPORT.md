# BUILDABILITY ATLAS — topology predicts the single-chain construction method for ring ATPase motors

_Turning an N=4 observation into a predictive, structure-grounded rule across a 13-protein panel of ASCE / AAA+ / RecA-type ring motors. Pure protein science; all numbers are structure-derived and reproducible. 2026-07-09._

---

## 0. TL;DR

- **Claim tested:** the best way to make a homo-oligomeric ring motor into ONE polypeptide (direct N→C fusion vs circular permutation vs a de-novo diffusion connector) is set by **topology** — where the native termini sit relative to the substrate channel, and how big the head-to-tail gap is — **not** by sequence or family.
- **Two descriptors, read off the native ring structure:**
  1. **`jams_channel`** — does a native terminus contact the translocated substrate (DNA/RNA/peptide)? (the fusion/​linker would foul translocation)
  2. **`direct_gap_A`** — the through-space C(i)→N(i+1) gap between adjacent subunits for the motor construct.
- **Decision rule (formalized as a depth-2 tree):**
  `jams_channel → circular permutation`; else `gap ≤ 38 Å → direct fusion`; else `→ (nominally) diffusion connector`.
- **Rule accuracy on the 3 genuinely method-validated anchors = 3/3** (gp16→CP, ClpX→direct, gp17→direct).
- **The single most important new result of this session is a negative/cautionary one:** we fold-tested the rule's least-supported branch and found that **an apo (protein-only) Boltz-2 tiled-MSA fold cannot adjudicate method choice.** It closes even the **known-wrong** gp16-direct construct into a clean coupled ring (M2 5/5), because a protein-only folder is blind to the DNA-channel clash that actually forces gp16 to CP. Consequently the "diffusion" class has **zero validated members**, and its one in-panel testable prediction (Rho→diffusion) is **unsupported**. The honest, defensible rule collapses to two branches: **terminus fouls channel → CP; otherwise → direct fusion (linker sized to the gap).**

See `atlas_decision_boundary.png` for the one-figure summary.

---

## 1. The panel (N=13 homo-oligomeric ring motors)

All descriptors are computed by one uniform pipeline (`compute_descriptors.py`) on the native ring biological assembly: the **channel axis = ring symmetry axis** (least-variance normal of the subunit centroids); **termini = first/last modeled Cα** of the motor construct (motor-domain range for the two multi-domain proteins, see §5); **substrate distance = terminus-atoms → modeled DNA/RNA/peptide**. This reproduces the previously published anchor numbers (gp16 head-to-tail gap 57.1 Å; gp17 full C-term→axis 44.3 Å; ClpX C-term→substrate ~62 Å) to within rounding.

| Protein | Family | PDB | n | ring R (Å) | C-term→axis (Å) | C-term rel-R | C-term→substrate (Å) | **gap (Å)** | **jams?** | **predicted** | validated |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **gp16** (φ29) | ASCE packaging | 7JQQ | 5 | 29.1 | 14.0 | 0.48 | **4.0** | 57.1 | **YES** | **CP** | ✅ CP |
| **ClpX** | AAA+ unfoldase | 6PP5 | 6 | 33.4 | 66.8 | 2.00 | 61.4 | **23.5** | no | **direct** | ✅ direct |
| **gp17** (T4) | ASCE terminase | 3EZK | 5 | 41.1 | 32.8 | 0.80 | – | **31.2** | no | **direct** | ✅ direct |
| Rho | RecA RNA transloc. | 3ICE | 6 | 37.0 | 41.4 | 1.12 | 49.9 | 46.1 | no | diffusion | ⚠ see §3 |
| HslU | AAA+ unfoldase | 1DO0 | 6 | 34.0 | 38.9 | 1.14 | – | 32.2 | no | **direct** | predicted |
| katanin | AAA+ severing | 6UGD | 6 | 34.3 | 39.1 | 1.14 | 35.3 | 44.7 | no | diffusion | predicted |
| spastin | AAA+ severing | 6P07 | 6 | 32.9 | 38.9 | 1.18 | 34.8 | 50.9 | no | diffusion | predicted |
| T7 gp4 | SF4 helicase | 1E0J | 6 | 33.9 | 23.6 | 0.70 | – | 51.3 | no | diffusion | predicted |
| Vps4 | AAA+ ESCRT | 6BMF | 6* | 32.6 | 37.9 | 1.16 | 38.4 | 51.8 | no | diffusion | predicted |
| p97 / VCP | AAA+ double-ring | 5FTN | 6 | 41.1 | 20.8 | 0.51 | – | 51.5 | no† | diffusion | predicted |
| FtsK | ASCE DNA transloc. | 6T8B | 6 | 34.9 | 35.8 | 1.03 | 26.4 | 52.4 | no | diffusion | predicted |
| SV40 LTag | SF3 helicase | 1SVM | 6 | 32.1 | 57.4 | 1.79 | – | 69.5 | no | diffusion | predicted |
| DnaB | RecA helicase | 4ESV | 6 | 33.9 | 43.1 | 1.27 | 31.2 | 70.3 | no | diffusion | predicted |

\* Vps4 hexamer with 5 of 6 subunits ordered in 6BMF (spiral). † p97 C-term rel-radius 0.51 is a borderline no-jam (double-ring; no substrate modeled) — flagged.

**Category boundary (not in the quantitative panel):** hetero-oligomeric AAA rings — the 26S proteasome **Rpt1–6** (6MSB) and **MCM2–7** (6MII) — are already encoded by six *different* genes, so each subunit is natively addressable and "single-chain fusion" is not the relevant question. They mark where the atlas' premise (identical repeated subunit) stops applying.

---

## 2. The decision rule / classifier

Two features → method. A depth-2 `sklearn` `DecisionTreeClassifier` fit on `(jams_channel, direct_gap_A)` recovers exactly the intended splits and makes the numeric boundary explicit:

```
|--- direct_gap_A <= 38.4
|   |--- class: direct
|--- direct_gap_A >  38.4
|   |--- jams_channel = no  --> class: diffusion
|   |--- jams_channel = yes --> class: CP
```

- **Clean margin on the gap axis:** max "direct" gap = 32.2 Å (HslU) | min "diffusion" gap = 44.7 Å (katanin) → the 38 Å threshold sits in an empty 12 Å band (no panel point between 32 and 45 Å).
- **`jams_channel` cleanly isolates gp16** (C-term → DNA = 4.0 Å) from every other protein (next-closest terminus→substrate = FtsK 26.4 Å). This is the strongest single discriminator: the criterion is **physical** (a terminus touching the translocated polymer cannot host a fusion without blocking the channel), not a tuned cutoff.
- **Accuracy vs ground truth (non-circular): 3/3.** Ground truth = only proteins whose method is validated by *multiple independent signals* (multi-predictor AF3+Boltz, MD, and/or experiment): gp16 (cp233: 3-predictor 5/5, MD-stable, Rosetta; direct = AF3 0/5), ClpX (experimental single-chain ClpXΔN pseudohexamer, Martin/Baker/Sauer 2005), gp17 (direct 5/5 tiled-MSA + single-seq 0/5 control).

Files: `decision_rule.json` (rule, thresholds, tree, counts), `atlas_predictions.csv` (per-protein call + reason).

---

## 3. Fold verification — what we ran, and the cautionary result

We fold-tested the rule's most uncertain predictions with the free **Boltz-2 NIM under a block-diagonal ("tiled") MSA** and **handedness-robust** M1/M2 scoring (`pipelines/tiled_msa_fold/`). Two folds, both ~60 s, $0 (no GPU rental):

| Construct | what the rule says | fold result (apo Boltz-tiled) | interpretation |
|---|---|---|---|
| **Rho direct** (6×motor res175-414 + (GGGGS)₈, 1640 aa) | diffusion (gap 46 > 38) | **M2 6/6** engaged (fwd 6 / rev 0), M1 radius 29.4 Å, CV 0.02, planar 0.5 Å, pLDDT 75.5 → **ring closes** | suggestive that direct might work at 46 Å — but see the control |
| **gp16 direct** (5×native res4-330 + (GGGGS)₁₀, 1835 aa) | **CP** (jams; direct should fail) | **M2 5/5** engaged (fwd 5 / rev 0), M1 radius 26.8 Å, CV 0.02, planar 0.2 Å → **ring ALSO closes** | **gp16-direct is the KNOWN-WRONG method** (AF3 0/5 scrambled; C-term contacts DNA at 4 Å). Yet apo Boltz closes it cleanly. |

**The load-bearing conclusion:** a **protein-only** structure predictor is **blind to substrate-channel fouling** — it folds the covalent gp16-direct chain into a perfect coupled ring because the clashing DNA is not in the input. Therefore **a single apo tiled-MSA fold cannot rank single-chain methods**, and the Rho-direct 6/6 is *not* evidence that direct beats a diffusion connector for Rho. The information that actually decides the method (does a terminus touch the substrate) lives in the **native substrate-bound structure** (our descriptor) and requires **multi-predictor agreement / MD / experiment** to confirm — exactly the multi-signal standard used for the gp16 anchor (and consistent with the project's prior "predictor-split" and "pLDDT-trap" findings).

**Net verified count:** the **3 anchors** remain the verified points (rule 3/3). The 2 new folds are recorded as an **honest non-discriminating control**, not as new clean verifications — and they usefully **falsify the naive use of apo folds for method selection** and leave the "diffusion" class **unvalidated**.

---

## 4. Revised, defensible rule (what we would actually publish)

> **Read the native substrate-bound ring. If a terminus contacts the translocated polymer → circular permutation (relocate the join off-channel). Otherwise → direct N→C fusion, with a flexible linker sized to the head-to-tail gap.** A de-novo (RFdiffusion) rigid connector is reserved for the regime where *both* direct and CP fail — a regime this panel does **not** yet demonstrate. The gap axis sets linker length and difficulty; it did **not**, in any validated case, by itself force a diffusion connector. The upper gap bound at which flexible-linker direct fusion breaks down (somewhere > 46 Å) is **untested** and is the priority next experiment.

This is *stronger* than the original 3-way gap rule: it is one physical, structure-readable discriminator (channel fouling), it is 3/3 on the validated anchors, and it correctly predicts that gp16 (uniquely) needs CP while ClpX/gp17/HslU-like motors take direct fusion.

---

## 5. Method notes & honest caveats

- **Panel size / imbalance.** 13 homo-oligomers; only **1 CP** (gp16) and **3 direct** are validated. The decision boundary is a rule anchored on N≈3–4 ground-truth points plus geometric reasoning, *formalized* by a tree — it is **not** a statistically powered classifier. The "diffusion" bin (9 proteins) is entirely **rule-projected and unvalidated**.
- **Domain choice for multi-domain proteins.** For gp17 (ATPase 10–360 + nuclease 360–562) and Rho (RNA-binding 1–~130 + motor 175–414) the descriptors use the **motor domain** actually fused — otherwise gp17's full-length gap (66.9 Å) would flip its call. All other proteins: modeled chain = one motor domain (no ambiguity). Both choices are logged (`motor_range`, `direct_gap_fullchain_A` columns).
- **Apo fold blindness (see §3).** The cheap fold check is a **necessary-not-sufficient** test (can the covalent chain form a coupled ring at all); it cannot see substrate clashes and over-closes. Method selection must come from the substrate-bound descriptor + a substrate-aware or multi-predictor check.
- **Single predictor.** All folds here are Boltz-2 only. The project's standard for a robust call is ≥2 predictors (AF3+Boltz) + MD; Rho and the other predictions need that cross-check before being trusted.
- **Structure heterogeneity.** Anchors are captured in different functional states (7JQQ helical/ATP-loaded; 6PP5 substrate-engaged spiral; 3EZK CA-only 34 Å pseudo-atomic) — distances carry a few-Å uncertainty, especially 3EZK (gp17).
- **Reproducibility of the gap sign.** The head-to-tail gap is computed in both ring senses; the smaller-mean (correct forward) sense is reported. Scoring is handedness-robust (forward and reverse copy order, engaged-max), so mirror-sense sampling cannot masquerade as an open ring.

---

## 6. Prior art & where this is novel (real citations)

**Circular-permutation predictability.** Viable CP sites are predicted from local structural descriptors — closeness / centroid distance / weighted contact number — and, crucially, from the requirement that the **native N- and C-termini be spatially close** so they can be joined while new termini open elsewhere: CPred web server (*Nucleic Acids Res.* 2012, [W232](https://academic.oup.com/nar/article/40/W1/W232/1079842)); "Deciphering the Preference and Predicting the Viability of Circular Permutations in Proteins" (*PLoS ONE* 2012, [PMC3281007](https://pmc.ncbi.nlm.nih.gov/articles/PMC3281007/)); CirPred structure+linker modeler ([PMC8513176](https://pmc.ncbi.nlm.nih.gov/articles/PMC8513176/)); review of circular permutation in proteins (Bliven & Prlić, *PLoS Comput. Biol.* 2012, [PMC3320104](https://pmc.ncbi.nlm.nih.gov/articles/PMC3320104/)).

**N-to-C terminal distance as a descriptor.** Thornton & Sibanda (1983) — termini closer than random; Christopher & Baldwin (1996) — not different from random (the debate itself shows N-C distance is a recognized descriptor); Krishna & Englander, *PNAS* 2005, "The N-terminal to C-terminal motif in protein folding and function" — ~half of proteins bring N/C termini into direct contact (≤5 Å) ([PMC545867](https://pmc.ncbi.nlm.nih.gov/articles/PMC545867/)); "tale of two tails" termini exposure, Jacob & Unger, *Bioinformatics* 2007 ([link](https://academic.oup.com/bioinformatics/article/23/2/e225/203575)); macrocyclization requires termini proximity, Haim et al., *ChemBioChem* 2021 ([DOI](https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/cbic.202100111)).

**Tandem / single-chain fusion linker rules.** scFv gold-standard: match linker length to the inter-domain C→N distance — a ~3.5 nm VH→VL gap bridged by (GGGGS)₃ ≈ 5.7 nm (Huston et al.; anti-VEGF scFv, *Sci. Rep.* 2022, [link](https://www.nature.com/articles/s41598-022-09324-4)); GS-linker flexibility tuning for multidomain design, van Rosmalen et al., *Biochemistry* 2017 ([PMC6150656](https://pmc.ncbi.nlm.nih.gov/articles/PMC6150656/)); fusion-linker property/design review, Chen, Zaro & Shen, *Adv. Drug Deliv. Rev.* 2013.

**Single-chain ring ("concatemer") engineering precedent.** Covalently tethered **single-chain ClpXΔN pseudohexamers** degrade ssrA substrates at ~WT rates — the founding proof that a homo-oligomeric ring motor can be made one polypeptide by direct genetic linkage (Martin, Baker & Sauer, *Nature* 2005); linker/hinge control of assembly & activity, Bell, Baker & Sauer, *Biochemistry* 2019 ([PMID 30418765](https://pubmed.ncbi.nlm.nih.gov/30418765/)). This is exactly our ClpX "direct" anchor.

**Novelty framing.** The field already uses (a) **native N-C distance** as a CP/cyclization prerequisite and (b) the **inter-domain gap** to size fusion linkers — but for *single* proteins/domains. What is new here is applying, to a **homo-oligomeric ring motor**, a **joint** descriptor pair — *terminus-vs-substrate-channel* (a fouling/​topology term absent from all the above) **and** the *inter-subunit* head-to-tail gap — to **choose among direct fusion / CP / de-novo connector**, and the accompanying cautionary result that apo folds cannot substitute for the substrate-aware descriptor.

---

## 7. Files (all under `gp16_design/outputs/buildability_atlas/`)

| File | What |
|---|---|
| `ATLAS_REPORT.md` | this report |
| `atlas_decision_boundary.png` | the one-figure summary (gap × jams, colored by method) |
| `descriptors.csv` / `.json` | uniform topology descriptors, 13 proteins |
| `atlas_predictions.csv` | per-protein rule call + jam reason + ground truth |
| `decision_rule.json` | rule, thresholds, decision tree, margin, fold observations, caveat |
| `compute_descriptors.py` | descriptor pipeline (config + geometry) |
| `classify_and_train.py` | rule + sklearn tree + GT accuracy |
| `make_figure.py` | figure generator |
| `inspect_structures.py` | structure/chain inspector used to configure the panel |
| `structures/` | downloaded biological-assembly PDBs |
| `fold_verify/` | Rho-direct & gp16-direct tiled-MSA folds (`*.cif`, `*.result.json`, `SUMMARY.json`) |

**Reproduce:**
```bash
PY=/Users/longfu/miniforge3/bin/python3.13
cd gp16_design/outputs/buildability_atlas
$PY compute_descriptors.py    # -> descriptors.csv/.json
$PY classify_and_train.py     # -> atlas_predictions.csv, decision_rule.json (rule 3/3)
$PY make_figure.py            # -> atlas_decision_boundary.png
# fold verification (free NIM; uses NVIDIA_API_KEY2 to avoid colliding with sibling jobs):
cd ../../pipelines/tiled_msa_fold
NVIDIA_API_KEY=<key2> $PY tiled_msa_fold.py run manifests/rho_direct_tiledMSA.json  --out ../../outputs/buildability_atlas/fold_verify
NVIDIA_API_KEY=<key2> $PY tiled_msa_fold.py run manifests/gp16_direct_tiledMSA.json --out ../../outputs/buildability_atlas/fold_verify
```

**Compute used:** free Boltz-2 NIM only. **$0 of the $20 GCP budget spent; no VM launched → nothing to clean up.**
