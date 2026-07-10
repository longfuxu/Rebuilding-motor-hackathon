# Does "single-chain > native coordination" generalize beyond gp16? — retrospective test on ClpX

**Question.** For gp16 cp233 we found the covalent single-chain ring is **~1.79×** more
inter-subunit coordinated than the native non-covalent ring by ENM perturbation-response
scanning (PRS), with Y129 a robustly more central force hub
(`gp16_design/outputs/coordination/`). Is this a **general principle of ring AAA+
ATPases**, or is it **gp16-specific**? We test it on a second, independent ring ATPase:
*E. coli* **ClpX** (hexameric AAA+ unfoldase), single-chain covalent design vs native.

**Headline answer.** **It does NOT generalize.** For ClpX the single-chain is at
**parity** with native by the primary proxy (PRS ATP→neighbour-ATP coupling design/native
= **0.91–1.08×**, cutoff-sensitive, straddling 1.0), and native is *higher* by two of the
other proxies (whole-ring inter/intra PRS 0.76×; arg-finger R307 betweenness 0.86×). The
large, robust single-chain advantage seen for gp16 (1.79×, one-directional across methods)
is **absent for ClpX**. This is a publishable negative/boundary result: the coordination
gain appears **gp16-specific**, not a universal consequence of covalent tethering.

---

## Setup (matched to the gp16 protocol, apples-to-apples)

- **Design (single-chain ClpX):** OpenFold3 prediction of the "good" covalent
  ClpXΔN pseudohexamer — 6× native residues 62–424 joined by 5× (GGGGS)×8 tethers, one
  chain, WT trans arginine finger R307, scored 6/6 coupler-engaged.
  `gp16_design/cs_sessions/session2_denovo_clpxvalidation/clpx_out/clpx_good__of3.cif`
- **Native ClpX:** cryo-EM hexamer **6PP5** (substrate-engaged spiral, 5×ATPγS + 1×ADP
  seam subunit F). `gp16_design/cs_sessions/session3_clpx_vs_gp16_topology/6PP5.cif`
- **Fair reduction:** both rings reduced to the **same residue set** — native residues
  present in all 6 native chains **and** all 6 design copies = **320 residues × 6 subunits
  = 1920 Cα** — with identical ANM/GNM parameters. Design construct→native mapping verified
  by aligning the Walker-A P-loop motif GPTGSGKT (`native_res = construct_pos − span_lo + 62`).
- **Landmarks (native ClpX / UniProt P0A6H1 numbering):** ATP P-loop / Walker-A **119–128**
  (catalytic K125); **trans arginine finger R307** (supplied in *trans* into the neighbour
  ATP site — the direct analog of gp16's R146); sensor-II R370; pore-1 loop aromatic
  **Y153** (GYVG substrate paddle — the analog of gp16's force residue Y129); substrate-grip
  set from 6PP5 chain-G proximity.
- **Ring topology recovered:** design = a clean **6/6** arg-finger cycle
  (0→5→4→3→2→1→0, all interfaces ~10.5 Å centroid, engaged); native = **5 engaged**
  interfaces (7.8–9.9 Å) **+ 1 disengaged ADP seam** (A↔F, 15.7 Å) — i.e. the real native
  5/6 spiral. The single chain closes the native seam.

Scripts: `clpx_coord_common.py` (mapping/landmarks), `clpx_prs_anm.py` (PRS/GNM),
`clpx_force_network.py` (betweenness), `clpx_make_figure.py`. Raw JSON alongside.

> **NO MD available for ClpX.** All numbers are ENM/PRS + GNM (elastic-network) only. The
> gp16 headline 1.79× was also an ENM-PRS number, so that axis is directly comparable; but
> gp16 additionally had a short-MD cross-check that ClpX lacks.

---

## (a) Coordination — ENM perturbation-response scanning (PRS)  ★ primary

ANM (cutoff 15 Å, all modes) → PRS; coordination = force-perturbation at one subunit's ATP
site (Walker-A) propagating to the **ring-neighbour's ATP site**, averaged around the
6-ring, normalised by the intra-subunit self-response. GNM (cutoff 10 Å) cross-correlation
for passive co-fluctuation.

| metric | design | native | design/native | gp16 design/native |
|---|---|---|---|---|
| **PRS ATP→neighbour-ATP coupling** ★ | 0.004975 | 0.004588 | **1.08×** | **1.79×** |
| &nbsp;&nbsp;↳ robustness across ANM cutoff | — | — | **0.91× (12 Å), 1.08× (15/18 Å)** | 1.79× (robust) |
| PRS ATP→neighbour, arg-finger direction | 0.004909 | 0.004579 | 1.07× | — |
| PRS ATP→ATP, ring-dist 2 | 0.003258 | 0.003091 | 1.05× | 2.28× |
| **PRS inter/intra (all residues, whole ring)** | 0.380 | 0.501 | **0.76× (native higher)** | 1.03× |
| GNM ATP↔neighbour-ATP cross-corr (passive) | 0.0364 | 0.0398 | 0.91× (native higher) | 0.67× |
| PRS ATP intra self-response | 0.1397 | 0.1428 | 0.98× | ≈1 |

**Read.** The primary proxy sits **at parity** (1.08× at 15 Å, but **0.91× at a tighter
12 Å cutoff** — the sign of the effect flips with a modelling choice, i.e. it is *within
noise of 1.0*). This is qualitatively different from gp16, whose 1.79× held one-directional
across cutoffs and methods. Worse for the "single-chain wins" hypothesis: the **whole-ring
inter/intra PRS coupling is clearly higher in native ClpX** (0.50 vs 0.38) — the native
substrate-engaged spiral is, if anything, the *more* inter-subunit-coupled object.

## (b) Force-transmission network (GNM-correlation network, no MD)

Residue network: nodes = Cα (1920, node-matched), edges = contacts (<9 Å), edge distance
= −log|C_ij| with C from GNM cross-correlation (co-moving contacts = "short" conduits).
Betweenness centrality = how central a residue is to force routing.

| node | design btw | native btw | design/native | gp16 analog (GNM) |
|---|---|---|---|---|
| **R307 trans arginine finger** (= gp16 R146) | 0.00327 | 0.00382 | **0.86× (native higher)** | — |
| Walker-A P-loop | 0.00591 | 0.00731 | 0.81× (native higher) | — |
| pore-grip 155/157/202 | 0.00212 | 0.00666 | 0.32× (native higher) | — |
| **Y153 pore-1 aromatic** (= gp16 Y129) †| 0.00292 | 0.00225 | **1.30×** (design higher) † | Y129 **1.65×** |
| intra conduit ATP→pore path length | 4.24 | 4.43 | design slightly shorter | — |

† *Y153 is disordered in the native ADP seam (F), so it cannot enter the node-matched set;
its betweenness is from each structure's full modeled Cα (design 6 subunits, native 5).
Native has fewer nodes there, which inflates its per-node betweenness — so design being
higher despite that handicap is the one design-favouring signal, but it is **not
node-matched** and should be read as suggestive only.*

**Read.** The one gp16-style "force hub" result (Y153 more central in the design, 1.30×)
survives weakly and only in the non-node-matched comparison; the **node-matched** trans
arginine finger R307 — the cleanest ClpX force-coupling residue — is **more** central in
**native** (0.86×), the opposite of gp16's Y129. So the "single-chain routes more force
through the functional hub" story does **not** robustly transfer.

---

## Verdict

**The gp16 finding does not generalize to ClpX — it is (as far as this test shows)
gp16-specific.** Against the available native structure, the covalent single-chain ClpX is
**at parity** with native in inter-subunit coordination (primary PRS proxy 0.91–1.08×,
straddling 1.0), and native is *higher* on the whole-ring inter/intra coupling (0.76×) and
on arg-finger centrality (0.86×). There is no ClpX analog of gp16's robust, one-directional
1.79× advantage.

**Interpretation — why gp16 ≠ ClpX (candidate causes, ranked):**
1. **Comparison/state mismatch (biggest confound).** gp16 was **apo design vs apo native**
   — matched, relaxed states. ClpX is **apo OF3 design vs a substrate-engaged native
   spiral** (6PP5): the native ClpX staircase is a substrate-locked, intimately stacked,
   functionally *coupled* state, which raises native coordination. A cleaner test needs an
   **apo native ClpX hexamer**. Until then the honest claim is conditional: *"vs the
   substrate-engaged native hexamer, the single-chain shows no coordination advantage."*
2. **Topology.** gp16's winner was an **optimised, MD-validated circular permutant (cp233)**
   whose new termini were engineered to reinforce the ring; the ClpX construct is a **direct
   head-to-tail fusion** with the (GGGGS)×8 tether that was gp16's own 2-of-3-predictor
   *failure* linker. The single-chain advantage in gp16 may come from the CP topology, not
   from covalency per se.
3. **Native seam.** 6PP5 carries a disengaged ADP seam (F) that *reduces* native coupling
   at one interface — yet native still matches/exceeds the design, which *strengthens* the
   conclusion that native ClpX is well-coordinated.
4. **No MD / ENM coarseness / apo OF3 design not MD-validated.**

**Both answers were publishable (general principle vs gp16-specific); the data say
gp16-specific.** The scientifically useful consequence: it **sharpens the gp16 claim** —
the coordination gain should be attributed to gp16's optimised circular-permutation
topology and/or the apo-vs-apo comparison, **not** advertised as a universal "covalent
single chains are more coordinated" law. The discriminating follow-up is a matched
**apo native ClpX vs the single-chain (ideally with MD)**; the ground truth remains the
single-molecule mixed-ring coordination experiment.
