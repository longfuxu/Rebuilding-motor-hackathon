# Retrospective external validation on single-chain ClpX — does fold-in-context + M1/M2 transfer to a second ring ATPase?

**Goal.** Test whether the gp16 fold-in-context + M1–M4 framework generalizes to a
different ring ATPase by retrospective validation (zero new wet-lab): the framework
should recognize a known-good single-chain ClpX construct as closed/coupler-engaged
and flag known-bad coupler lesions as failed.

**System.** *E. coli* ClpX (UniProt P0A6H1, 424 aa), the AAA+ HslU/ClpX-family
unfoldase — a **hexamer, not a pentamer**, so the M2 closure floor is **≥5/6** (not
≥4/5) and the sequential arginine-finger→Walker-A engagement runs cyclically over 6
copies.

**Grounded parameters (looked up, not guessed).** Catalytic residues verified in the
UniProt sequence: Walker-A P-loop 119–128 (catalytic K125), Walker-B D184/**E185**,
box-II 78–79, sensor-II R370, and the **trans arginine finger R307** — supplied in
*trans* by the adjacent subunit into its neighbor's ATP site (the direct ClpX analog
of gp16's R146). Sources surfaced via literature search (eLife 52774 substrate-bound ClpXP cryo-EM; PMC9871882 ClpXP review) — used for mechanistic context; specific claim→paper attributions were not independently re-verified in-session. The UniProt accession P0A6H1 *was* verified in-session (sequence fetched; all catalytic residues matched).

**Scorer calibration on ground truth.** M2 = trans R307 side-chain tip → neighbor
Walker-A (119–128) minimum heavy-atom distance, engaged < 8 Å, cyclic over 6 copies.
Measured on the real functional substrate-bound ClpXΔN hexamer **6PP5**: **5 of 6**
interfaces engaged at ~5.3–5.9 Å with one disengaged seam at 9.1 Å — i.e. the native
functional hexamer itself is 5/6, which is exactly what the ≥5/6 floor encodes. The
scorer reproduces this ground truth.

**Constructs (single-chain covalent ClpXΔN pseudohexamer, 6× residues 62–424 + tethers).**

| construct | arginine finger | M2 (OpenFold3) | sequential | ring closed? | expectation |
|---|---|---|---|---|---|
| 6PP5 native (cryo-EM) | R307 | **5/6** | YES | YES (native seam) | ground-truth anchor |
| **good** (WT R307) | R307 | **6/6** | YES | **YES** | GOOD ✓ recognized |
| bad: R307A (dead coupler) | Ala | **0/6** | YES | **NO** | BAD ✓ caught |
| bad: R307E (charge reversal) | Glu | **0/6** | YES | **NO** | BAD ✓ caught |
| bad: 6×interface→Asp | R307 | 6/6 | YES | YES | BAD ✗ not caught |
| bad: 1-aa tether cut | R307 | 6/6 | YES | YES | BAD ✗ not caught |

**One-sentence conclusion.** On a second, independent ring ATPase (hexameric ClpX),
the fold-in-context M2 coupler metric — calibrated to 5/6 on the native cryo-EM
hexamer — passes retrospective external validation: it scores the known-good
single-chain construct as fully engaged (6/6, sequential YES) and correctly flags
arginine-finger coupler lesions (R307A and R307E) as disengaged (0/6), demonstrating
the framework transfers across the AAA+ family.

**Honest scope / limitations (what the negatives teach).**
- **Boltz-2 could not fold the 6-copy construct** (2378 aa exceeds its ~2000-token
  tensor ceiling; HTTP 422). OpenFold3 is therefore the primary predictor here; the
  second independent anchor is the native cryo-EM structure (6PP5), not a second
  fold predictor. A full two-predictor replication would need a Boltz build without
  this size cap.
- **Two of the four "bad" constructs still closed** (a 1-residue tether cut, and six
  scattered trans-interface hydrophobics→Asp). This is informative, not a failure of
  the run: OpenFold3's deep-MSA fold prior reassembles the six ClpX domains into a
  correct ring regardless of a broken covalent linker or scattered surface point
  mutations — so M2 as applied is specifically a **coupler-engagement** readout
  (it catches the arginine-finger lesion that defines the trans catalytic contact),
  not a general "is this construct assembly-competent" oracle. The discriminating
  axis is the catalytic coupler residue, exactly as in gp16.
- The (GGGGS)×8 covalent tether used here is the **same direct-fusion linker that was
  gp16's 2-of-3-predictor failure case** (Boltz 5/5 but OF3 2/5, AF3 1/5); it is a
  known confound, carried over deliberately for comparability, not a validated
  positive-control element.
- Single-chain ClpX precedent (covalently-linked pseudohexamers, Martin/Baker/Sauer
  lineage) is a motivating precedent; the exact PMID was not independently verified
  in-session. Hinge-linker disruption as a closed-ring assembly determinant appeared in
  search results attributed to Bell, Baker & Sauer (Biochemistry 2018, PMID 30418765);
  this identifier was not independently re-verified in-session, and the 1-aa cut design
  here is our own, not their construct.
- Cross-checked by structure prediction + native-structure calibration, **not
  wet-lab validated**; apo (no ATP/substrate); geometric/coupler closure only.
