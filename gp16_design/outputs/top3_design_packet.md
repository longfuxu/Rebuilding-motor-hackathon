# Top-3 Design Packet — Position-Addressable φ29 gp16 Ring
_Generated 2026-07-07 · source structure PDB 7JQQ (cryo-EM, 4.1 Å) · sequence UniProt P11014 (332 aa)_

## Design principle (key correction incorporated)
7JQQ is the **helical / lock-washer** state of the motor — ATP-loaded and mid-opening
(subunits A,B,C carry ATP-γ-S+Mg²⁺; D,E are apo; 9.5 Å axial staircase). The resting ring
is hypothesized **planar**. Therefore every inter-subunit linker is sized for the **helical
(open) maximum** with flexible excess, so it is *compliant* — not load-bearing — in the planar
state and does **not splint** the planar↔lock-washer stroke that ATP loading drives.

## Progressive build strategy (as directed)
Validate topology on the **2-subunit CP dimer first**, then extend 3 → 4 → 5. Each step reuses
the identical CP site and linker rule, so a success at n=2 de-risks the full ring.

## Lead circular-permutation sites
| CP site | Loop | Junction (helical max) | Stroke (planar↔helical) | DNA / pRNA / ATP clearance |
|---|---|---|---|---|
| **280** | 279–282 | 27.0 Å | **3.0 Å** (lowest) | 14.0 / 23.3 / 46.9 Å — clear |
| **297** | 297–298 | 24.7 Å (shortest) | 6.2 Å | 11.2 / 20.2 / 52.4 Å — clear |

Both cut in exposed loops far from Walker A (24–31), Walker B (115–119), the trans arginine
finger **R146**, the DNA channel, and the pRNA interface. Native head-to-tail fusion is rejected:
the native C-terminus (res330) is itself a DNA-channel contact (6.2 Å) and the A→E seam spans 66 Å.

## Top 3 synthesis candidates
1. **gp16_CP280_2mer** (719 aa) — first build. CP@280, internal GS₅ bridge (330→4), inter-subunit GS₃.
   Lowest stroke; validates that CP topology folds and the junction closes before scaling.
2. **gp16_CP280_3mer** (1086 aa) — second build. Adds one junction; tests cooperative closure of a partial arc.
3. **gp16_CP280_5mer** (1820 aa) — target ring. Full P1–P5 addressable single chain; synthesize only after n=2/3 pass.

CP@297 variants are the parallel backup series (shortest junction) should CP@280 fold poorly.

## Module architecture (per subunit, CP@280)
`[res 280→330] – GS₅ internal bridge – [res 4→279]`  → new N-terminus at 280, new C-terminus at 279.
Subunits joined by `GS₃` inter-subunit linkers. Native disordered res 1–3 omitted.

## Addressability
Chain order fixes ring positions P1→Pn unambiguously. Each position independently accepts one
point mutation (dead-subunit scan) or one label site (fluorescence), because the positions are now
genetically distinct sequence blocks rather than identical copies.
