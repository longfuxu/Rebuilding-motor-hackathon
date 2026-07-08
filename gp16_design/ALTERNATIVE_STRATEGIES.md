# gp16 single-chain ring — design strategies & results record

**Goal (primary):** design a complete, buildable single-chain (or covalently-defined) phi29 gp16
ring in which one ring position is genetically addressable — the way single-chain ClpX places a
catalytically-dead subunit at a defined position (Martin, Baker & Sauer 2005, *Nature*, PMID 16237435,
DOI 10.1038/nature04031). Winning the hackathon is a by-product; the deliverable is the construct.

All results below are **computational and cross-checked by two predictors (Boltz-2 + Chai-1)** — not
wet-lab validated. Folds were made **apo** (no pRNA/DNA/ATP unless stated). 7JQQ is the **helical
(ATP-loaded, lock-washer)** state; the resting ring may be planar. These constructs address a
**geometric ring seat**, not a guaranteed fixed functional role.

---

## The design arc (why we are where we are)

| Cycle | What was tried | Result | Rule learned |
|---|---|---|---|
| 1 | Single-chain, circular permutation at res 280/297 | **Failed** — cut sits *inside* the C-domain and fragments a folding unit; trans-R146 lost | **R1: cut only at folding-unit boundaries** |
| 2 | Domain-aware cut-penalty filter (free); hinge-CP (228/233); native-order 2-mer | 2-mer interface reformation is **stochastic** (Boltz 0/5, Chai 7/10) | **R2: the gp16 interface is a ≥3-subunit property** — each C-domain contacts its *two neighbors'* C-domains, not its own N-domain. A 2-mer under-determines it. |
| 3 | **B1**: native-order covalent dimer + 3 WT = 5-chain ring | **Works** — reforms trans-R146 as reliably as the native ring | **R3: test in ring context; tether-don't-permute is the minimal design that satisfies R1–R3** |

### Cycle-3 headline (cross-checked)
| construct | Boltz-2 R146 engaged | Chai-1 R146 engaged | tether-flanking seats |
|---|---|---|---|
| native ring (reference) | 15/15 | 50/50 | — |
| **B1_L40 (winner)** | **15/15** | **47/50** | Boltz 6/6 · Chai 20/20 |
| B1_L32 | 15/15 | 45/50 | Boltz 6/6 · Chai 19/20 |
| isolated 2-mer (cycle 2) | 0/5 | 7/10 | under-determined |

Bonus observation: the **apo** ring folds *symmetric* (all 5 arginine fingers ~6.6 Å); 7JQQ (+ATP) is
*asymmetric* (tight at ATP interfaces, open at the seam) — the planar↔helical hypothesis rendered by
the fold. Needs the +ATP refold (Cycle 4) to quantify.

---

## Where B1 sits vs the primary goal

B1 (covalent **dimer** + 3 WT) is the **minimal position-addressable unit** and it is proven in ring
context. It is **not yet the complete single-chain ring.** The route to the full single chain is a
**progressive fusion ladder**, now testable because R3 gives the correct test context:

- dimer ✅ (Cycle 3) → **trimer → tetramer → pentamer**, native order, channel-clearing linkers,
  scored on whether **all 5 interfaces** reform in ring context (each internal subunit gains both
  neighbors as the chain closes).
- Open risk = the **seam**: the A→E "special/open" interface spans 66 Å in the helical state and the
  native C-terminus (res330) contacts the DNA channel — so the ring-closing linker is the hard one.

---

## Strategy menu (with verified prior art)

**B1 — covalent sub-oligomer embedded in a WT ring → full single chain.** *Primary route.*
Precedent that this works for ring ATPases: single-chain ClpX with dead subunits at defined positions
(Martin, Baker & Sauer 2005, *Nature*, PMID 16237435; 2008 *NSMB*, PMID 18931677). Single-chain
oligomer for asymmetric mutation / one-site labeling: HIV-1 protease tethered dimer (Cheng et al. 1990,
*PNAS*, PMID 2263618; Torbeev & Kent 2012, *Org Biomol Chem*, DOI 10.1039/c2ob25569c).

**B2 — folding-unit graph (method spine).** Nodes = folding units (N-dom, C-dom per subunit) with
inter-subunit contact weights from 7JQQ; a valid single-chain design = a path that never severs a unit
and joins only spatially-adjacent units; score = linker span + broken inter-subunit contacts. Explains
cycle 1–3 in one framework. Tolerance-at-boundaries basis: Coyote-Maestas et al. 2021, *Nat Commun*,
PMID 34880224; Atkinson et al. 2019, *PEDS*, PMID 32626892; CP design rules: Yu & Lutz 2011,
*Trends Biotechnol*, PMID 21087800.

**B3 — pRNA-scaffold addressing (phi29-unique alternative).** Each gp16 subunit contacts exactly one
pRNA chain (A↔K … E↔O). Engineer pRNA interlocking-loop complementarity to template a defined subunit
order without a giant single chain (Guo et al. 1998, *Mol Cell*, PMID 9734359; Zhang F et al. 1998,
*Mol Cell*, PMID 9702199).

**B4 — orthogonal-interface heterodimer register.** Charge-swapped subunits that pair only in a defined
order, doped into WT — a defined register without covalent fusion (β-clamp analog: Heltzel et al. 2009,
*PNAS*, PMID 19617571).

---

## Standing caveats (keep in every writeup)
- Ring cycles planar↔helical during translocation: Woodson et al. 2021, *Sci Adv*, PMID 33990327,
  DOI 10.1126/sciadv.abc1955.
- The "special/regulatory" subunit is a dynamically reassigned role, not a fixed seat: Tafoya et al.
  2018, *PNAS*, PMID 29987018, DOI 10.1073/pnas.1801709115; Chistol et al. 2012, *Cell*, PMID 23178120.
- The field currently perturbs subunits only by **stochastic doping** + optical-tweezers statistics
  (Moffitt et al. 2009, *Nature*, PMID 19092931) — the unmet need a position-defined construct fills.

**Citation provenance:** all PMIDs/DOIs above were verified against Europe PMC + Crossref during the
project. The Cofactor literature MCP was offline; verification used the public gateways directly.
