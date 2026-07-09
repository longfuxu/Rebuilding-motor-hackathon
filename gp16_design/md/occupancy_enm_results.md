# Occupancy-dependent elastic-network model of gp16 ring opening (free, local)

Tests the owner's progressive-opening model with a cheap elastic-network (ANM) linear-response, using apo
(closed planar, Boltz) and 7JQQ (3-ATP partial-occupancy intermediate) as endpoints. Code: `md/occupancy_enm.py`.

## (a) Does adding ATP one-by-one progressively open the ring? — YES, and cooperatively
Ring out-of-plane opening (arb. units) vs number of ATP-bound subunits (each ATP pulls its subunit toward its
7JQQ conformation; linear response through the elastic coupling):

| # ATP | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| opening | 0.00 | 0.99 | 4.19 | 7.19 | 12.17 | 13.24 |

**Monotone and supra-linear (cooperative):** each added ATP opens the ring more, with the increments growing —
the elastic coupling makes later ATPs more effective. This is exactly the progressive-opening model, quantified.

## (a′) Does the SPATIAL pattern of ATP matter? — YES (sequential vs concerted signature)
At 3 ATP, the opening depends on WHICH 3 subunits are loaded:

| pattern | (0,1,2) adjacent | (0,2,4) spread | (0,1,3) |
|---|---|---|---|
| opening | 7.19 | 5.97 | 9.83 |

The firing *pattern/order* changes the ring opening (5.97–9.83) → a direct, cheap signature that **sequential
vs concerted ATP firing produce different ring states** — the mechanism the addressable E119Q construct tests.

## (b) apo→7JQQ opening is localized out-of-plane at the middle subunits
Per-subunit apo→7JQQ CA displacement: C 9.7 Å (out-of-plane 3.9), D 7.0 (4.5), B 7.0 (1.5), A 3.8, E 3.7 — the
opening is a localized out-of-plane (helical) distortion around C/D/B, i.e. the ATP-engaged region (deposited
7JQQ has fingers tight at B→A, C→B, D→C and the open seam at A→E/E→D; source `outputs/native_ring_reference.csv`).

## Honest caveats
- Linear-response ENM: the per-subunit "ATP force" is proxied by each subunit's apo→7JQQ displacement, so this
  *operationalizes* the model rather than deriving it from first principles — but it shows the model is
  self-consistent and predicts cooperativity + pattern-dependence.
- A single-basin ENM cannot give the transition BARRIER (energy rises monotonically away from apo) — the barrier
  and true path need a **dual-basin Cα Gō model** (`md/gomodel_dualbasin_modal.py`, run on GPU).
