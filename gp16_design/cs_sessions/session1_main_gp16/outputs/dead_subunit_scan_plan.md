# Dead-Subunit Scan Plan — Position-Specific Catalytic Perturbation
_Generated 2026-07-07 · residue numbers in UniProt P11014 numbering (matches 7JQQ)_

## Rationale
A homomeric ring cannot assign a catalytic-dead subunit to a defined position. In the ordered
single chain, one mutation is placed in exactly one of P1–P5, one position at a time.

## Perturbation menu (choose one per position)
| Target | Residue(s) | Effect | Basis |
|---|---|---|---|
| Walker A / P-loop | **K30A** | abolish ATP binding | Walker A = 24–31 (GARGIGKS), catalytic Lys30 |
| Walker B | **E119Q** | block hydrolysis, retain binding | Walker B = 115–119 (IVFDE), catalytic D118/E119 |
| Arginine finger (trans) | **R146A** | break inter-subunit catalytic coupling | R146 of neighbor reaches into the ATP pocket (computed, <6 Å to AGS) |

## Scan design (5 constructs on the CP@280 5-mer backbone)
P1-dead, P2-dead, P3-dead, P4-dead, P5-dead — identical except the position of a single E119Q (recommended
first: retains nucleotide binding, isolates the hydrolysis step).

## Predicted readouts (single-molecule optical tweezers)
- packaging velocity, dwell-time distribution, burst size, step size (if resolvable), pause/slip frequency,
  ATP-dependence, initiation probability.

## Interpretation logic (tests the special-subunit hypothesis)
- **One position gives timer/regulatory-like failure, others tolerate** → special subunit is *position-defined*.
- **All positions fail equivalently** → special role *rotates* or emerges dynamically.
- **Adjacent positions graded** → coordination propagates through nearest-neighbor interfaces.
- Because 7JQQ shows the ATP-block (A,B,C) contiguous and apo-block (D,E) contiguous, the informative
  comparison is a dead subunit placed *at the ATP/apo boundary* (the A→E special interface) vs. within a block.
