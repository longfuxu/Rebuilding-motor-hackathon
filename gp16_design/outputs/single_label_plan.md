# Single-Label Fluorescence Plan — Position-Specific Labeling
_Generated 2026-07-07 · UniProt P11014 numbering_

## Rationale
Genetically identical subunits give mixed labeling stoichiometry and unknown label geometry —
fatal for FRET / single-fluorophore ring-opening assays. The ordered chain lets one label be placed
at one defined position.

## Label-site selection rule
Prefer solvent-exposed loops FAR from Walker A/B, R146, DNA channel, and pRNA interface — the same
exposed loops identified as CP candidates double as label sites. Introduce a unique Cys (or ybbR/SNAP tag)
at one position; keep all others label-free by construction.

## Candidate label positions (exposed, functionally clear loops)
- around res **248–249**, **284–286**, **169–171** — exposed, clear of ATP/DNA/pRNA (computed clearances all >18 Å to ATP).
- avoid res 280 and 297 themselves in a given subunit if used as that subunit's CP junction.

## Readouts
- photobleaching-step stoichiometry (confirms single label), packaging-correlated intensity/environment change,
  position-dependent FRET vs. a DNA/pRNA/capsid marker, compatibility with combined optical-tweezers + confocal.

## What each geometry reports
- A label near the **A→E special interface** reports ring opening / the broken interface.
- A label on a **DNA-proximal loop** reports translocation-coupled motion.
- Realizes the FRET-during-packaging experiment proposed as future work in the 2021 inchworm model.
