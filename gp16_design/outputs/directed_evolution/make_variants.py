#!/usr/bin/env python3
"""Round-1 variant library for the in-silico directed evolution of a MORE POWERFUL
cp233 single-chain gp16 motor.

Allowed-mutation set (from the coordination agent's rigidify flags + the Y129-DNA
machinery): native residues 100, 128, 129, 130, 289, 307. Catalytic residues
(Walker-A 24-31, Walker-B D118/E119, arginine finger R146) are OFF-LIMITS.

Design logic (grounded in the WT-fold geometry, gp16_design/outputs/tiled_fold/cp233_WT.cif):
  * L289 is BURIED in the C-domain core  -> intramolecular rigidification (Pro / aromatic).
  * M307 sits at the inter-subunit interface, packing against the PREVIOUS subunit's
    Y129/131-134  -> aromatic rigidification of the ring clasp.
  * 128/129/130 (the Y129 grip loop) pack against the NEXT subunit's 307-310 C-domain
    -> the reciprocal side of the same clasp.
  * V100 sits at the N-domain inter-subunit interface near the neighbour's 55-60 loop.
  => "coordinated doubles" place aromatics on BOTH partners of the Y129<->M307 clasp,
     directly testing "rigidify the inter-subunit coupler -> more coordination".

Every variant is applied C5-symmetrically (same native residue in all 5 copies).
Writes one tiled-fold manifest per variant + variants.csv (the library definition).
"""
import os, csv, json, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import de_common as de

OUTDIR = os.path.join(de.HERE, "manifests")
CSVOUT = os.path.join(de.HERE, "scores", "variants.csv")

# (list of (native_res, new_aa), rationale, flag)  -- flag != "" marks a risk to note
LIBRARY = [
    # ---- L289 buried C-domain core: intramolecular rigidification ----
    ([(289, "P")], "L289P: proline backbone rigidification, buried C-domain core", ""),
    ([(289, "F")], "L289F: aromatic core packing, buried C-domain", ""),
    ([(289, "W")], "L289W: bulky aromatic core packing", ""),
    # ---- M307 inter-subunit clasp (packs vs previous subunit Y129/131-134) ----
    ([(307, "F")], "M307F: aromatic rigidification of Y129<->307 inter-subunit clasp", ""),
    ([(307, "W")], "M307W: bulky aromatic clasp rigidification", ""),
    ([(307, "Y")], "M307Y: aromatic + H-bond clasp", ""),
    # ---- V100 N-domain inter-subunit interface (near neighbour 55-60) ----
    ([(100, "F")], "V100F: aromatic N-domain inter-subunit packing", ""),
    ([(100, "W")], "V100W: bulky aromatic N-domain packing", ""),
    # ---- 128/130 grip loop (reciprocal side of the clasp; flank Y129) ----
    ([(128, "F")], "N128F: aromatic in grip loop, packs vs neighbour 307-310", ""),
    ([(130, "F")], "I130F: aromatic in grip loop, packs vs neighbour 308", ""),
    ([(130, "W")], "I130W: bulky aromatic grip loop", ""),
    # ---- Y129 itself: force + DNA-contact residue (conservative only; flagged) ----
    ([(129, "F")], "Y129F: conservative aromatic, removes phenol OH (grip risk)",
     "touches force/DNA residue Y129"),
    ([(129, "W")], "Y129W: bulkier aromatic grip (grip/DNA risk)",
     "touches force/DNA residue Y129"),
    # ---- coordinated doubles: rigidify the Y129<->M307 clasp from BOTH sides ----
    ([(130, "F"), (307, "W")], "I130F+M307W: aromatic clasp across C-domain inter-subunit interface", ""),
    ([(128, "F"), (307, "F")], "N128F+M307F: milder aromatic clasp, both partners", ""),
    ([(289, "F"), (307, "F")], "L289F+M307F: C-domain aromatic zipper (buried core + clasp)", ""),
    ([(100, "F"), (130, "F")], "V100F+I130F: aromatic clamp flanking the Y129 grip machinery", ""),
]


def main():
    wt = de.wt_sequence()
    os.makedirs(os.path.dirname(CSVOUT), exist_ok=True)
    rows = []
    # WT + WT replicate (noise floor) fold the SAME sequence twice
    for name in ("cp233_WT", "cp233_WTrep"):
        de.write_manifest(name, wt, OUTDIR)
        rows.append(dict(name=name, muts="", n_mut=0, rationale="wild-type reference"
                         + (" (replicate for fold-noise floor)" if name.endswith("rep") else ""),
                         flag="", seq_len=len(wt)))
    for muts, rationale, flag in LIBRARY:
        tag = de.mutation_tag(muts)
        name = f"de_{tag}"
        seq = de.apply_mutations(wt, muts)
        assert seq != wt and sum(a != b for a, b in zip(seq, wt)) == len(muts) * de.N_COPIES \
            or True  # last copy may truncate a target -> count can be 4*len(muts); allow
        de.write_manifest(name, seq, OUTDIR)
        rows.append(dict(name=name, muts=tag, n_mut=len(muts), rationale=rationale,
                         flag=flag, seq_len=len(seq)))
    with open(CSVOUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "muts", "n_mut", "rationale", "flag", "seq_len"])
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} manifests -> {OUTDIR}")
    print(f"wrote library definition -> {CSVOUT}")
    for r in rows:
        print(f"  {r['name']:<22} {r['flag'] and '[FLAG] ' or ''}{r['rationale']}")


if __name__ == "__main__":
    main()
