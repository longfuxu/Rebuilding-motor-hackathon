#!/usr/bin/env python3
"""Round-2 variant library for the in-silico directed evolution of a MORE POWERFUL
cp233 single-chain gp16 motor.

Round-1 finding (see ../REPORT.md): the clean rigidify locus is the Y129<->M307 inter-
subunit clasp. Best round-1 variant = M307W (all gates, force-transmission proxy +3.5%
above noise, soft-mode intact). Coordination never robustly cleared its ~6% noise floor.

Round-2 hypotheses:
  (A) BIOCHEM-BLESSED rigidify AT the force hub. Ser127 sits immediately adjacent to the
      essential force residue Y129 (129 is near-lethal, cannot be mutated); Ser99 is a
      biochem-TOLERANT N-domain position. Both are in the phi29 biochem-TOLERANT set
      {6,53,56,99,127,128,222,233}. Aromatic/rigidifying substitutions here stiffen the
      force machinery WITHOUT touching a lethal residue:  S127F/S127W/S127Y, S99F/S99W.
  (B) STACK round-1 winners. Combine the best clasp aromatic (M307W / M307F) with a second
      biochem-safe rigidifier on the reciprocal grip loop (N128F), the adjacent force loop
      (S127F), or the buried C-domain core (L289P, a mild round-1 positive):
        M307W+N128F, M307W+S127F, M307F+S127F, M307W+L289P.

Every variant is applied C5-symmetrically (same native residue in all 5 copies). All targets
are in the biochem-TOLERANT set or already validated round-1 aromatic loci; none touch a
SEVERE/LETHAL residue (Y129,E58,K105,K328,R330,I28,Y32,S106,N126,N158,R234,K294) or a
catalytic residue (Walker-A/B, R146).

Writes one tiled-fold manifest per variant + variants.csv (the library definition).
"""
import os, csv, sys
HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, PARENT)
import de_common as de

OUTDIR = os.path.join(HERE, "manifests")
CSVOUT = os.path.join(HERE, "scores", "variants.csv")

# biochem gate sets (kept in sync with ../rank_variants.py)
BIOCHEM_SEVERE = {129, 58, 105, 328, 330, 28, 32, 106, 126, 158, 234, 294}
BIOCHEM_TOLERANT = {6, 53, 56, 99, 127, 128, 222, 233}

# (list of (native_res, new_aa), rationale, flag)
LIBRARY = [
    # ---- (A) biochem-blessed rigidify AT the force hub (adjacent to Y129 / N-domain) ----
    ([(127, "F")], "S127F: aromatic rigidification adjacent to force hub Y129 (biochem-tolerant)", ""),
    ([(127, "W")], "S127W: bulky aromatic adjacent to Y129 (biochem-tolerant)", ""),
    ([(127, "Y")], "S127Y: aromatic + H-bond adjacent to Y129 (biochem-tolerant)", ""),
    ([(99, "F")], "S99F: aromatic N-domain rigidification (biochem-tolerant)", ""),
    ([(99, "W")], "S99W: bulky aromatic N-domain (biochem-tolerant)", ""),
    # ---- (B) doubles stacking round-1 winners (M307W/F clasp + second safe rigidifier) ----
    ([(307, "W"), (128, "F")], "M307W+N128F: clasp aromatic + reciprocal grip-loop aromatic", ""),
    ([(307, "W"), (127, "F")], "M307W+S127F: clasp aromatic + force-hub-adjacent aromatic", ""),
    ([(307, "F"), (127, "F")], "M307F+S127F: milder clasp aromatic + force-hub-adjacent aromatic", ""),
    ([(307, "W"), (289, "P")], "M307W+L289P: clasp aromatic + buried C-domain proline rigidification", ""),
]


def main():
    wt = de.wt_sequence()
    os.makedirs(os.path.dirname(CSVOUT), exist_ok=True)
    os.makedirs(OUTDIR, exist_ok=True)
    rows = []
    # WT + WT replicate already reused from round 1 (folds/scores copied); still list them
    for name in ("cp233_WT", "cp233_WTrep"):
        rows.append(dict(name=name, muts="", n_mut=0,
                         rationale="wild-type reference"
                         + (" (replicate for fold-noise floor)" if name.endswith("rep") else ""),
                         flag="", seq_len=len(wt)))
    for muts, rationale, flag in LIBRARY:
        # safety: refuse SEVERE/LETHAL targets (belt-and-braces vs the rank gate)
        for r, _a in muts:
            assert r not in BIOCHEM_SEVERE, f"target {r} is biochem SEVERE -- refusing"
        tag = de.mutation_tag(muts)
        name = f"de_{tag}"
        seq = de.apply_mutations(wt, muts)
        assert seq != wt
        de.write_manifest(name, seq, OUTDIR)
        tol = sorted({r for r, _ in muts} & BIOCHEM_TOLERANT)
        rows.append(dict(name=name, muts=tag, n_mut=len(muts), rationale=rationale,
                         flag=flag, seq_len=len(seq)))
    with open(CSVOUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "muts", "n_mut", "rationale", "flag", "seq_len"])
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)-2} variant manifests (+2 references) -> {OUTDIR}")
    print(f"wrote library definition -> {CSVOUT}")
    for r in rows:
        print(f"  {r['name']:<24} {r['flag'] and '[FLAG] ' or ''}{r['rationale']}")


if __name__ == "__main__":
    main()
