#!/usr/bin/env python3
"""GREEDY EPISTATIC round: fix the round-1 winner M307W as the BACKGROUND and scan a BROAD
set of allowed second sites X, building M307W+X for every X. Directly tests whether ANY
second mutation stacks on top of the best single (two mutations better than one), instead of
the 4 pre-chosen doubles of round 2 (under-sampled).

Second sites = biochem-TOLERANT set {56,99,127,128,222,233} + untested-rigidify interface/core
loci {100,130,289}, each mutated to the rigidifying residue(s) that make structural sense
(aromatics F/W/Y; Pro for the buried C-domain core L289). 307 is NOT re-mutated. All SEVERE/
LETHAL residues (129,58,105,328,330,28,32,106,126,158,234,294) are excluded, as are catalytic
residues (Walker-A/B, R146) via de_common.FORBIDDEN.

Every variant is C5-symmetric. 3 of these (M307W+N128F/S127F/L289P) were already folded+scored
in round2/ and are reused (their folds/scores were copied into round2_greedy/). Writes one
tiled-fold manifest per NEW variant + variants.csv (the full library incl. reused ones)."""
import os, csv, sys
HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, PARENT)
import de_common as de

OUTDIR = os.path.join(HERE, "manifests")
CSVOUT = os.path.join(HERE, "scores", "variants.csv")

BIOCHEM_SEVERE = {129, 58, 105, 328, 330, 28, 32, 106, 126, 158, 234, 294}
BIOCHEM_TOLERANT = {6, 53, 56, 99, 127, 128, 222, 233}

BG = (307, "W")   # fixed background = round-1 winner M307W

# second sites -> rigidifying substitutions (native residue : list of new AAs)
SECOND_SITES = {
    56:  ["F", "W"],          # K56  tolerant, N-domain
    99:  ["F", "W"],          # S99  tolerant, N-domain (S99W = round-2 force line)
    100: ["F", "W"],          # V100 untested, N-domain inter-subunit interface
    127: ["F", "W", "Y"],     # S127 tolerant, adjacent to force hub Y129
    128: ["F", "W"],          # N128 tolerant, grip loop (reciprocal clasp side)
    130: ["F", "W"],          # I130 untested, grip loop
    222: ["F", "W"],          # Q222 tolerant
    233: ["F", "W"],          # K233 tolerant
    289: ["P", "F", "W"],     # L289 untested, buried C-domain core (Pro rigidifier)
}


def main():
    wt = de.wt_sequence()
    os.makedirs(OUTDIR, exist_ok=True)
    os.makedirs(os.path.dirname(CSVOUT), exist_ok=True)
    rows = [dict(name="cp233_WT", muts="", n_mut=0, rationale="wild-type reference", flag="", seq_len=len(wt)),
            dict(name="cp233_WTrep", muts="", n_mut=0, rationale="wild-type replicate (noise floor)", flag="", seq_len=len(wt)),
            dict(name="de_M307W", muts="M307W", n_mut=1, rationale="round-1 winner = the fixed BACKGROUND", flag="background", seq_len=len(wt))]
    n_new = 0
    for res in sorted(SECOND_SITES):
        assert res not in BIOCHEM_SEVERE and res != 307
        for aa in SECOND_SITES[res]:
            muts = [BG, (res, aa)]
            tag = de.mutation_tag(muts)          # e.g. "M307W_S127F"
            name = f"de_{tag}"
            seq = de.apply_mutations(wt, muts)
            assert seq != wt
            de.write_manifest(name, seq, OUTDIR)  # (re)writes manifest; fold_all skips if folded
            n_new += 1
            tol = res in BIOCHEM_TOLERANT
            rows.append(dict(name=name, muts=tag, n_mut=2,
                             rationale=f"M307W background + {de.wt_identity(res)}{res}{aa} "
                                       f"({'biochem-tolerant' if tol else 'untested-rigidify'} second site)",
                             flag="", seq_len=len(seq)))
    with open(CSVOUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "muts", "n_mut", "rationale", "flag", "seq_len"])
        w.writeheader(); w.writerows(rows)
    print(f"wrote {n_new} M307W+X manifests -> {OUTDIR}")
    print(f"library ({len(rows)} rows incl. WT/WTrep/M307W-background) -> {CSVOUT}")
    for r in rows:
        print(f"  {r['name']:<22} {r['rationale']}")


if __name__ == "__main__":
    main()
