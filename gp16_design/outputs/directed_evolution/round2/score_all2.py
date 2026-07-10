#!/usr/bin/env python3
"""Score every folded round-2 variant with score_variant.score_cif, caching each metric dict
to round2/scores/<name>.json. Resumable. The round-1 WT fold (copied into round2/folds) is the
M5-pocket geometry reference. We repoint score_variant.FOLDS at round2/folds so read_gate()
reads the round-2 result.json M1/M2 verdicts (not the round-1 folds dir)."""
import os, sys, json, glob, time
HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, PARENT)
import score_variant as sv

FOLDS = os.path.join(HERE, "folds")
SCORES = os.path.join(HERE, "scores")
# critical: read_gate() uses sv.FOLDS to find <name>.result.json -> point it at round2/folds
sv.FOLDS = FOLDS
sv.WT_CIF = os.path.join(FOLDS, "cp233_WT.cif")


def main(force=False):
    os.makedirs(SCORES, exist_ok=True)
    wt_cif = os.path.join(FOLDS, "cp233_WT.cif")
    assert os.path.exists(wt_cif), "WT fold missing in round2/folds"
    print("loading WT pocket reference ...", flush=True)
    wt_S = sv.setup(wt_cif)
    cifs = sorted(glob.glob(os.path.join(FOLDS, "*.cif")))
    for cif in cifs:
        name = os.path.basename(cif).replace(".cif", "")
        out = os.path.join(SCORES, f"{name}.json")
        if os.path.exists(out) and not force:
            print(f"[skip] {name}", flush=True)
            continue
        t0 = time.time()
        rec, _ = sv.score_cif(name, cif, wt_S=wt_S)
        json.dump(rec, open(out, "w"), indent=2)
        print(f"[score] {name}: PRS_NN={rec['prs_ATP_NN_coupling']} "
              f"conduit={rec['len_ATP_Y129_DNA']} grip={rec['DNA_inward_total']} "
              f"iface={rec['interSub_heavy_contacts']} breath20={rec['breathing_top20']} "
              f"pocket={rec['M5_pocket_rmsd_vs_wt']} m2={rec.get('m2_engaged')} "
              f"({round(time.time()-t0)}s)", flush=True)
    print("ALL SCORED", flush=True)


if __name__ == "__main__":
    main(force="--force" in sys.argv)
