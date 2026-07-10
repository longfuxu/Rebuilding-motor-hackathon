#!/usr/bin/env python3
"""Score every folded directed-evolution variant with score_variant.score_cif and cache the
per-variant metric dict to scores/<name>.json. Resumable (skips existing). The WT fold is
loaded once and reused as the M5-pocket geometry reference for all variants."""
import os, sys, json, glob, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score_variant as sv

HERE = os.path.dirname(os.path.abspath(__file__))
FOLDS = os.path.join(HERE, "folds")
SCORES = os.path.join(HERE, "scores")


def main(force=False):
    os.makedirs(SCORES, exist_ok=True)
    wt_cif = os.path.join(FOLDS, "cp233_WT.cif")
    assert os.path.exists(wt_cif), "WT fold missing"
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
        # WT itself: pocket ref = self -> 0 by construction; use wt_S so pocket=0 exactly
        rec, _ = sv.score_cif(name, cif, wt_S=wt_S)
        json.dump(rec, open(out, "w"), indent=2)
        print(f"[score] {name}: PRS_NN={rec['prs_ATP_NN_coupling']} "
              f"Ybtw={rec['Y129_betweenness']} grip={rec['DNA_inward_total']} "
              f"iface={rec['interSub_heavy_contacts']} breath20={rec['breathing_top20']} "
              f"pocket={rec['M5_pocket_rmsd_vs_wt']} ({round(time.time()-t0)}s)", flush=True)
    print("ALL SCORED", flush=True)


if __name__ == "__main__":
    main(force="--force" in sys.argv)
