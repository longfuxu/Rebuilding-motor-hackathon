#!/usr/bin/env python3
"""MM-GBSA interface-energy validation (amber ff14SB + OBC2, single-trajectory) on WT + the
round-2 top-3 (clean gate-passing) variants. Orthogonal, non-DL, non-ENM read on the M2
coupler interface: more negative dG_bind = tighter inter-subunit coupling.

Reuses ../mmgbsa_top.py's runner and ../../independent_validation/interface_energy/
mmgbsa_interface.py, repointed at round2/. Targets default to ranking_summary.json top3 but
can be overridden on the command line (space-separated variant names)."""
import os, sys, json, csv
HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, PARENT)
import mmgbsa_top as mm

FOLDS = os.path.join(HERE, "folds")
SCORES = os.path.join(HERE, "scores")
MMDIR = os.path.join(SCORES, "mmgbsa")
# repoint the round-1 runner at round-2 dirs
mm.FOLDS = FOLDS
mm.SCORES = SCORES
mm.MMDIR = MMDIR


def main():
    os.makedirs(MMDIR, exist_ok=True)
    if len(sys.argv) > 1:
        targets = ["cp233_WT"] + [a for a in sys.argv[1:] if a != "cp233_WT"]
    else:
        summ = json.load(open(os.path.join(SCORES, "ranking_summary.json")))
        targets = ["cp233_WT"] + summ["top3"]
    print("MM-GBSA targets:", targets, flush=True)
    out = {}
    for name in targets:
        csvp = mm.run_one(name)
        out[name] = mm.mean_dg(csvp)
    wt = out.get("cp233_WT")
    rows = []
    print("\n==== MM-GBSA mean interface dG_bind (kcal/mol; more negative = tighter) ====")
    for name in targets:
        dg = out[name]
        dd = (dg - wt) if (dg is not None and wt is not None) else None
        rows.append(dict(name=name, mean_dG_bind=round(dg, 2) if dg is not None else None,
                         delta_vs_WT=round(dd, 2) if dd is not None else None))
        print(f"  {name:<20} mean dG_bind = {dg}  (dd vs WT = {dd})")
    json.dump(dict(mean_dG_bind=out, rows=rows),
              open(os.path.join(MMDIR, "mmgbsa_summary.json"), "w"), indent=2)
    with open(os.path.join(MMDIR, "mmgbsa_summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "mean_dG_bind", "delta_vs_WT"])
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {os.path.join(MMDIR,'mmgbsa_summary.csv')}")


if __name__ == "__main__":
    main()
