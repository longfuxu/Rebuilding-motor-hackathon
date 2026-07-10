#!/usr/bin/env python3
"""Independent interface-energy validation (MM-GBSA, amber ff14SB + OBC2, single-trajectory)
on WT + the top-ranked gate-passing variants. This is the orthogonal, non-deep-learning,
non-ENM read on the M2 coupler interface -- more negative dG_bind = tighter inter-subunit
coupling. It replaces the fold-noisy heavy-atom contact-count surrogate for the winners.

Reuses gp16_design/outputs/independent_validation/interface_energy/mmgbsa_interface.py.
Copies/landmarks are read from each fold's result.json (constant for the cp233 layout).
"""
import os, sys, json, csv, glob, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
FOLDS = os.path.join(HERE, "folds")
SCORES = os.path.join(HERE, "scores")
MMDIR = os.path.join(SCORES, "mmgbsa")
MMGBSA = os.path.abspath(os.path.join(HERE, "..", "independent_validation",
                                     "interface_energy", "mmgbsa_interface.py"))
PY = sys.executable


def run_one(name):
    cif = os.path.join(FOLDS, f"{name}.cif")
    rj = json.load(open(os.path.join(FOLDS, f"{name}.result.json")))
    copies = rj["copies"]; r146 = rj["r146_incopy"]; walk = rj["walker_incopy"]
    outdir = os.path.join(MMDIR, name)
    os.makedirs(outdir, exist_ok=True)
    csvp = os.path.join(outdir, f"{name}_interface_energy.csv")
    if os.path.exists(csvp):
        print(f"[skip] {name} MM-GBSA exists", flush=True)
        return csvp
    print(f"[mmgbsa] {name} ...", flush=True)
    cmd = [PY, MMGBSA, "--design", cif, "--copies", copies,
           "--r146_incopy", str(r146), "--walker_incopy", walk,
           "--label", name, "--out", outdir, "--platform", "CPU"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"[mmgbsa] {name} FAILED:\n{r.stderr[-800:]}", flush=True)
        return None
    return csvp


def mean_dg(csvp):
    if not csvp or not os.path.exists(csvp):
        return None
    vals = [float(row["dG_bind"]) for row in csv.DictReader(open(csvp))]
    return sum(vals) / len(vals) if vals else None


def main():
    summ = json.load(open(os.path.join(SCORES, "ranking_summary.json")))
    targets = ["cp233_WT"] + summ["top3"]
    print("MM-GBSA targets:", targets)
    out = {}
    for name in targets:
        csvp = run_one(name)
        out[name] = mean_dg(csvp)
    wt = out.get("cp233_WT")
    rows = []
    print("\n==== MM-GBSA mean interface dG_bind (kcal/mol; more negative = tighter) ====")
    for name in targets:
        dg = out[name]
        dd = (dg - wt) if (dg is not None and wt is not None) else None
        rows.append(dict(name=name, mean_dG_bind=round(dg, 2) if dg else None,
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
