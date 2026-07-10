#!/usr/bin/env python3
"""Rank the round-2 variants with the round-1 gated scoring function (../rank_variants.py),
repointed at round2/scores, then add an explicit head-to-head vs the round-1 winner M307W.

Outputs (round2/scores/):
  ranked_variants.csv, ranking_summary.json   (from rank_variants, round-2 panel)
  vs_m307w.csv, vs_m307w.json                 (per-proxy gains of each round-2 variant
                                               minus M307W's gains, both vs cp233_WT)
"""
import os, sys, json, csv
HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, PARENT)
import rank_variants as rv

SCORES = os.path.join(HERE, "scores")
R1_SCORES = os.path.join(PARENT, "scores")   # round-1 scores (for the M307W benchmark)

# repoint the round-1 ranker at the round-2 score dir
rv.SCORES = SCORES
rv.OUT_CSV = os.path.join(SCORES, "ranked_variants.csv")

# round2/scores also holds non-score json (esm_plausibility.json = a list, vs_m307w.json,
# ranking_summary.json). Filter load_scores to per-variant score records (dict with "name").
import glob as _glob, json as _json
def _safe_load_scores():
    d = {}
    for f in _glob.glob(os.path.join(SCORES, "*.json")):
        try:
            r = _json.load(open(f))
        except Exception:
            continue
        if isinstance(r, dict) and "name" in r and "prs_ATP_NN_coupling" in r:
            d[r["name"]] = r
    return d
rv.load_scores = _safe_load_scores


def m307w_benchmark():
    """M307W (round-1 winner) proxy gains vs cp233_WT, using the same PROXY defs + WT ref
    as the round-2 panel. Returns (m307w_gains, wt_metrics) or (None, None) if absent."""
    wt_p = os.path.join(SCORES, "cp233_WT.json")
    m_p = os.path.join(R1_SCORES, "de_M307W.json")
    if not (os.path.exists(wt_p) and os.path.exists(m_p)):
        return None, None
    wt = json.load(open(wt_p)); m = json.load(open(m_p))
    gains = {pk: rv.frac_gain(m.get(mk), wt.get(mk), hi) for pk, (mk, hi) in rv.PROXY.items()}
    metrics = {rv.PROXY[pk][0]: m.get(rv.PROXY[pk][0]) for pk in rv.PROXY}
    return dict(gains=gains, metrics=metrics,
                soft_ret=(m.get("breathing_top20") or 0) / (wt.get("breathing_top20") or 1),
                pocket=m.get("M5_pocket_rmsd_vs_wt")), wt


def main():
    rv.main()   # writes round2/scores/ranked_variants.csv + ranking_summary.json
    summ = json.load(open(os.path.join(SCORES, "ranking_summary.json")))
    m307w, _wt = m307w_benchmark()

    print("\n" + "=" * 72)
    if m307w is None:
        print("M307W round-1 score not found -- skipping head-to-head.")
        return
    mg = m307w["gains"]
    print("HEAD-TO-HEAD vs round-1 winner M307W (proxy gains vs WT; delta = round2 - M307W)")
    print(f"  M307W gains vs WT:  coord={mg['coord']*100:+.1f}%  force={mg['force']*100:+.1f}%  "
          f"grip={mg['grip']*100:+.1f}%  soft_ret={m307w['soft_ret']:.3f}")
    print(f"{'variant':<18} {'pass':<5} {'d_coord':>8} {'d_force':>8} {'d_grip':>8} "
          f"{'beats_M307W?':>26}")
    rows = []
    pv = summ["per_variant"]
    order = summ["ranked"]
    for n in order:
        R = pv[n]
        g = R["gains"]
        d = {pk: g[pk] - mg[pk] for pk in rv.PROXY}
        passed = all(R["gate"].values())
        # "beats M307W" = gate-passing AND >= M307W on the two low-noise scored proxies
        # (force is the round-1 robust axis; coord is the axis we tried to move this round)
        beats_force = g["force"] >= mg["force"]
        beats_coord = g["coord"] >= mg["coord"]
        verdict = []
        if passed and beats_force and beats_coord:
            verdict.append("BOTH coord+force")
        elif passed and beats_force:
            verdict.append("force only")
        elif passed and beats_coord:
            verdict.append("coord only")
        elif not passed:
            verdict.append("GATED")
        else:
            verdict.append("no")
        rows.append(dict(name=n, muts=R.get("muts", ""), pass_all=passed,
                         d_coord=round(d["coord"], 4), d_force=round(d["force"], 4),
                         d_grip=round(d["grip"], 4), d_soft=round(g["soft"], 4),
                         beats_force=beats_force, beats_coord=beats_coord,
                         verdict=";".join(verdict),
                         v_coord_gain=round(g["coord"], 4), v_force_gain=round(g["force"], 4),
                         m307w_coord_gain=round(mg["coord"], 4),
                         m307w_force_gain=round(mg["force"], 4)))
        print(f"{n:<18} {str(passed):<5} {d['coord']*100:>+7.1f}% {d['force']*100:>+7.1f}% "
              f"{d['grip']*100:>+7.1f}% {';'.join(verdict):>26}")

    with open(os.path.join(SCORES, "vs_m307w.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    json.dump(dict(m307w_gains=mg, m307w_soft_ret=m307w["soft_ret"],
                   m307w_pocket=m307w["pocket"], rows=rows),
              open(os.path.join(SCORES, "vs_m307w.json"), "w"), indent=2)
    print(f"\nwrote {os.path.join(SCORES, 'vs_m307w.csv')}")


if __name__ == "__main__":
    main()
