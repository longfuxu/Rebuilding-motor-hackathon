#!/usr/bin/env python3
"""Greedy-epistasis ranking. Background = M307W. The decisive question:
    does ANY M307W+X beat the M307W BACKGROUND on force OR coordination, by more than the
    fold-noise floor, with all gates (M1,M2,M5,soft,biochem) intact?

Two views:
  (1) round-1 gated pipeline vs WT (../rank_variants.py, repointed here) -> ranked_variants.csv,
      ranking_summary.json. de_M307W sits in this panel as the labelled BACKGROUND row, so its
      own gains vs WT (coord +6.4%, force +3.5%) are the bar every M307W+X must clear.
  (2) head-to-head vs the M307W BACKGROUND -> vs_m307w.csv (the decisive table): per-variant
      fractional gain of coord (PRS_NN, higher=better) and force (conduit length, lower=better)
      RELATIVE TO M307W, compared to the WT-replicate noise floor. STACKS = gate-pass AND
      (Δcoord_vs_M307W > coord_noise OR Δforce_vs_M307W > force_noise).
"""
import os, sys, json, csv
HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, PARENT)
import rank_variants as rv

SCORES = os.path.join(HERE, "scores")
rv.SCORES = SCORES
rv.OUT_CSV = os.path.join(SCORES, "ranked_variants.csv")

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

BG = "de_M307W"   # background


def frac_vs(v, bg, higher):
    if bg in (None, 0) or v is None:
        return 0.0
    g = (v - bg) / abs(bg)
    return g if higher else -g


def main():
    rv.main()   # ranked_variants.csv + ranking_summary.json (vs WT)
    summ = json.load(open(os.path.join(SCORES, "ranking_summary.json")))
    S = _safe_load_scores()
    noise = summ["noise"]                       # {coord, force, grip, soft} fractional (|WT-WTrep|)
    bg = S.get(BG)
    assert bg is not None, "M307W background score missing"

    coord_key, force_key = rv.PROXY["coord"][0], rv.PROXY["force"][0]
    coord_noise, force_noise = noise["coord"], noise["force"]

    rows = []
    pv = summ["per_variant"]
    for n in summ["ranked"]:
        if n == BG:
            continue
        r = S[n]
        d_coord = frac_vs(r.get(coord_key), bg.get(coord_key), True)    # higher PRS = better
        d_force = frac_vs(r.get(force_key), bg.get(force_key), False)   # lower conduit = better
        gate = pv[n]["gate"]; passed = all(gate.values())
        stack_coord = passed and d_coord > coord_noise
        stack_force = passed and d_force > force_noise
        stacks = stack_coord or stack_force
        rows.append(dict(
            name=n, muts=pv[n].get("muts", ""), pass_all=passed,
            d_coord_vs_M307W=round(d_coord, 4), d_force_vs_M307W=round(d_force, 4),
            coord_noise=round(coord_noise, 4), force_noise=round(force_noise, 4),
            stacks_coord=stack_coord, stacks_force=stack_force, STACKS=stacks,
            soft_retention=round(pv[n]["soft_ret"], 4),
            gates_failed=";".join(k for k, v in gate.items() if not v) or "none",
            # absolute metrics + gains vs WT for context
            gain_coord_vs_WT=round(pv[n]["gains"]["coord"], 4),
            gain_force_vs_WT=round(pv[n]["gains"]["force"], 4),
            prs_ATP_NN_coupling=r.get(coord_key), len_ATP_Y129_DNA=r.get(force_key)))

    # sort: stackers first, then by best relative margin over M307W (coord or force / its noise)
    def margin(row):
        return max(row["d_coord_vs_M307W"] / coord_noise, row["d_force_vs_M307W"] / force_noise)
    rows.sort(key=lambda x: (x["STACKS"], x["pass_all"], margin(x)), reverse=True)

    cols = ["name", "muts", "pass_all", "STACKS", "stacks_coord", "stacks_force",
            "d_coord_vs_M307W", "d_force_vs_M307W", "coord_noise", "force_noise",
            "soft_retention", "gain_coord_vs_WT", "gain_force_vs_WT",
            "prs_ATP_NN_coupling", "len_ATP_Y129_DNA", "gates_failed"]
    with open(os.path.join(SCORES, "vs_m307w.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)

    stackers = [r for r in rows if r["STACKS"]]
    print("\n" + "=" * 78)
    print(f"GREEDY EPISTASIS: does any M307W+X beat the M307W background above noise?")
    print(f"  M307W background: PRS_NN={bg.get(coord_key):.6f}  conduit={bg.get(force_key):.4f}")
    print(f"  noise floor: coord={coord_noise*100:.2f}%  force={force_noise*100:.2f}%")
    print(f"  gate-passing M307W+X: {sum(1 for r in rows if r['pass_all'])}/{len(rows)}")
    print(f"\n{'variant':<18} {'pass':<5} {'dCoord%':>8} {'dForce%':>8} {'softRet':>7} {'STACKS?':>8}")
    for r in rows:
        print(f"{r['name']:<18} {str(r['pass_all']):<5} {r['d_coord_vs_M307W']*100:>+7.1f}% "
              f"{r['d_force_vs_M307W']*100:>+7.1f}% {r['soft_retention']:>7.3f} "
              f"{'YES' if r['STACKS'] else ('coord' if r['stacks_coord'] else ('force' if r['stacks_force'] else 'no')):>8}")
    print("\n" + "=" * 78)
    if stackers:
        b = stackers[0]
        print(f"STACKING EXISTS: {len(stackers)} M307W+X beat M307W above noise. "
              f"Best accumulated variant = {b['name']} "
              f"(dCoord {b['d_coord_vs_M307W']*100:+.1f}%, dForce {b['d_force_vs_M307W']*100:+.1f}%).")
    else:
        print("NO STACKING: no gate-passing M307W+X beats M307W above the noise floor on either "
              "coord or force. The 'singles do not stack, even on the best background' conclusion "
              "is now robust across a broad second-site scan.")
    json.dump(dict(noise=noise, background=BG,
                   bg_metrics={coord_key: bg.get(coord_key), force_key: bg.get(force_key)},
                   stackers=[r["name"] for r in stackers], rows=rows),
              open(os.path.join(SCORES, "greedy_summary.json"), "w"), indent=2)
    print(f"\nwrote {os.path.join(SCORES,'vs_m307w.csv')}\nwrote {os.path.join(SCORES,'greedy_summary.json')}")


if __name__ == "__main__":
    main()
