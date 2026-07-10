#!/usr/bin/env python3
"""Aggregate the directed-evolution variant scores into gates + a weighted power score, rank,
and write the ranked CSV. The scoring function (see DELIVERABLES_AND_NORTHSTAR.md):

    score = GATE(M1 & M2 & M5-pocket & soft-mode-floor)
            x [ w1*z(coordination) + w2*z(force) + w3*z(grip)
                + w4*z(interface) + w5*z(soft-mode-retention) ]

  GATES (hard, reject if fail) -- an over-rigidified motor that fails any of these is a
  dead brick, not a stronger motor:
    G1 M1  : closed compact sequential ring (fold)
    G2 M2  : >=4/5 handedness-robust arginine fingers engaged (fold)
    G3 M5  : ATP-pocket CA geometry preserved (pocket RMSD vs WT fold < POCKET_MAX)
    G4 soft: ring open/close (breathing) overlap with the softest modes retained
             (>= SOFT_MIN x WT) -- the power-stroke competence floor.

  POWER PROXIES (maximize, expressed as fractional gain vs WT, then z-scored across the
  gate-passing panel so weights are comparable). We use ONLY the proxies that survive the
  WT-replicate noise test (|WT - WTrep|):
    coordination = PRS ATP-site -> ring-neighbour ATP-site coupling   (noise ~6%)
    force        = ATP->Y129->DNA conduit path length, LOWER = stiffer (noise ~1.5%)
    grip         = # DNA-contact sidechains lining the pore inward     (noise ~6%)
    soft-mode    = breathing overlap retention (also the G4 gate)      (noise ~1%)
  This is exactly the northstar's 4-term formula. Y129 node betweenness (noise ~26%) and the
  raw inter-subunit heavy-atom contact count (noise ~25%) are too fold-noisy to rank on, so
  they are REPORTED ONLY (flagged); interface energy is validated on the winners by MM-GBSA.

  A WT replicate fold (cp233_WTrep) gives the no-mutation NOISE FLOOR: a proxy gain smaller
  than |WT - WTrep| is within fold noise and is flagged (not trusted).

  Caveat: static proxies RANK candidates; absolute force/power is not computable statically
  (ground truth = single-molecule, follow-on program). Never judged on pLDDT.
"""
import os, sys, json, glob, csv
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SCORES = os.path.join(HERE, "scores")
OUT_CSV = os.path.join(SCORES, "ranked_variants.csv")

# ---- gate thresholds ----
POCKET_MAX = 2.0        # A, M5 pocket CA RMSD vs WT fold
SOFT_MIN = 0.85         # breathing-overlap retention floor (fraction of WT)

# ---- biochem gate (phi29 bulk mutation data; gp16_design/outputs/excel_mutation_analysis) ----
# A variant that mutates any residue biochem-flagged SEVERE/LETHAL is rejected regardless of
# how well it folds/scores. Y129 is near-lethal (Y129A = severe packaging defect) -- its
# indispensability CONFIRMS it is the force hub our proxy already flags (high Y129 betweenness),
# so we rigidify AROUND it, never mutate it.
BIOCHEM_SEVERE = {129, 58, 105, 328, 330, 28, 32, 106, 126, 158, 234, 294}
BIOCHEM_TOLERANT = {6, 53, 56, 99, 127, 128, 222, 233}   # biochem-blessed safe residues

def muts_residues(muts_tag):
    import re
    return [int(x) for x in re.findall(r"\d+", muts_tag or "")]

# ---- weights (transparent, a-priori; power-proxy weights not yet biochem-calibrated) ----
WEIGHTS = dict(coord=0.30, force=0.30, grip=0.20, soft=0.20)

# proxy -> (metric key, higher_is_better). Only low-noise proxies enter the score.
PROXY = dict(
    coord=("prs_ATP_NN_coupling", True),
    force=("len_ATP_Y129_DNA", False),      # lower path length = stiffer/stronger transmission
    grip=("DNA_inward_total", True),
    soft=("breathing_top20", True),
)
# reported but NOT scored (too fold-noisy per WT replicate)
REPORT_ONLY = dict(
    betweenness=("Y129_betweenness", True),
    iface=("interSub_heavy_contacts", True),
)


def load_scores():
    d = {}
    for f in glob.glob(os.path.join(SCORES, "*.json")):
        r = json.load(open(f))
        d[r["name"]] = r
    return d


def load_library():
    lib = {}
    p = os.path.join(SCORES, "variants.csv")
    if os.path.exists(p):
        for row in csv.DictReader(open(p)):
            lib[row["name"]] = row
    return lib


def frac_gain(v, wt, higher):
    if wt in (None, 0) or v is None:
        return 0.0
    g = (v - wt) / abs(wt)
    return g if higher else -g


def main():
    S = load_scores()
    lib = load_library()
    assert "cp233_WT" in S, "WT score missing -- run score_all.py"
    wt = S["cp233_WT"]
    wtrep = S.get("cp233_WTrep")

    # --- noise floor (|WT - WTrep| fractional) per proxy ---
    noise = {}
    for pk, (mk, hi) in PROXY.items():
        if wtrep and wt.get(mk) not in (None, 0):
            noise[pk] = abs(frac_gain(wtrep[mk], wt[mk], hi))
        else:
            noise[pk] = None

    variants = [n for n in S if n not in ("cp233_WT", "cp233_WTrep")]

    # --- gains + gates ---
    rec = {}
    for n in variants:
        r = S[n]
        gains = {pk: frac_gain(r.get(mk), wt.get(mk), hi) for pk, (mk, hi) in PROXY.items()}
        soft_ret = (r.get("breathing_top20") or 0) / (wt.get("breathing_top20") or 1)
        pocket = r.get("M5_pocket_rmsd_vs_wt")
        mres = muts_residues(lib.get(n, {}).get("muts", n.replace("de_", "")))
        severe_hit = sorted(set(mres) & BIOCHEM_SEVERE)
        gate = dict(
            G1_M1=bool(r.get("m1_pass")),
            G2_M2=bool(r.get("m2_pass")),
            G3_M5=bool(pocket is not None and pocket < POCKET_MAX),
            G4_soft=bool(soft_ret >= SOFT_MIN),
            G5_biochem=bool(len(severe_hit) == 0),
        )
        rec[n] = dict(r=r, gains=gains, soft_ret=soft_ret, pocket=pocket,
                      gate=gate, pass_all=all(gate.values()),
                      severe_hit=severe_hit,
                      tolerant_hit=sorted(set(mres) & BIOCHEM_TOLERANT))

    # --- z-score each proxy gain across GATE-PASSING variants ---
    passers = [n for n in variants if rec[n]["pass_all"]]
    z = {pk: {} for pk in PROXY}
    for pk in PROXY:
        vals = np.array([rec[n]["gains"][pk] for n in passers]) if passers else np.array([])
        mu, sd = (vals.mean(), vals.std()) if len(vals) else (0.0, 1.0)
        sd = sd if sd > 1e-9 else 1.0
        for n in variants:
            z[pk][n] = (rec[n]["gains"][pk] - mu) / sd

    for n in variants:
        score = sum(WEIGHTS[pk] * z[pk][n] for pk in PROXY)
        rec[n]["z"] = {pk: z[pk][n] for pk in PROXY}
        rec[n]["power_score"] = score if rec[n]["pass_all"] else float("-inf")

    ranked = sorted(variants, key=lambda n: rec[n]["power_score"], reverse=True)

    # --- write CSV ---
    cols = ["rank", "name", "muts", "pass_all", "G1_M1", "G2_M2", "G3_M5", "G4_soft",
            "G5_biochem", "biochem_note", "power_score", "coord_gain", "force_gain", "grip_gain",
            "soft_retention", "prs_ATP_NN_coupling", "len_ATP_Y129_DNA",
            "DNA_inward_total", "breathing_top20", "min_pore_radius", "M5_pocket_rmsd_vs_wt",
            "m2_engaged", "Y129_betweenness_noisy", "interSub_heavy_contacts_noisy",
            "flag", "beats_wt_coord", "beats_wt_force",
            "above_noise_coord", "above_noise_force", "rationale"]
    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        # WT reference row first
        w.writerow(dict(rank=0, name="cp233_WT", muts="(reference)", pass_all="ref",
                        power_score=0.0,
                        prs_ATP_NN_coupling=wt.get("prs_ATP_NN_coupling"),
                        len_ATP_Y129_DNA=wt.get("len_ATP_Y129_DNA"),
                        DNA_inward_total=wt.get("DNA_inward_total"),
                        breathing_top20=wt.get("breathing_top20"),
                        min_pore_radius=wt.get("min_pore_radius"),
                        Y129_betweenness_noisy=wt.get("Y129_betweenness"),
                        interSub_heavy_contacts_noisy=wt.get("interSub_heavy_contacts"),
                        M5_pocket_rmsd_vs_wt=0.0, soft_retention=1.0))
        for i, n in enumerate(ranked, 1):
            R = rec[n]; r = R["r"]; g = R["gains"]
            w.writerow(dict(
                rank=i, name=n, muts=lib.get(n, {}).get("muts", n.replace("de_", "")),
                pass_all=R["pass_all"], **{k: R["gate"][k] for k in R["gate"]},
                biochem_note=("SEVERE:" + ",".join(map(str, R["severe_hit"]))) if R["severe_hit"]
                else ("tolerant:" + ",".join(map(str, R["tolerant_hit"])) if R["tolerant_hit"]
                      else "untested"),
                power_score=round(R["power_score"], 4) if np.isfinite(R["power_score"]) else "GATED",
                coord_gain=round(g["coord"], 4), force_gain=round(g["force"], 4),
                grip_gain=round(g["grip"], 4),
                soft_retention=round(R["soft_ret"], 4),
                prs_ATP_NN_coupling=r.get("prs_ATP_NN_coupling"),
                len_ATP_Y129_DNA=r.get("len_ATP_Y129_DNA"),
                DNA_inward_total=r.get("DNA_inward_total"),
                breathing_top20=r.get("breathing_top20"),
                min_pore_radius=r.get("min_pore_radius"),
                Y129_betweenness_noisy=r.get("Y129_betweenness"),
                interSub_heavy_contacts_noisy=r.get("interSub_heavy_contacts"),
                M5_pocket_rmsd_vs_wt=R["pocket"], m2_engaged=r.get("m2_engaged"),
                flag=lib.get(n, {}).get("flag", ""),
                beats_wt_coord=g["coord"] > 0, beats_wt_force=g["force"] > 0,
                above_noise_coord=(noise["coord"] is None or abs(g["coord"]) > noise["coord"]),
                above_noise_force=(noise["force"] is None or abs(g["force"]) > noise["force"]),
                rationale=lib.get(n, {}).get("rationale", "")))

    # --- console summary + return structured result for the report ---
    print(f"\nNOISE FLOOR (|WT - WTrep| fractional): "
          + ", ".join(f"{k}={('%.4f'%v) if v is not None else 'NA'}" for k, v in noise.items()))
    print(f"WT: PRS_NN={wt['prs_ATP_NN_coupling']} conduitLen={wt['len_ATP_Y129_DNA']} "
          f"grip={wt['DNA_inward_total']} breath20={wt['breathing_top20']}")
    print(f"\n{'rank':<4} {'variant':<16} {'pass':<5} {'score':>7} {'coord%':>8} "
          f"{'force%':>8} {'grip%':>7} {'softRet':>7} gates")
    for i, n in enumerate(ranked, 1):
        R = rec[n]; g = R["gains"]
        gates_str = "PASS" if R["pass_all"] else "FAIL:" + ",".join(
            k for k, v in R["gate"].items() if not v)
        sc = f"{R['power_score']:.3f}" if np.isfinite(R["power_score"]) else "GATED"
        print(f"{i:<4} {n:<16} {str(R['pass_all']):<5} {sc:>7} {g['coord']*100:>7.1f}% "
              f"{g['force']*100:>7.1f}% {g['grip']*100:>6.1f}% "
              f"{R['soft_ret']:>7.3f} {gates_str}")

    top3 = [n for n in ranked if rec[n]["pass_all"]][:3]
    print(f"\nTOP-3 gate-passing: {top3}")
    json.dump(dict(ranked=ranked, top3=top3, noise=noise,
                   wt={k: wt.get(PROXY[k][0]) for k in PROXY},
                   passers=passers, weights=WEIGHTS,
                   thresholds=dict(POCKET_MAX=POCKET_MAX, SOFT_MIN=SOFT_MIN),
                   per_variant={n: dict(power_score=rec[n]["power_score"] if np.isfinite(rec[n]["power_score"]) else None,
                                        gains=rec[n]["gains"], soft_ret=rec[n]["soft_ret"],
                                        pocket=rec[n]["pocket"], gate=rec[n]["gate"],
                                        severe_hit=rec[n]["severe_hit"], tolerant_hit=rec[n]["tolerant_hit"],
                                        muts=lib.get(n, {}).get("muts", ""),
                                        flag=lib.get(n, {}).get("flag", ""),
                                        metrics={PROXY[k][0]: rec[n]["r"].get(PROXY[k][0]) for k in PROXY})
                                for n in variants}),
              open(os.path.join(SCORES, "ranking_summary.json"), "w"), indent=2)
    print(f"\nwrote {OUT_CSV}\nwrote {os.path.join(SCORES,'ranking_summary.json')}")


if __name__ == "__main__":
    main()
