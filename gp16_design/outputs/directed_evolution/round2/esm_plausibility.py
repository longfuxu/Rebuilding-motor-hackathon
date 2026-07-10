#!/usr/bin/env python3
"""ESM-2 (t30 150M, cached locally) masked-marginal plausibility filter for the round-2
variants. For each mutated native position we mask it on the gp16 MONOMER native sequence
(native 4..330, the biologically meaningful context -- NOT the 1750 aa fusion) and score
    LLR = log P(mutant AA | context) - log P(WT AA | context)
This is the standard ESM variant-effect proxy. A strongly negative LLR flags an implausible
substitution; ~0 or positive means ESM finds it as/more natural than WT. Sum over both sites
for doubles. Writes round2/scores/esm_plausibility.{json,csv}. Best-effort: if the model can't
load, the pipeline still ranks on the structural gates/proxies without it.
"""
import os, sys, json, csv, re
HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, PARENT)
import de_common as de

SCORES = os.path.join(HERE, "scores")


def native_monomer():
    """native gp16 residues 4..330 in native order, from construct copy 0."""
    wt = de.wt_sequence()
    return "".join(wt[de.native_to_copy_positions(r)[0]] for r in range(4, 331)), 4


def variant_muts():
    """{variant_name: [(native_res, wt_aa, mut_aa), ...]} from round2/scores/variants.csv."""
    out = {}
    p = os.path.join(SCORES, "variants.csv")
    for row in csv.DictReader(open(p)):
        tag = row["muts"]
        if not tag:
            continue
        muts = []
        for m in re.finditer(r"([A-Z])(\d+)([A-Z])", tag):
            wt_aa, res, mut_aa = m.group(1), int(m.group(2)), m.group(3)
            muts.append((res, wt_aa, mut_aa))
        out[row["name"]] = muts
    return out


def main():
    import torch, esm
    seq, first = native_monomer()
    model, alphabet = esm.pretrained.esm2_t30_150M_UR50D()
    model.eval()
    bc = alphabet.get_batch_converter()
    tok = alphabet.tok_to_idx

    def masked_logprobs(seq, mask_pos0):
        """log-softmax over vocab at a single masked sequence position (0-based in seq)."""
        _, _, toks = bc([("wt", seq)])
        toks = toks.clone()
        toks[0, mask_pos0 + 1] = alphabet.mask_idx        # +1 for BOS
        with torch.no_grad():
            logits = model(toks)["logits"][0, mask_pos0 + 1]
        return torch.log_softmax(logits, dim=-1)

    muts_by_var = variant_muts()
    # verify WT identities align with the monomer string
    rows = []
    for name, muts in muts_by_var.items():
        total_llr = 0.0
        per = []
        for res, wt_aa, mut_aa in muts:
            pos0 = res - first
            assert seq[pos0] == wt_aa, f"{name}: monomer[{res}]={seq[pos0]} != WT {wt_aa}"
            lp = masked_logprobs(seq, pos0)
            llr = float(lp[tok[mut_aa]] - lp[tok[wt_aa]])
            total_llr += llr
            per.append(dict(res=res, wt=wt_aa, mut=mut_aa, llr=round(llr, 3),
                            wt_logp=round(float(lp[tok[wt_aa]]), 3),
                            mut_logp=round(float(lp[tok[mut_aa]]), 3)))
        rows.append(dict(name=name, muts="_".join(f"{w}{r}{m}" for r, w, m in muts),
                         esm_llr_sum=round(total_llr, 3),
                         esm_plausible=bool(total_llr > -6.0),   # lenient structural filter
                         per_site=per))
        print(f"{name:<18} LLR_sum={total_llr:+.2f}  {'plausible' if total_llr>-6.0 else 'IMPLAUSIBLE'}"
              f"   {['%s%d%s:%.2f'%(p['wt'],p['res'],p['mut'],p['llr']) for p in per]}")

    json.dump(rows, open(os.path.join(SCORES, "esm_plausibility.json"), "w"), indent=2)
    with open(os.path.join(SCORES, "esm_plausibility.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "muts", "esm_llr_sum", "esm_plausible"])
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in ["name", "muts", "esm_llr_sum", "esm_plausible"]})
    print(f"\nwrote {os.path.join(SCORES,'esm_plausibility.csv')}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"\n[esm] SKIPPED (best-effort): {repr(e)[:200]}")
        sys.exit(0)
