#!/usr/bin/env python3
"""ESM-2 naturalness (sequence-level, orthogonal to structure predictors): per-residue
log-likelihood under a protein language model. Compares the de-novo RFdiffusion connector
region against the native gp16 motif in the same chain. Low LL = 'unnatural' sequence."""
import torch, esm, numpy as np

# small-ish model, CPU-friendly single forward pass
MODELS = ["esm2_t30_150M_UR50D", "esm2_t12_35M_UR50D"]
model = None
for m in MODELS:
    try:
        model, alphabet = getattr(esm.pretrained, m)(); print("loaded", m, flush=True); break
    except Exception as e:
        print("skip", m, repr(e)[:60], flush=True)
model.eval(); bc = alphabet.get_batch_converter()


def per_res_ll(seq):
    _, _, toks = bc([("x", seq)])
    with torch.no_grad():
        lp = torch.log_softmax(model(toks)["logits"][0], dim=-1)
    ll = []
    for i, aa in enumerate(seq):
        ti = alphabet.get_idx(aa)
        ll.append(float(lp[i + 1, ti]))  # +1 for BOS
    return np.array(ll)


def main():
    GP16 = open("/Users/longfu/.claude/jobs/48b7a543/tmp/gp16_4_330.seq").read().strip()
    print(f"NATIVE gp16 (327aa) mean LL: {per_res_ll(GP16).mean():.3f}", flush=True)
    recs = []; cur = None
    for l in open("/Users/longfu/.claude/jobs/48b7a543/tmp/mpnn_work/all_designs.fasta"):
        if l.startswith(">"): cur = l[1:].strip(); recs.append([cur, ""])
        elif recs: recs[-1][1] += l.strip()
    for name, s in [r for r in recs if r[0].startswith("sal_L50")]:
        ll = per_res_ll(s); N = len(s)
        motif = np.concatenate([ll[:327], ll[N - 327:]])
        conn = ll[327:N - 327]
        print(f"{name}: motif meanLL {motif.mean():.3f} | de-novo connector meanLL {conn.mean():.3f} "
              f"(Δ {conn.mean()-motif.mean():+.3f})", flush=True)


if __name__ == "__main__":
    main()
