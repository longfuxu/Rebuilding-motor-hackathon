#!/usr/bin/env python3
"""Fold gp17 constructs via the free NVIDIA Boltz-2 NIM, with a TILED (block-diagonal) MSA.

Tiling: the single ATPase-domain subunit MSA (ColabFold) is placed on the diagonal, one
block per copy, gaps everywhere else (incl. linker columns). This gives each copy its own
coevolution signal without cross-copy pairing (the project's fix for the single-seq
confound; see DESIGN_METHODS_AND_CREDIBILITY.md §4). The hosted NIM accepts a precomputed
a3m via polymers[].msa.<db>.a3m.alignment.
"""
import sys, os, json, time, urllib.request, urllib.error

OUT = os.path.dirname(os.path.abspath(__file__))
KEY = [l.split('=',1)[1].strip().strip('"').strip("'")
       for l in open('/Users/longfu/Developer/claude-science-hackthon/.env') if l.startswith('NVIDIA_API_KEY')][0]
URL = "https://health.api.nvidia.com/v1/biology/mit/boltz2/predict"
MAX_DEPTH = 256   # homologs per block (payload/size cap)


ALLOWED = set("ACDEFGHIKLMNPQRSTVWY-")   # 20 aa + gap; map X/U/*/null etc. to '-'

def _matchstate(r):
    """Keep match columns (drop lowercase insertions); map any non-standard char to '-'."""
    return "".join((c if c in ALLOWED else "-") for c in r if not c.islower() and c != "\x00")

def load_a3m(path):
    """Return (query_seq, [homolog_matchstate_strings]). Sanitised for the NIM a3m parser."""
    recs = []
    hdr, seq = None, []
    for L in open(path):
        L = L.rstrip("\n")
        if L.startswith(">"):
            if hdr is not None:
                recs.append("".join(seq))
            hdr, seq = L, []
        elif L:
            seq.append(L)
    if hdr is not None:
        recs.append("".join(seq))
    query = "".join(c for c in recs[0] if c not in "-\x00" and not c.islower())
    homologs = [_matchstate(r) for r in recs[1:]]
    return query, homologs

def clean_a3m(query, homologs):
    """Rebuild a plain (untiled) a3m string from sanitised match-state rows."""
    lines = [">query", query]
    for i, h in enumerate(homologs):
        if len(h) == len(query):
            lines += [f">h{i}", h]
    return "\n".join(lines) + "\n"


def tiled_a3m(construct_seq, subunit_query, homologs, copy_starts, copy_len, depth=MAX_DEPTH):
    """Block-diagonal a3m string for a single-chain construct."""
    assert len(subunit_query) == copy_len, (len(subunit_query), copy_len)
    total = len(construct_seq)
    homs = [h for h in homologs if len(h) == copy_len][:depth]
    lines = [">query", construct_seq]
    for k, s in enumerate(copy_starts):     # s is 1-based
        pre, post = s - 1, total - (s - 1) - copy_len
        for i, h in enumerate(homs):
            lines.append(f">c{k}_{i}")
            lines.append("-" * pre + h + "-" * post)
    return "\n".join(lines) + "\n", len(homs)


def fold(name, seq, a3m_text=None, label="", timeout=1800, polymers=None):
    if polymers is None:
        poly = {"id": "A", "molecule_type": "protein", "sequence": seq}
        if a3m_text is not None:
            poly["msa"] = {"colabfold_env": {"a3m": {"alignment": a3m_text, "format": "a3m"}}}
        polymers = [poly]
    payload = json.dumps({"polymers": polymers}).encode()
    print(f"[{name}] {len(seq)} aa | msa={'tiled' if a3m_text else 'single-seq'} "
          f"| payload {len(payload)/1e6:.1f} MB ...", flush=True)
    t0 = time.time()
    req = urllib.request.Request(URL, data=payload,
            headers={"Authorization": f"Bearer {KEY}", "Content-Type": "application/json"})
    try:
        resp = json.load(urllib.request.urlopen(req, timeout=timeout))
    except urllib.error.HTTPError as e:
        print(f"[{name}] HTTP {e.code}: {e.read()[:500].decode(errors='replace')}", flush=True)
        return None
    dt = time.time() - t0
    cif = resp["structures"][0]["structure"]
    metrics = {k: resp.get(k) for k in ("confidence_scores", "ptm_scores", "iptm_scores",
                                        "complex_plddt_scores", "complex_pde_scores")}
    metrics.update(name=name, len_aa=len(seq), wall_s=round(dt, 1),
                   predictor="boltz2-nim", msa=label or ("tiled" if a3m_text else "single-seq"))
    open(f"{OUT}/fold_results/{name}.cif", "w").write(cif)
    json.dump(metrics, open(f"{OUT}/fold_results/{name}_metrics.json", "w"), indent=2)
    print(f"[{name}] OK {dt:.0f}s | conf={metrics['confidence_scores']} "
          f"ptm={metrics['ptm_scores']} plddt={metrics['complex_plddt_scores']} -> fold_results/{name}.cif", flush=True)
    return metrics


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    constructs = json.load(open(f"{OUT}/seqs/constructs.json"))
    subq, homs = load_a3m(f"{OUT}/seqs/gp17_subunit.a3m")
    print(f"subunit MSA: query {len(subq)} aa, {len(homs)} homologs "
          f"({sum(1 for h in homs if len(h)==len(subq))} full-length)", flush=True)

    if which in ("all", "monomer"):
        # positive control: monomer with its own (untiled, sanitised) MSA
        mono_a3m = clean_a3m(subq, homs)
        fold("gp17_atpase_monomer", constructs["gp17_atpase_monomer"]["seq"],
             a3m_text=mono_a3m, label="colabfold-msa (monomer control)")

    if which in ("all", "df"):
        c = constructs["gp17_directfusion_L20"]
        a3m, d = tiled_a3m(c["seq"], subq, homs, c["copy_starts"], c["copy_len"])
        print(f"tiled a3m: {d} homologs/block x5 blocks", flush=True)
        open(f"{OUT}/fold_results/gp17_directfusion_L20_tiled.a3m", "w").write(a3m)
        fold("gp17_directfusion_L20", c["seq"], a3m_text=a3m, label="tiled block-diagonal (colabfold)")

    if which in ("df_ss",):
        # negative control: same construct, SINGLE-SEQ (no MSA) -> expected confounded (framework 4)
        c = constructs["gp17_directfusion_L20"]
        fold("gp17_directfusion_L20_singleseq", c["seq"], a3m_text=None, label="single-seq (control)")

    if which in ("all", "native"):
        # native apo homo-pentamer reference: 5 separate chains, each with the subunit MSA.
        sub = constructs["gp17_atpase_monomer"]["seq"]
        mono_a3m = clean_a3m(subq, homs)
        polymers = [{"id": cid, "molecule_type": "protein", "sequence": sub,
                     "msa": {"colabfold_env": {"a3m": {"alignment": mono_a3m, "format": "a3m"}}}}
                    for cid in "ABCDE"]
        fold("gp17_native_pentamer", sub, label="native 5-chain ref (per-chain MSA)", polymers=polymers)
