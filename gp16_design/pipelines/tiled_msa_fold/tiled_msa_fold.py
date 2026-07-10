#!/usr/bin/env python3
"""tiled_msa_fold — fold a single-chain ring (N tandem copies of one domain + linkers)
with a BLOCK-DIAGONAL ("tiled") MSA via the free NVIDIA Boltz-2 NIM, then score M1/M2.

WHY THIS EXISTS
    Single-sequence structure prediction cannot evaluate a large single-chain ring:
    even the validated lead cp233 folds to the closed ring 0/5 under single-seq Boltz,
    but 5/5 with a tiled MSA. A tiled MSA gives every copy the MONOMER's evolutionary
    alignment on its own diagonal block (gaps off-block), so each subunit folds like the
    real domain while the covalent order enforces the ring. This is the only way to
    fold these constructs correctly on a predictor.

METHOD (mirrors the CS session's cofactor.fold: monomer MSA -> clean -> tile -> NIM)
    1. Monomer MSA: reuse a cached monomer a3m (e.g. handoff/core.a3m, 327-aa gp16 core)
       or fetch one from the free ColabFold MMseqs2 server for any new monomer.
    2. clean_a3m: drop lowercase insertion columns -> fixed-width match-state rows.
    3. derive_blocks: greedily map each construct copy back onto the monomer by exact
       substring match (handles circular permutations: a copy = M[cut:] + linker + M[:cut];
       and native-order tandems: a copy = M + linker). Auto-derives the score_m2 landmarks.
    4. build_tiled_a3m: query row = full construct; each monomer homolog row is placed
       (permuted) into ONE copy's residue columns, gaps elsewhere, header tagged _c{k}.
    5. fold via NIM /predict with the custom MSA (schema confirmed: polymers[i].msa =
       {"<db>": {"a3m": {"alignment": <text>, "format": "a3m"}}}).
    6. score with reproduce/score_m2.py (M2 = trans-R146 -> neighbour Walker-A < 8 A).

USAGE
    export NVIDIA_API_KEY=...   # from the main-checkout .env
    python tiled_msa_fold.py run manifests/cp233_WT.json --out ../../outputs/tiled_fold
    python tiled_msa_fold.py tile manifests/cp233_WT.json -o /tmp/cp233.a3m   # build a3m only
    python tiled_msa_fold.py msa <monomer_seq> -o mono.a3m                    # fetch monomer MSA

MANIFEST (json)
    { "name": "cp233_WT",
      "sequence": "<full single-chain construct>",
      "n_copies": 5,
      "monomer_a3m": "/abs/path/core.a3m",     # OR "monomer_ref": "<monomer seq>" to fetch
      "score": {"native_r146": 146, "walker": [24,31], "monomer_first_native_res": 4} }
    The `score` block is optional; gp16 defaults are used and the within-copy positions
    (r146_incopy, walker_incopy) + copy spans are DERIVED automatically from the tiling.
"""
import sys, os, json, time, hashlib, tarfile, io, argparse, subprocess
import urllib.request, urllib.error

NIM_URL = "https://health.api.nvidia.com/v1/biology/mit/boltz2/predict"
COLABFOLD = "https://api.colabfold.com"
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))          # gp16_design/
SCORE_M2 = os.path.join(REPO, "reproduce", "score_m2.py")
CACHE = os.path.join(HERE, "msa_cache")
AA = set("ACDEFGHIKLMNPQRSTVWY")


# ---------------------------------------------------------------- a3m utilities
def clean_a3m(a3m_text):
    """Parse a3m -> (query_seq, [(header, matchstr), ...]) with insertions dropped.

    Every returned matchstr has length == query match length (uppercase + '-' only),
    so the alignment is a fixed-width matrix ready to tile."""
    lines = [l.rstrip("\n") for l in a3m_text.splitlines() if l.strip()]
    recs, i = [], 0
    while i < len(lines):
        if lines[i].startswith(">"):
            hdr, seq = lines[i], lines[i + 1] if i + 1 < len(lines) else ""
            match = "".join(c for c in seq if c == "-" or c.isupper())
            recs.append((hdr, match)); i += 2
        else:
            i += 1
    qlen = len(recs[0][1])
    recs = [(h, m) for (h, m) in recs if len(m) == qlen]     # drop malformed rows
    return recs[0][1], recs


def fetch_colabfold_msa(seq, mode="env", timeout=1800):
    """Fetch a monomer a3m from the free ColabFold MMseqs2 API (disk-cached by seq hash)."""
    os.makedirs(CACHE, exist_ok=True)
    key = hashlib.sha1((mode + "|" + seq).encode()).hexdigest()[:16]
    cpath = os.path.join(CACHE, f"{key}.a3m")
    if os.path.exists(cpath):
        return open(cpath).read()
    def post(path, data):
        r = urllib.request.urlopen(urllib.request.Request(
            f"{COLABFOLD}{path}", data=urllib.parse.urlencode(data).encode()), timeout=120)
        return json.load(r)
    import urllib.parse
    sub = post("/ticket/msa", {"q": f">query\n{seq}", "mode": mode})
    tid = sub["id"]
    t0 = time.time()
    while time.time() - t0 < timeout:
        st = json.load(urllib.request.urlopen(f"{COLABFOLD}/ticket/{tid}", timeout=60))
        if st["status"] == "COMPLETE":
            break
        if st["status"] in ("ERROR", "UNKNOWN"):
            raise RuntimeError(f"colabfold MSA {st}")
        time.sleep(8)
    raw = urllib.request.urlopen(f"{COLABFOLD}/result/download/{tid}", timeout=300).read()
    tf = tarfile.open(fileobj=io.BytesIO(raw))
    a3m = ""
    for m in tf.getmembers():
        if m.name.endswith(".a3m"):
            a3m += tf.extractfile(m).read().decode(errors="replace")
    open(cpath, "w").write(a3m)
    return a3m


# --------------------------------------------------------- construct <-> monomer
def _extend(C, M, i, j, win=4, need=3):
    """Length of the run at (C[i:], M[j:]) tolerating isolated substitutions (point
    mutations, e.g. E119Q) but stopping at linkers (long non-matching stretches)."""
    L = 0
    while i + L < len(C) and j + L < len(M):
        if C[i + L] == M[j + L]:
            L += 1; continue
        m = w = 0                                            # peek past a mismatch
        while w < win and i + L + 1 + w < len(C) and j + L + 1 + w < len(M):
            m += C[i + L + 1 + w] == M[j + L + 1 + w]; w += 1
        if w and m >= need:
            L += 1                                           # isolated substitution -> bridge
        else:
            break
    return L


def derive_blocks(C, M, n_copies, min_run=12):
    """Map each of the n_copies onto the monomer M by greedy substring matching.

    Returns copies = [[(c_start, length, mono_start), ...], ...] (0-based). Each copy is
    a contiguous group of runs whose monomer columns tile M exactly once; unmatched
    stretches (linkers) are skipped. Works for circular permutations (2 runs/copy),
    native-order tandems (1 run/copy), and constructs carrying point mutations (the run
    is anchored by an exact >=min_run stretch, then extended across substitutions)."""
    runs, i, LM = [], 0, len(M)
    while i < len(C):
        bestL, bestj = 0, None
        for j in range(LM):                                  # anchor by longest EXACT match
            L = 0
            while i + L < len(C) and j + L < LM and C[i + L] == M[j + L]:
                L += 1
            if L > bestL:
                bestL, bestj = L, j
        if bestL >= min_run:
            full = _extend(C, M, i, bestj)                   # then span point mutations
            runs.append((i, full, bestj)); i += full
        else:
            i += 1
    copies, cur, covered = [], [], set()
    for r in runs:
        cur.append(r); covered |= set(range(r[2], r[2] + r[1]))
        if len(covered) == LM:
            copies.append(cur); cur, covered = [], set()
    if len(copies) != n_copies:
        raise RuntimeError(f"derive_blocks: found {len(copies)} full copies, expected {n_copies} "
                           f"(runs={len(runs)}); check monomer_ref / min_run")
    return copies


def score_args(C, copies, M, native_r146=146, walker=(24, 31), first_native=4):
    """Derive score_m2.py args from the tiling: --copies spans and the within-copy
    positions of R146 and Walker-A (in COPY 1, mirrored across copies by construction)."""
    def mono0(r):
        return r - first_native                              # native res -> 0-based monomer col
    starts = [min(r[0] for r in cp) for cp in copies]
    spans, r146_incopy, wlo, whi = [], None, None, None
    for ci, cp in enumerate(copies):
        m2c = {}
        for (cst, L, j) in cp:
            for t in range(L):
                m2c[j + t] = cst + t
        lo = starts[ci]
        hi = (starts[ci + 1] - 1) if ci + 1 < len(copies) else len(C) - 1
        spans.append((lo, hi))
        if ci == 0:
            r146_incopy = m2c[mono0(native_r146)] - lo + 1
            wpos = sorted(m2c[mono0(r)] for r in range(walker[0], walker[1] + 1))
            assert wpos == list(range(wpos[0], wpos[-1] + 1)), "Walker-A not contiguous in copy"
            wlo, whi = wpos[0] - lo + 1, wpos[-1] - lo + 1
    copystr = ",".join(f"A:{lo + 1}-{hi + 1}" for lo, hi in spans)
    return copystr, r146_incopy, f"{wlo}-{whi}"


def build_tiled_a3m(C, recs, copies, max_depth=None):
    """Block-diagonal tile: query row = construct C; each monomer homolog row (recs[1:])
    is placed into ONE copy's columns (permuted per the copy's runs), gaps elsewhere."""
    homologs = recs[1:]
    if max_depth:
        homologs = homologs[:max_depth]
    out = [">query", C]
    n = len(C)
    for k, cp in enumerate(copies):
        for hdr, match in homologs:
            row = ["-"] * n
            any_res = False
            for (cst, L, j) in cp:
                seg = match[j:j + L]
                for t, ch in enumerate(seg):
                    row[cst + t] = ch
                    if ch != "-":
                        any_res = True
            if any_res:
                out.append(f"{hdr.split()[0]}_c{k + 1}")
                out.append("".join(row))
    return "\n".join(out) + "\n"


# ------------------------------------------------------------------------- NIM
def fold_nim(seq, a3m_text, key, samples=1, recycling=3, retries=6):
    payload = json.dumps({
        "polymers": [{"id": "A", "molecule_type": "protein", "sequence": seq,
                      "msa": {"tiled": {"a3m": {"alignment": a3m_text, "format": "a3m"}}}}],
        "recycling_steps": recycling, "sampling_steps": 50, "diffusion_samples": samples,
    }).encode()
    for attempt in range(retries):
        try:
            req = urllib.request.Request(NIM_URL, data=payload,
                    headers={"Authorization": f"Bearer {key}", "Content-Type": "application/json"})
            t0 = time.time()
            resp = json.load(urllib.request.urlopen(req, timeout=1800))
            return resp, round(time.time() - t0, 1)
        except urllib.error.HTTPError as e:
            body = e.read()[:300].decode(errors="replace")
            if e.code in (429, 500, 502, 503, 504):
                wait = 20 * (attempt + 1)
                print(f"  HTTP {e.code}, backoff {wait}s ({body[:80]})", flush=True)
                time.sleep(wait); continue
            return {"__err": f"HTTP {e.code}: {body}"}, 0
        except Exception as e:
            print(f"  {repr(e)[:100]}, retry", flush=True); time.sleep(15)
    return {"__err": "exhausted retries"}, 0


def run_score(cif, copies_str, r146_incopy, walker_incopy):
    sc = subprocess.run(["python3", SCORE_M2, cif, "--copies", copies_str,
                         "--copy_start_res", "1", "--r146_incopy", str(r146_incopy),
                         "--walker_incopy", walker_incopy], capture_output=True, text=True)
    return sc.stdout, sc.stderr


def _engaged(out):
    import re
    m = re.search(r"M2:\s*(\d+)/(\d+)", out)
    return (int(m.group(1)), int(m.group(2))) if m else (0, 0)


def score_ring(cif, copies_str, r146_incopy, walker_incopy):
    """Handedness-robust M2. A near-C5 single-chain ring can wind either way in one
    diffusion sample, so R146 may engage the DESIGNED neighbour (k->k+1) or its mirror
    (k->k-1). Both are equally-closed rings; only the metric's assumed direction differs.
    Score BOTH directions and take the coherent maximum, reporting which handedness won."""
    fwd, _ = run_score(cif, copies_str, r146_incopy, walker_incopy)
    rev, _ = run_score(cif, ",".join(reversed(copies_str.split(","))), r146_incopy, walker_incopy)
    ef, n = _engaged(fwd); er, _ = _engaged(rev)
    if er > ef:
        eng, hand, best = er, "mirror(k->k-1)", rev
    else:
        eng, hand, best = ef, "designed(k->k+1)", fwd
    verdict = {"engaged": eng, "n": n, "handedness": hand,
               "m2_forward": ef, "m2_reverse": er,
               "score_forward": fwd, "score_reverse": rev}
    return verdict, best


# ------------------------------------------------------------------------- run
def load_manifest(path):
    m = json.load(open(path))
    seq = "".join(m["sequence"].split())
    assert set(seq) <= AA, f"non-AA chars in sequence: {set(seq)-AA}"
    m["sequence"] = seq
    return m


def get_monomer(m):
    if m.get("monomer_a3m"):
        return open(m["monomer_a3m"]).read()
    ref = "".join(m["monomer_ref"].split())
    return fetch_colabfold_msa(ref)


def do_tile(m, max_depth=None):
    C = m["sequence"]
    mono_a3m = get_monomer(m)
    Mseq, recs = clean_a3m(mono_a3m)
    copies = derive_blocks(C, Mseq, m["n_copies"])
    s = m.get("score", {})
    copies_str, r146, wlk = score_args(
        C, copies, Mseq,
        native_r146=s.get("native_r146", 146),
        walker=tuple(s.get("walker", [24, 31])),
        first_native=s.get("monomer_first_native_res", 4))
    a3m = build_tiled_a3m(C, recs, copies, max_depth=max_depth)
    meta = {"copies": copies_str, "r146_incopy": r146, "walker_incopy": wlk,
            "monomer_len": len(Mseq), "msa_depth": len(recs) - 1,
            "tiled_rows": a3m.count("\n>") + 1, "construct_len": len(C)}
    return a3m, meta


def do_run(manifest_path, outdir, max_depth=None, samples=1):
    key = os.environ["NVIDIA_API_KEY"]
    m = load_manifest(manifest_path)
    name = m["name"]
    os.makedirs(outdir, exist_ok=True)
    print(f"[{name}] construct {len(m['sequence'])} aa x {m['n_copies']} copies", flush=True)
    a3m, meta = do_tile(m, max_depth=max_depth)
    print(f"[{name}] tiled a3m: depth {meta['msa_depth']} -> {meta['tiled_rows']} rows; "
          f"copies={meta['copies']} r146_incopy={meta['r146_incopy']} walker={meta['walker_incopy']}",
          flush=True)
    open(os.path.join(outdir, f"{name}.tiled.a3m"), "w").write(a3m)
    resp, dt = fold_nim(m["sequence"], a3m, key, samples=samples)
    if "__err" in resp:
        res = {"name": name, "ok": False, "err": resp["__err"], **meta}
        json.dump(res, open(os.path.join(outdir, f"{name}.result.json"), "w"), indent=2)
        print(f"[{name}] FOLD FAILED: {resp['__err']}", flush=True)
        return res
    cif = os.path.join(outdir, f"{name}.cif")
    open(cif, "w").write(resp["structures"][0]["structure"])
    metrics = {k: resp.get(k) for k in ("confidence_scores", "ptm_scores", "iptm_scores",
                                        "complex_plddt_scores", "complex_pde_scores")}
    verdict, out = score_ring(cif, meta["copies"], meta["r146_incopy"], meta["walker_incopy"])
    m2 = (f"# M2(ring): {verdict['engaged']}/{verdict['n']} arginine fingers engaged "
          f"[{verdict['handedness']}]  (forward {verdict['m2_forward']}/{verdict['n']}, "
          f"reverse {verdict['m2_reverse']}/{verdict['n']})")
    m1 = next((l for l in out.splitlines() if l.startswith("# M1:")), "")
    m4 = next((l for l in out.splitlines() if l.startswith("# M4:")), "")
    res = {"name": name, "ok": True, "wall_s": dt, "metrics": metrics,
           "M2": m2, "M1": m1, "M4": m4, "ring": verdict, "score_full": out, **meta}
    json.dump(res, open(os.path.join(outdir, f"{name}.result.json"), "w"), indent=2)
    print(f"[{name}] {dt}s  {m2}", flush=True)
    print(f"[{name}]      {m1}", flush=True)
    return res


def do_rescore(outdir):
    """Re-score existing <name>.cif in outdir with the handedness-robust ring M2, using
    each result.json's stored tiling landmarks. No re-folding."""
    for rj in sorted(glob_result(outdir)):
        d = json.load(open(rj))
        cif = os.path.join(outdir, f"{d['name']}.cif")
        if not (d.get("ok") and os.path.exists(cif)):
            continue
        verdict, out = score_ring(cif, d["copies"], d["r146_incopy"], d["walker_incopy"])
        d["ring"] = verdict
        d["M2"] = (f"# M2(ring): {verdict['engaged']}/{verdict['n']} arginine fingers engaged "
                   f"[{verdict['handedness']}]  (forward {verdict['m2_forward']}/{verdict['n']}, "
                   f"reverse {verdict['m2_reverse']}/{verdict['n']})")
        d["M1"] = next((l for l in out.splitlines() if l.startswith("# M1:")), d.get("M1", ""))
        d["M4"] = next((l for l in out.splitlines() if l.startswith("# M4:")), d.get("M4", ""))
        d["score_full"] = out
        json.dump(d, open(rj, "w"), indent=2)
        print(f"[{d['name']}] {d['M2']}")


def glob_result(outdir):
    import glob as _g
    return _g.glob(os.path.join(outdir, "*.result.json"))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run"); r.add_argument("manifest")
    r.add_argument("--out", default=os.path.join(REPO, "outputs", "tiled_fold"))
    r.add_argument("--max_depth", type=int, default=None)
    r.add_argument("--samples", type=int, default=1)
    t = sub.add_parser("tile"); t.add_argument("manifest"); t.add_argument("-o", required=True)
    t.add_argument("--max_depth", type=int, default=None)
    mm = sub.add_parser("msa"); mm.add_argument("seq"); mm.add_argument("-o", required=True)
    rs = sub.add_parser("rescore")
    rs.add_argument("--out", default=os.path.join(REPO, "outputs", "tiled_fold"))
    a = ap.parse_args()
    if a.cmd == "run":
        do_run(a.manifest, a.out, max_depth=a.max_depth, samples=a.samples)
    elif a.cmd == "rescore":
        do_rescore(a.out)
    elif a.cmd == "tile":
        m = load_manifest(a.manifest)
        a3m, meta = do_tile(m, max_depth=a.max_depth)
        open(a.o, "w").write(a3m)
        print(json.dumps(meta, indent=2)); print("wrote", a.o)
    elif a.cmd == "msa":
        open(a.o, "w").write(fetch_colabfold_msa("".join(a.seq.split())))
        print("wrote", a.o)


if __name__ == "__main__":
    main()
