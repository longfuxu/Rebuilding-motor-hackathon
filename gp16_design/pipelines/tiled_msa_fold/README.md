# tiled_msa_fold

Fold a **single-chain ring** (N tandem copies of one domain joined by linkers) with a
**block-diagonal ("tiled") MSA** on the free **NVIDIA Boltz-2 NIM**, then score the ring
geometry (M1) and the trans arginine-finger interface (M2).

## Why this exists (the load-bearing reason)

Single-sequence structure prediction **cannot** evaluate a large single-chain ring. Even
the validated lead **cp233** folds to the closed ring **0/5** under single-seq Boltz, but
**5/5** with a tiled MSA. A tiled MSA gives every copy the *monomer's* evolutionary
alignment on its own diagonal block (gaps off-block), so each subunit folds like the real
domain while the covalent order enforces the ring. **Any big single-chain ring must be
folded this way** — do not judge these constructs on single-seq folds or on global pTM.

## Method (mirrors the CS session `cofactor.fold`: monomer MSA -> clean -> tile -> NIM)

1. **Monomer MSA** — reuse a cached monomer a3m (`gp16_core.a3m`, the 327-aa gp16 core =
   native residues 4-330, 1167 homologs) or fetch a fresh one for any new monomer from the
   free ColabFold MMseqs2 server (`msa` subcommand; disk-cached under `msa_cache/`).
2. **`clean_a3m`** — drop lowercase a3m insertion columns -> a fixed-width match-state
   matrix (every row == monomer length).
3. **`derive_blocks`** — greedily map each construct copy back onto the monomer by exact
   substring match. Handles **circular permutations** (a copy = `M[cut:] + linker + M[:cut]`,
   2 runs/copy) and **native-order tandems** (a copy = `M + linker`, 1 run/copy). From the
   mapping it auto-derives the `score_m2.py` landmarks (copy spans, `r146_incopy`,
   `walker_incopy`) — no hand-tuning per construct.
4. **`build_tiled_a3m`** — query row = full construct; each monomer homolog row is placed
   (permuted) into exactly one copy's residue columns, gaps everywhere else, header tagged
   `_c{k}`. Result for a 5-copy gp16 ring = 5836 rows (1 query + 1167×5), matching the CS
   session's reference `pentamer.a3m`.
5. **`fold_nim`** — POST to `https://health.api.nvidia.com/v1/biology/mit/boltz2/predict`
   with the custom MSA.
6. **score** — `../../reproduce/score_m2.py`. **M2** = for each subunit, min heavy-atom
   distance from its R146 guanidinium to the *sequential-neighbour* copy's Walker-A
   (res 24-31); engaged < 8 Å; pass = ≥4/5. **M1** = ring-closure geometry
   (compact + planar + sequential register). Never uses global pTM.

## Usage

```bash
export NVIDIA_API_KEY=...          # from the main-checkout .env (health.api.nvidia.com)

# End-to-end: tile -> fold -> score, writes <name>.cif / .tiled.a3m / .result.json
python tiled_msa_fold.py run manifests/cp233_WT.json --out ../../outputs/tiled_fold

# Build the tiled a3m only (inspect / feed elsewhere); prints derived score args
python tiled_msa_fold.py tile manifests/cp233_WT.json -o /tmp/cp233.a3m

# Fetch a monomer MSA for a brand-new domain from ColabFold
python tiled_msa_fold.py msa SLFYNPQK...FRKMR -o new_core.a3m
```

Optional flags: `--max_depth N` (cap monomer MSA rows before tiling), `--samples N`
(diffusion samples). Full depth (1167) posts a ~10 MB payload and is accepted by the NIM.

### Manifest format (`manifests/*.json`)

```json
{ "name": "cp233_WT",
  "sequence": "<full single-chain construct>",
  "n_copies": 5,
  "monomer_a3m": "/abs/path/gp16_core.a3m",
  "score": {"native_r146": 146, "walker": [24, 31], "monomer_first_native_res": 4} }
```

`monomer_a3m` (cached) or `monomer_ref` (a monomer sequence to fetch). The `score` block is
optional — the gp16 defaults shown are used, and the within-copy landmark positions +
copy spans are **derived** from the tiling, so the same manifest shape works for any CP site.

## Confirmed Boltz-2 NIM custom-MSA schema (empirically validated; the key gotcha)

The NIM `/predict` request accepts a per-polymer `msa` field. Shape (validated by 422
probes — the field IS parsed and consumed, not ignored):

```json
{ "polymers": [ { "id": "A", "molecule_type": "protein", "sequence": "<query>",
    "msa": { "<db_name>": { "a3m": { "alignment": "<a3m text>", "format": "a3m" } } } } ] }
```

- `format` must be one of `{a3m, csv, fasta, sto}` (case-sensitive).
- `<db_name>` is an arbitrary key (this pipeline uses `"tiled"`).
- The a3m's first row must equal the polymer `sequence` (the query).
- Response: `structures[0].structure` (mmCIF text), plus `complex_plddt_scores`,
  `confidence_scores`, `ptm_scores`, `pae`, etc. HTTP 429 = rate limit (back off, retry).

## Files

- `tiled_msa_fold.py` — library + CLI (`run` / `tile` / `msa`).
- `gp16_core.a3m` — cached gp16 monomer MSA (327 aa, 1167 homologs); mirror of
  `cs_sessions/session1_main_gp16/handoff/core.a3m`.
- `manifests/` — per-construct manifests (cp233 WT + E119Q ×2, cp285, cp297, RFdiff rings).
- `msa_cache/` — ColabFold MSAs fetched for new monomers.
- Outputs -> `../../outputs/tiled_fold/` (`<name>.cif`, `<name>.tiled.a3m`,
  `<name>.result.json`) + `SUMMARY.md` / `summary.csv`.
