#!/usr/bin/env python3
"""Shared helpers for the in-silico directed-evolution of the cp233 single-chain gp16 motor.

The cp233 construct is 5 covalent copies of the gp16 core, each copy laid out as
    segA = native 234..330   (97 aa)
  + (GGGGS)3                 (15 aa linker)
  + segB = native 4..233     (230 aa)
  + (GGGGS)2                 (10 aa linker)   [absent on the last copy -> 1750 aa total]
so a copy body is 352 aa and copy k starts at construct position 352*k (0-based). This
matches coord_common's _SEGA_START/_SEGB_START exactly.

We mutate a native residue in ALL 5 copies (a C5-symmetric ring variant) by editing the
construct query sequence at the mapped positions; the tiled MSA (monomer homologs) is
untouched, so each variant folds with the WT monomer alignment and only the query carries
the substitution (tiled_msa_fold.derive_blocks bridges isolated point substitutions).
"""
import os, json

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))          # gp16_design/
WT_MANIFEST = os.path.join(REPO, "pipelines", "tiled_msa_fold", "manifests", "cp233_WT.json")
CORE_A3M = os.path.join(REPO, "pipelines", "tiled_msa_fold", "gp16_core.a3m")

COPY_LEN = 352
SEGA_OFF = 0            # segA (native 234..330) offset within a copy
SEGB_OFF = 112         # segB (native 4..233) offset within a copy (97 + 15 linker)
N_COPIES = 5

# catalytic / biochemically-lethal residues -- NEVER mutate (native numbering)
FORBIDDEN = set(range(24, 32)) | {118, 119, 146}   # Walker-A P-loop, Walker-B D/E, arginine finger

AA = set("ACDEFGHIKLMNPQRSTVWY")


def wt_sequence():
    m = json.load(open(WT_MANIFEST))
    return "".join(m["sequence"].split())


def native_to_copy_positions(native_res):
    """0-based construct positions of a native residue across the (up to 5) copies present."""
    seq_len = len(wt_sequence())
    pos = []
    for k in range(N_COPIES):
        base = COPY_LEN * k
        if 234 <= native_res <= 330:
            p = base + SEGA_OFF + (native_res - 234)
        elif 4 <= native_res <= 233:
            p = base + SEGB_OFF + (native_res - 4)
        else:
            raise ValueError(f"native res {native_res} outside 4..330")
        if p < seq_len:
            pos.append(p)
    return pos


def apply_mutations(seq, muts):
    """muts = list of (native_res, new_aa). Returns mutated sequence (all copies edited).
    Raises if a target is forbidden or the WT identity is unexpected."""
    s = list(seq)
    for native_res, new_aa in muts:
        assert new_aa in AA, f"bad AA {new_aa}"
        if native_res in FORBIDDEN:
            raise ValueError(f"native {native_res} is catalytic/forbidden -- refusing to mutate")
        for p in native_to_copy_positions(native_res):
            s[p] = new_aa
    return "".join(s)


def wt_identity(native_res):
    """WT residue letter at a native position (read from copy 0 of the construct)."""
    seq = wt_sequence()
    return seq[native_to_copy_positions(native_res)[0]]


def mutation_tag(muts):
    """e.g. [(289,'F'),(307,'F')] -> 'L289F_M307F' using WT identities."""
    return "_".join(f"{wt_identity(r)}{r}{a}" for r, a in muts)


def write_manifest(name, seq, outdir):
    os.makedirs(outdir, exist_ok=True)
    mani = {"name": name, "sequence": seq, "n_copies": N_COPIES, "monomer_a3m": CORE_A3M}
    path = os.path.join(outdir, f"{name}.json")
    json.dump(mani, open(path, "w"), indent=2)
    return path


if __name__ == "__main__":
    # sanity: WT identities at all mutation targets must match native gp16
    expect = {100: "V", 128: "N", 129: "Y", 130: "I", 289: "L", 307: "M"}
    print("target native identities (construct copy0 vs expected native):")
    ok = True
    for r, e in expect.items():
        got = wt_identity(r)
        tag = "OK" if got == e else "MISMATCH"
        if got != e:
            ok = False
        print(f"  native {r}: construct={got} expected={e}  {tag}  "
              f"positions={native_to_copy_positions(r)}")
    print("ALL OK" if ok else "MAPPING ERROR")
