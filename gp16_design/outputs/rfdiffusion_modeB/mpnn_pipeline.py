#!/usr/bin/env python3
"""Thread native gp16 seq onto RFdiffusion connector backbones, then ProteinMPNN
design ONLY the connector (catalytic motif frozen). Emits full designed sequences.

Layout: out res 1..327 = gp16 4-330 (subunit A), connector = 328..(N-327),
subunit B = (N-326)..N = gp16 4-330. Connector length = N-654."""
import os, sys, glob, subprocess, json

MPNN = "/Users/longfu/.claude/jobs/48b7a543/tmp/ProteinMPNN"
PY = "/Users/longfu/miniforge3/bin/python3.13"
BB = "/Users/longfu/Developer/claude-science-hackthon/.claude/worktrees/export-gp16-results/gp16_design/outputs/rfdiffusion_modeB/backbones"
WORK = "/Users/longfu/.claude/jobs/48b7a543/tmp/mpnn_work"
os.makedirs(WORK, exist_ok=True)
GP16 = open("/Users/longfu/.claude/jobs/48b7a543/tmp/gp16_4_330.seq").read().strip()  # 327 aa, gp16 res4-330
MOTIF = 327
AA1TO3 = {"A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","Q":"GLN","E":"GLU","G":"GLY","H":"HIS",
          "I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE","P":"PRO","S":"SER","T":"THR","W":"TRP",
          "Y":"TYR","V":"VAL"}


def thread(pdb_in, pdb_out):
    """Relabel motif residue names to native gp16; connector stays as-is (GLY)."""
    # find N
    ns = [int(l[22:26]) for l in open(pdb_in) if l[:4] == "ATOM" and l[12:16].strip() == "CA"]
    N = max(ns); L = N - 2 * MOTIF
    def native_name(rn):
        if 1 <= rn <= MOTIF:              # subunit A -> gp16 pos rn
            return AA1TO3[GP16[rn - 1]]
        if rn >= N - MOTIF + 1:           # subunit B -> gp16 pos (rn-(N-MOTIF))
            return AA1TO3[GP16[rn - (N - MOTIF) - 1]]
        return None                       # connector: leave unchanged
    out = []
    for l in open(pdb_in):
        if l[:4] == "ATOM":
            rn = int(l[22:26]); nm = native_name(rn)
            if nm:
                l = l[:17] + nm + l[20:]
        out.append(l)
    open(pdb_out, "w").write("".join(out))
    return N, L


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("ERR:", " ".join(cmd[:3]), r.stderr[-400:])
    return r


def main():
    results = {}
    for pdb in sorted(glob.glob(f"{BB}/*.pdb")):
        tag = os.path.splitext(os.path.basename(pdb))[0]
        d = f"{WORK}/{tag}"; os.makedirs(d, exist_ok=True)
        threaded = f"{d}/threaded.pdb"
        N, L = thread(pdb, threaded)
        conn = list(range(MOTIF + 1, N - MOTIF + 1))      # connector positions (1-indexed in chain A)
        # parse
        run([PY, f"{MPNN}/helper_scripts/parse_multiple_chains.py",
             "--input_path", d, "--output_path", f"{d}/parsed.jsonl"])
        # fixed positions: specify the connector as NON-fixed (design), everything else fixed
        pos = " ".join(str(p) for p in conn)
        run([PY, f"{MPNN}/helper_scripts/make_fixed_positions_dict.py",
             "--input_path", f"{d}/parsed.jsonl", "--output_path", f"{d}/fixed.jsonl",
             "--chain_list", "A", "--position_list", pos, "--specify_non_fixed"])
        # design
        run([PY, f"{MPNN}/protein_mpnn_run.py",
             "--jsonl_path", f"{d}/parsed.jsonl", "--fixed_positions_jsonl", f"{d}/fixed.jsonl",
             "--out_folder", d, "--num_seq_per_target", "4", "--sampling_temp", "0.1",
             "--seed", "37", "--batch_size", "1", "--omit_AAs", "C",
             "--use_soluble_model", "--model_name", "v_48_020"])
        # collect designed seqs (skip the first record = input seq)
        fa = f"{d}/seqs/threaded.fa"
        seqs = []
        if os.path.isfile(fa):
            cur = None
            for line in open(fa):
                if line.startswith(">"):
                    cur = line.strip(); seqs.append(["", cur])
                elif seqs:
                    seqs[-1][0] += line.strip()
            # first entry is the input (native+Gly connector); rest are designs
            designs = seqs[1:]
            # verify catalytic residues preserved in designs (R146 -> pos146, K30, E119)
            ok = []
            for s, hdr in designs:
                cat_ok = len(s) == N and s[145] == "R" and s[29] == "K" and s[118] == "E"
                ok.append(cat_ok)
            results[tag] = {"N": N, "L": L, "n_designs": len(designs),
                            "cat_preserved": sum(ok), "connector_positions": [conn[0], conn[-1]]}
            # write full designed sequences to a combined fasta
            with open(f"{WORK}/all_designs.fasta", "a") as out:
                for i, (s, hdr) in enumerate(designs):
                    out.write(f">{tag}_d{i}\n{s}\n")
        else:
            results[tag] = {"error": "no mpnn output", "N": N, "L": L}
        print(tag, results[tag])
    json.dump(results, open(f"{WORK}/mpnn_summary.json", "w"), indent=2)
    print("SUMMARY", json.dumps(results))


if __name__ == "__main__":
    open("/Users/longfu/.claude/jobs/48b7a543/tmp/mpnn_work/all_designs.fasta", "w").close()
    main()
