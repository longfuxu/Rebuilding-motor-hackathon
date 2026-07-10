#!/usr/bin/env python3
"""Fetch a subunit MSA (a3m) from the free ColabFold MMseqs2 server (api.colabfold.com).
Same MSA source that `boltz --use_msa_server` uses. Saves seqs/gp17_subunit.a3m."""
import sys, time, os, io, tarfile, urllib.request, urllib.parse, urllib.error

OUT = os.path.dirname(os.path.abspath(__file__))
SEQ = "".join(l.strip() for l in open(f"{OUT}/seqs/gp17_subunit.fasta") if not l.startswith(">"))
HOST = "https://api.colabfold.com"

def post(path, data):
    req = urllib.request.Request(HOST+path, data=urllib.parse.urlencode(data).encode(),
            headers={"User-Agent":"gp17-designer"})
    return urllib.request.urlopen(req, timeout=120).read()

def get(path):
    req = urllib.request.Request(HOST+path, headers={"User-Agent":"gp17-designer"})
    return urllib.request.urlopen(req, timeout=120).read()

import json
print(f"submitting {len(SEQ)} aa to ColabFold MSA server (mode=env)...", flush=True)
r = json.loads(post("/ticket/msa", {"q": f">101\n{SEQ}\n", "mode": "env"}))
tid = r["id"]; print("ticket", tid, "status", r.get("status"), flush=True)

t0 = time.time()
while time.time()-t0 < 1200:
    st = json.loads(get(f"/ticket/{tid}"))
    s = st.get("status")
    if s in ("COMPLETE","ERROR"): print("status", s, flush=True); break
    time.sleep(8)
assert s == "COMPLETE", f"MSA failed: {st}"

blob = get(f"/result/download/{tid}")
tar = tarfile.open(fileobj=io.BytesIO(blob))
names = tar.getnames(); print("archive:", names, flush=True)
# uniref.a3m is the main alignment; may also have bfd/env
a3m_txt = None
for cand in ("uniref.a3m", "bfd.mgnify30.metaeuk30.smag30.a3m"):
    if cand in names:
        a3m_txt = tar.extractfile(cand).read().decode(); break
if a3m_txt is None:  # take the biggest .a3m
    a3ms = [n for n in names if n.endswith(".a3m")]
    a3m_txt = tar.extractfile(sorted(a3ms, key=lambda n: tar.getmember(n).size)[-1]).read().decode()

open(f"{OUT}/seqs/gp17_subunit.a3m","w").write(a3m_txt)
depth = a3m_txt.count("\n>")  # number of homolog records after query
print(f"saved seqs/gp17_subunit.a3m  depth~{a3m_txt.count(chr(62))} records", flush=True)
