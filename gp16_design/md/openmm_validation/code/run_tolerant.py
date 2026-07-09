#!/usr/bin/env python3
"""Disconnect-tolerant runner for the remaining gp16 systems.
Detached VM jobs survive transient CLI websocket drops, so this NEVER relaunches a
running job on a single 'not found' — it retries. Only after the session is confirmed
gone (empty `colab sessions`) does it recreate + relaunch. Downloads each system on DONE.

Usage: python run_tolerant.py B C     (systems to run, in order)
Env: GP16_NS (3), GP16_EQ (100)."""
import sys, os, time, tarfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import colab_helpers as C

HERE = os.path.dirname(os.path.abspath(__file__))
INP = os.path.abspath(os.path.join(HERE, "..", "inputs"))
RES = os.path.abspath(os.path.join(HERE, "..", "results"))
os.makedirs(RES, exist_ok=True)
INPUTS = {"A": "A_apo.pdb", "B": "B_7jqq_helical.pdb", "C": "C_design.cif"}
NS = float(os.environ.get("GP16_NS", "3"))
EQ = float(os.environ.get("GP16_EQ", "100"))
SYS = [s for s in sys.argv[1:] if s in INPUTS] or ["B", "C"]

def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

def _dead(o, e):
    b = ((o or "") + (e or "")).lower()
    return ("not found" in b) or ("no active session" in b) or ("no such session" in b)

def session_gone():
    rc, o, e = C.sessions(timeout=40)
    return "gp16md" not in (o or "")

def upload_all():
    for f in ("md_driver.py", "analyze.py"):
        C.upload(os.path.join(HERE, f), f"/content/{f}")
    for v in INPUTS.values():
        C.upload(os.path.join(INP, v), f"/content/{v}")

def ensure_session():
    if not session_gone():
        try:
            rc, o, e = C.exec_code("import openmm,mdtraj,pdbfixer;print('DEPS_OK')", timeout=60)
            if "DEPS_OK" in (o or ""): log("session alive+deps"); return
        except Exception: pass
    log("recreating fresh A100 session ...")
    C.new_session("A100", timeout=300)
    C.exec_code("import subprocess;subprocess.Popen(\"pip install -q openmm pdbfixer mdtraj "
                ">/content/pip.txt 2>&1 && python -c 'import openmm,pdbfixer,mdtraj' && "
                "touch /content/PIP_DONE || touch /content/PIP_FAIL\",shell=True);print('pip')", timeout=45)
    t0 = time.time()
    while time.time() - t0 < 25*60:
        rc, o, e = C.exec_code("import os;print('OK' if os.path.exists('/content/PIP_DONE') "
                               "else('FAIL' if os.path.exists('/content/PIP_FAIL') else 'w'))", timeout=40)
        if "OK" in (o or ""): break
        if "FAIL" in (o or ""): raise RuntimeError("pip failed")
        time.sleep(20)
    upload_all(); log("session ready")

STATUS = ("import os\n"
          "print('DONE' if os.path.exists('/content/mdout/{s}.DONE') else "
          "('FAIL' if os.path.exists('/content/mdout/{s}.FAIL') else 'run'))\n"
          "p='/content/mdout/{s}.log'\n"
          "print(open(p).read()[-140:].replace(chr(10),' | ') if os.path.exists(p) else 'no log')\n")

def vm_state(sysid):
    """Return ('DONE'|'FAIL'|'run'|'absent', tail) or raise on transient (caller retries)."""
    rc, o, e = C.exec_code(
        "import os,glob\n"
        f"d='/content/mdout/{sysid}'\n"
        f"print('DONE' if os.path.exists('/content/mdout/{sysid}.DONE') else "
        f"('FAIL' if os.path.exists('/content/mdout/{sysid}.FAIL') else "
        f"('run' if (os.path.isdir(d) and os.path.exists(d+'/.lock')) else 'absent')))\n"
        f"p='/content/mdout/{sysid}.log'\n"
        f"print(open(p).read()[-140:].replace(chr(10),' | ') if os.path.exists(p) else '')",
        timeout=45, retries=1)
    if _dead(o, e): raise ConnectionError("transient")
    lines = (o or "").strip().splitlines()
    st = lines[0].strip() if lines else "absent"
    tail = lines[1] if len(lines) > 1 else ""
    return st, tail

def launch(sysid):
    inp = INPUTS[sysid]
    cmd = (f"cd /content && python md_driver.py --system {sysid} --input {inp} --mode implicit "
           f"--ns {NS} --equil_ps {EQ} --report_ps 20 --gb_cutoff_nm 2.0 --min_iters 2000 "
           f"--timestep_fs 4 --hmr 1 --out /content/mdout/{sysid} > /content/mdout/{sysid}.log 2>&1 && "
           f"python analyze.py --system {sysid} --top /content/mdout/{sysid}/{sysid}_start.pdb "
           f"--traj /content/mdout/{sysid}/{sysid}_prod.dcd --out /content/mdout/{sysid} "
           f"--report_ps 20 >> /content/mdout/{sysid}.log 2>&1 && touch /content/mdout/{sysid}.DONE "
           f"|| touch /content/mdout/{sysid}.FAIL")
    C.exec_code(f"import os,subprocess;os.makedirs('/content/mdout',exist_ok=True)\n"
                f"for p in ('/content/mdout/{sysid}.DONE','/content/mdout/{sysid}.FAIL'):\n"
                f"    (os.path.exists(p) and os.remove(p))\n"
                f"subprocess.Popen({cmd!r},shell=True);print('launched {sysid}')", timeout=45)
    log(f"{sysid}: launched (ns={NS})")

def fetch(sysid):
    C.exec_code(f"import subprocess;subprocess.run('cd /content/mdout && tar czf /content/{sysid}.tgz "
                f"{sysid}/{sysid}_timeseries.csv {sysid}/{sysid}_contacts.csv {sysid}/{sysid}_rmsf.csv "
                f"{sysid}/{sysid}_summary.json {sysid}/{sysid}_start.pdb {sysid}/{sysid}_final.pdb "
                f"{sysid}/{sysid}_prod.csv {sysid}/{sysid}_equil.csv {sysid}.log',shell=True);"
                f"print('tarred')", timeout=60)
    local = os.path.join(RES, f"{sysid}.tgz")
    rc, o, e = C.download(f"/content/{sysid}.tgz", local)
    if rc == 0 and os.path.exists(local):
        with tarfile.open(local) as t: t.extractall(RES)
        log(f"{sysid}: downloaded+extracted"); return True
    log(f"{sysid}: download failed"); return False

def have_local(sysid):
    return os.path.exists(os.path.join(RES, sysid, f"{sysid}_summary.json"))

def run_one(sysid):
    if have_local(sysid): log(f"{sysid}: already local"); return True
    ensure_session()
    # launch only if not already done/running on VM
    try:
        st, _ = vm_state(sysid)
    except ConnectionError:
        st = "absent"
    if st == "DONE": return fetch(sysid)
    if st != "run":
        launch(sysid)
    streak = 0
    t0 = time.time()
    while time.time() - t0 < 90*60:
        time.sleep(30)
        try:
            st, tail = vm_state(sysid)
            streak = 0
            log(f"  {sysid}: {st} | {tail[:120]}")
            if st == "DONE": return fetch(sysid)
            if st == "FAIL": log(f"{sysid}: driver FAIL"); return False
            if st == "absent":
                # process gone but no DONE/FAIL -> likely killed; relaunch
                log(f"{sysid}: job absent w/o DONE — relaunching"); launch(sysid)
        except ConnectionError:
            streak += 1
            log(f"  {sysid}: transient drop (streak {streak})")
            if streak >= 5 and session_gone():
                log(f"{sysid}: session confirmed gone — recreating"); ensure_session(); launch(sysid); streak = 0
    return False

def main():
    log(f"=== tolerant runner: systems {SYS}, NS={NS} ===")
    for s in SYS:
        for attempt in range(3):
            if run_one(s): break
            log(f"{s}: attempt {attempt} failed; retrying")
    done = [s for s in "ABC" if have_local(s)]
    log(f"DONE_LOCAL={done}")
    if all(have_local(s) for s in "ABC"):
        try: C.stop(); log("session stopped (all done)")
        except Exception: pass
    print(f"RUN_TOLERANT_COMPLETE done={''.join(done)}", flush=True)

if __name__ == "__main__":
    main()
