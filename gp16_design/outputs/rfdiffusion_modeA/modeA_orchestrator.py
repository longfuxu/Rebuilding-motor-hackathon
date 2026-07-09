#!/usr/bin/env python3
"""Reliable Mode-A orchestrator: partial diffusion to rigidify cp233 (1750 aa).
Feasibility (1 design, partial_T=10) first to check 40GB OOM; if it produces a
backbone, run production (partial_T 8/12/20). Hard local timeouts, tarball capture,
immediate colab stop."""
import subprocess, time, os, json

COLAB = "/Users/longfu/miniforge3/bin/colab"
S = "rfdA"
TMP = "/Users/longfu/.claude/jobs/48b7a543/tmp"
DL = f"{TMP}/dl_modeA"; os.makedirs(DL, exist_ok=True)
CIF = ("/Users/longfu/Developer/claude-science-hackthon/.claude/worktrees/export-gp16-results/"
       "gp16_design/outputs/structures/af3_sweep/cp233/inter10/fold_2026_07_08_cp233_int15_inter10_model_0.cif")


def run(args, timeout, inp=None, retries=2):
    for k in range(retries + 1):
        try:
            r = subprocess.run([COLAB] + args, input=inp, capture_output=True, text=True, timeout=timeout)
            return r.returncode, r.stdout, r.stderr
        except subprocess.TimeoutExpired:
            print(f"  [timeout {timeout}s {args[0]} try {k}]", flush=True)
    return 124, "", "hard-timeout"


def cexec(code, timeout=60):
    return run(["exec", "-s", S, "--timeout", str(max(10, timeout - 8))], timeout, inp=code)[1]


def stop():
    print("stopping session:", run(["stop", "-s", S], 60)[0], flush=True)


def newsession():
    run(["stop", "-s", S], 40)
    for a in range(3):
        rc, out, err = run(["new", "-s", S, "--gpu", "A100"], 180)
        if "READY" in out:
            print("session READY", flush=True); return True
        print(f"[new try {a}] {out.strip()[-80:]}", flush=True); time.sleep(15)
    return False


def launch_and_wait(spec_remote, tag, max_polls=40, poll_s=30):
    """Launch driver on spec, poll to DONE/DEAD, return progress dict."""
    cexec(f'''import os,subprocess,time
os.chdir("/content")
subprocess.run("pkill -9 -f rfdiff_driver.py 2>/dev/null",shell=True);time.sleep(2)
for f in ["{tag}.DONE","{tag}.progress.json","{tag}.rfinf.log","{tag}.driver.log","{tag}.pid","{tag}.lock"]:
    try: os.remove("/content/"+f)
    except: pass
os.system("cd /content && nohup python -u rfdiff_driver.py {spec_remote} > {tag}.driver.log 2>&1 & echo $!>{tag}.pid")
time.sleep(6)
print("LAUNCH", subprocess.run("pgrep -af \\"[p]ython.*rfdiff_driver\\"",shell=True,capture_output=True,text=True).stdout.strip() or "NONE")''', timeout=45)
    for i in range(max_polls):
        out = cexec(f'''import os,glob,subprocess
done=os.path.isfile("/content/{tag}.DONE")
alive=subprocess.run("pgrep -fc \\"[p]ython.*rfdiff_driver\\"",shell=True,capture_output=True,text=True).stdout.strip()
npdb=len(glob.glob("/content/outputs/{tag.replace('.json','')}*_*.pdb"))+len(glob.glob("/content/outputs/*"+"{tag}"[:5]+"*_*.pdb"))
oom=""
if os.path.isfile("/content/{tag}.rfinf.log"):
    t=open("/content/{tag}.rfinf.log").read()
    if "out of memory" in t.lower() or "CUDA" in t and "memory" in t.lower(): oom="OOM"
print("PROG",("DONE" if done else ("DEAD" if alive=="0" else "RUN")),"pdbs="+str(npdb),oom)''', timeout=45)
        print(f"  {tag} poll {i}: {out.strip()[:70] or '<timeout>'}", flush=True)
        if "DONE" in out or "DEAD" in out:
            break
        time.sleep(poll_s)
    # fetch progress + rfinf tail for diagnosis
    prog = cexec(f'import os;print(open("/content/{tag}.progress.json").read() if os.path.isfile("/content/{tag}.progress.json") else "{{}}")', timeout=40)
    tail = cexec(f'import os;print(open("/content/{tag}.rfinf.log").read()[-800:] if os.path.isfile("/content/{tag}.rfinf.log") else "")', timeout=40)
    return prog, tail


def capture(tag):
    cexec(f'import subprocess;subprocess.run("cd /content && tar czf {tag}_results.tgz outputs/*{tag[:5]}* {tag}.progress.json 2>/dev/null",shell=True);print("tar")', timeout=60)
    rc, out, err = run(["download", "-s", S, f"/content/{tag}_results.tgz", f"{DL}/{tag}_results.tgz"], 120, retries=3)
    if os.path.isfile(f"{DL}/{tag}_results.tgz"):
        subprocess.run(f"cd {DL} && tar xzf {tag}_results.tgz", shell=True)
    print(f"captured {tag} rc={rc}", flush=True)


def main():
    t0 = time.time()
    if not newsession():
        print("MODEA_ABORT no session", flush=True); return
    # upload
    for f, r in [(f"{TMP}/rfdiff_setup.py", "/content/rfdiff_setup.py"),
                 (f"{TMP}/rfdiff_driver.py", "/content/rfdiff_driver.py"),
                 (CIF, "/content/cp233_inter10_model0.cif")]:
        print("upload", os.path.basename(f), run(["upload", "-s", S, f, r], 90)[0], flush=True)
    # specs
    feas = '{"tag":"modeAf","iterations":50,"jobs":[["modeAf_pT10","A1-1750",1,"10","/content/cp233_inter10.pdb","A"]]}'
    prod = ('{"tag":"modeAp","iterations":50,"jobs":['
            '["modeAp_pT8","A1-1750",2,"8","/content/cp233_inter10.pdb","A"],'
            '["modeAp_pT12","A1-1750",2,"12","/content/cp233_inter10.pdb","A"],'
            '["modeAp_pT20","A1-1750",2,"20","/content/cp233_inter10.pdb","A"]]}')
    open(f"{TMP}/modeAf.json", "w").write(feas); open(f"{TMP}/modeAp.json", "w").write(prod)
    for f in ["modeAf.json", "modeAp.json"]:
        run(["upload", "-s", S, f"{TMP}/{f}", f"/content/{f}"], 60)
    # install
    print("installing...", flush=True)
    cexec('import subprocess;print(subprocess.run("nohup python -u /content/rfdiff_setup.py > /content/setup.log 2>&1 &",shell=True))', timeout=40)
    ready = False
    for i in range(40):
        out = cexec('import os,sys\ntry:\n sys.path.append("RFdiffusion");from inference.utils import parse_pdb;import dgl\n print("OK",os.path.isdir("RFdiffusion/models") and len(os.listdir("RFdiffusion/models")))\nexcept Exception as e:print("NO",repr(e)[:50])', timeout=45)
        print(f"  install {i}: {out.strip()[:70]}", flush=True)
        if "OK" in out and "OK 0" not in out and "OK False" not in out:
            ready = True; break
        time.sleep(15)
    if not ready:
        print("INSTALL_FAIL", flush=True); stop(); return
    # convert cp233 CIF->PDB
    conv = cexec('import subprocess,os\nsubprocess.run("pip install -q gemmi",shell=True)\nimport gemmi\nst=gemmi.read_structure("/content/cp233_inter10_model0.cif");st.setup_entities()\nopen("/content/cp233_inter10.pdb","w").write(st.make_pdb_string())\nprint("PDB",os.path.getsize("/content/cp233_inter10.pdb"))', timeout=90)
    print("convert:", conv.strip()[-60:], flush=True)

    # feasibility
    print("=== FEASIBILITY partial_T=10, 1 design (OOM check) ===", flush=True)
    prog, tail = launch_and_wait("modeAf.json", "modeAf", max_polls=50, poll_s=30)
    print("feas progress:", prog.strip()[:400], flush=True)
    feas_ok = '"ok": true' in prog.lower() and '"n_pdbs": 0' not in prog
    capture("modeAf")
    if not feas_ok:
        print("FEASIBILITY_FAILED (likely OOM or error). tail:", tail.strip()[-400:], flush=True)
        print(f"MODEA_DONE feasible=False in {int(time.time()-t0)}s", flush=True); stop(); return
    print("FEASIBILITY OK -> production", flush=True)

    # production
    prog2, tail2 = launch_and_wait("modeAp.json", "modeAp", max_polls=80, poll_s=40)
    print("prod progress:", prog2.strip()[:600], flush=True)
    capture("modeAp")
    print(f"MODEA_DONE feasible=True in {int(time.time()-t0)}s", flush=True)
    stop()


if __name__ == "__main__":
    main()
