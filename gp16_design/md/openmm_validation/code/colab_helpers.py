"""Robust Colab-CLI driver: every call has a hard subprocess timeout + retries.
Never wedges (memory: colab exec/download websockets intermittently hang past --timeout)."""
import subprocess, time, os, tempfile

COLAB = "/Users/longfu/miniforge3/bin/colab"
SESSION = os.environ.get("GP16_SESSION", "gp16md")

def _run(args, timeout, retries=2, quiet=False):
    last = None
    for att in range(retries+1):
        try:
            r = subprocess.run([COLAB]+args, capture_output=True, text=True, timeout=timeout)
            if not quiet:
                tag = " ".join(args[:2])
                print(f"  [colab {tag}] rc={r.returncode} ({att})", flush=True)
            return r.returncode, r.stdout, r.stderr
        except subprocess.TimeoutExpired as e:
            last = e
            print(f"  [colab {' '.join(args[:2])}] TIMEOUT {timeout}s (attempt {att})", flush=True)
            time.sleep(3)
    return 124, "", f"timeout after {retries+1} attempts: {last}"

def new_session(gpu="A100", timeout=240):
    return _run(["new","-s",SESSION,"--gpu",gpu], timeout=timeout, retries=1)

def sessions(timeout=40):
    return _run(["sessions"], timeout=timeout, retries=1)

def stop(timeout=60):
    return _run(["stop","-s",SESSION], timeout=timeout, retries=2)

def upload(local, remote, timeout=180, retries=3):
    return _run(["upload","-s",SESSION,local,remote], timeout=timeout, retries=retries)

def download(remote, local, timeout=300, retries=4):
    return _run(["download","-s",SESSION,remote,local], timeout=timeout, retries=retries)

def exec_code(code, timeout=60, retries=2):
    """Write code to a temp file, colab exec -f it. Returns (rc,out,err)."""
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False,
                                     dir="/Users/longfu/.claude/jobs/7babac4f/tmp") as fh:
        fh.write(code); path = fh.name
    try:
        return _run(["exec","-s",SESSION,"-f",path,"--timeout",str(timeout)],
                    timeout=timeout+30, retries=retries)
    finally:
        try: os.remove(path)
        except OSError: pass

def exec_file(path, timeout=60, retries=2):
    return _run(["exec","-s",SESSION,"-f",path,"--timeout",str(timeout)],
                timeout=timeout+30, retries=retries)

def poll_for(sentinel_code, ok_token, fail_token="__FAIL__", interval=25,
             max_min=90, tail_lines=6):
    """Repeatedly exec sentinel_code (prints status incl ok_token/fail_token) until seen."""
    t0 = time.time()
    while (time.time()-t0) < max_min*60:
        rc,out,err = exec_code(sentinel_code, timeout=45, retries=1, )
        blob = (out or "")+(err or "")
        for ln in blob.strip().splitlines()[-tail_lines:]:
            print("   |", ln, flush=True)
        if ok_token in blob:
            return True, blob
        if fail_token in blob:
            return False, blob
        time.sleep(interval)
    return False, "poll timeout"
