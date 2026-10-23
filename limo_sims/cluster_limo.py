"""LIMO spatial clustering on Python t^2 maps (MCC=2, minnbchan=2)."""
from __future__ import annotations

import json
import shutil
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
from _paths import DATA, LIMO_TOOLS

CHANLOCS = DATA / "limo_chanlocs.mat"


def find_matlab():
    exe = shutil.which("matlab")
    if exe:
        return exe
    apps = Path("/Applications")
    if apps.is_dir():
        for p in sorted(apps.glob("MATLAB_*.app/bin/matlab")):
            return str(p)
        p = apps / "MATLAB.app" / "bin" / "matlab"
        if p.is_file():
            return str(p)
    return None


def run_matlab(work):
    matlab = find_matlab()
    if matlab is None:
        return {"ok": False, "error": "matlab_not_on_path"}
    cmd = [
        matlab, "-batch",
        f"addpath('{HERE}'); run_limo_ttest('{work}', '{LIMO_TOOLS}', '{CHANLOCS}');",
    ]
    t0 = time.time()
    proc = subprocess.run(cmd, cwd=str(work), capture_output=True, text=True)
    return {
        "ok": proc.returncode == 0,
        "error": None if proc.returncode == 0 else (proc.stderr or proc.stdout)[-4000:],
        "seconds": time.time() - t0,
        "cmd": " ".join(cmd),
    }


def cluster_maps(work):
    matlab = run_matlab(work)
    cluster = {}
    js = Path(work) / "limo_cluster_summary.json"
    if js.is_file():
        cluster = json.loads(js.read_text())
    return matlab, cluster
