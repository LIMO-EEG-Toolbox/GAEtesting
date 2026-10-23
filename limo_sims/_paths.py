"""Folders that travel with this package.

Copy `limo_sims/` into another repo (e.g. LIMO-EEG-Toolbox/GAEtesting) as a
unit. Maps live in `data/P3_limo/` *inside this folder*. Do not edit
gae_weighting.py.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def _first_dir(cands, *, marker=None):
    leftover = None
    for p in cands:
        if p is None:
            continue
        p = Path(p)
        if leftover is None:
            leftover = p
        hit = (p / marker).is_file() if marker else p.is_dir()
        if hit:
            return p
    return leftover


DATA = _first_dir(
    (
        HERE / "data" / "P3_limo",
        ROOT / "data" / "P3_limo",  # old layout, this working copy
    ),
    marker="Beta_files.txt",
)

env_tools = os.environ.get("LIMO_TOOLS")
LIMO_TOOLS = _first_dir(
    (
        Path(env_tools) if env_tools else None,
        ROOT / "limo_tools",
        HERE / "limo_tools",
    ),
)

POOL = next(
    (
        p
        for p in (
            HERE / "cache" / "mixed_pool.npz",
            ROOT / "gae_audit" / "cache" / "mixed_pool.npz",
        )
        if p.is_file()
    ),
    HERE / "cache" / "mixed_pool.npz",
)

for p in (HERE,):
    s = str(p)
    if s not in sys.path:
        sys.path.insert(0, s)
