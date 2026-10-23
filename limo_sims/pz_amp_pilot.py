#!/usr/bin/env python3
"""Pz AMP detect-floor: 10/10 mean LIMO FWER, bisection, no GAE.

    python limo_sims/pz_amp_pilot.py
"""
from __future__ import annotations

import json
import os
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

HERE = Path(__file__).resolve().parent
import _paths  # noqa: F401

from match_limo import MATCH_SEEDS, N_HIT, N_MATCH_MC, detect_floor, frozen_split
from pool import load_mixed_pad

N_A = 5
N_B = 20
MINNBCHAN = 0
CHAN_NAME = "Pz"
OUT = HERE / "results" / "pz_amp_pilot"
JSON_PATH = OUT / "pz_amp_pilot.json"
PNG_PATH = OUT / "focal_pz_amp_pilot.png"


def main():
    t0 = time.time()
    cohort = load_mixed_pad()
    names = [str(x).strip() for x in cohort["labels"]]
    pz = [names.index(CHAN_NAME)]
    ctxs = [
        frozen_split(cohort, N_A, N_B, seed, minnbchan=MINNBCHAN, verbose=(seed == MATCH_SEEDS[0]))
        for seed in MATCH_SEEDS
    ]
    print(f"floor Pz  ({N_HIT}/{N_MATCH_MC} mean FWER)", flush=True)
    a0, rec, mass_of = detect_floor(ctxs, pz)
    cache = dict(sorted(mass_of.cache.items()))
    OUT.mkdir(parents=True, exist_ok=True)
    report = dict(
        design=(
            f"{N_B} vs {N_A} Pz-only Gaussian, mean Welch t, no GAE. "
            f"Floor is mean LIMO FWER (minnbchan={MINNBCHAN}) on {N_MATCH_MC} "
            f"calibration splits (seeds {list(MATCH_SEEDS)}): smallest amp with "
            f"{N_HIT}/{N_MATCH_MC} hits, bisection in [0, 256]."
        ),
        recommended_amp=float(a0),
        n_hit=int(rec["n_hit"]),
        n_mc=int(rec["n_mc"]),
        mass=float(rec["mass"]),
        max_th=float(rec["max_th"]),
        minnbchan=MINNBCHAN,
        n_A=N_A,
        n_B=N_B,
        match_seeds=list(MATCH_SEEDS),
        evals=[
            dict(amp=float(a), n_hit=int(v["n_hit"]), hit_rate=float(v["hit_rate"]),
                 mass=float(v["mass"]))
            for a, v in cache.items()
        ],
        seconds=time.time() - t0,
    )
    JSON_PATH.write_text(json.dumps(report, indent=2) + "\n")

    amps = [float(a) for a in cache]
    hits = [cache[a]["hit_rate"] for a in cache]
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    ax.plot(amps, hits, "o-", color="#1f77b4", lw=1.6, ms=6)
    ax.axvline(a0, color="0.4", ls=":", lw=1.2)
    ax.axhline(1.0, color="0.75", ls="--", lw=0.8)
    ax.set_xlabel("bump amplitude")
    ax.set_ylabel("mean LIMO FWER (10 draws)")
    ax.set_ylim(-0.05, 1.08)
    ax.set_title(f"Pz detect-floor  20 vs 5  a={a0:.2f}")
    fig.savefig(PNG_PATH, dpi=160, bbox_inches="tight")
    plt.close(fig)

    print(json.dumps({"recommended_amp": a0, "n_hit": rec["n_hit"]}, indent=2))
    print(f"wrote {JSON_PATH}")
    print(f"wrote {PNG_PATH}")


if __name__ == "__main__":
    main()
