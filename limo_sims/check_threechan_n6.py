#!/usr/bin/env python3
"""One Monte Carlo: three-channel STUDY triangle, 6 vs 30, LIMO cluster.

    python limo_sims/check_threechan_n6.py
"""
from __future__ import annotations

import json

from run_mc import HERE, run_one

N_PLANTED = 3
N_A = 6
N_B = 30
SPLIT_SEED = 0
AMP = 64.0
NBOOT = 101


def main():
    report, path = run_one(
        n_planted=N_PLANTED,
        n_a=N_A,
        n_b=N_B,
        amp=AMP,
        nboot=NBOOT,
        split_seed=SPLIT_SEED,
        work=HERE / "work" / "threechan_n6",
        out=HERE / "results" / "check_threechan_n6",
    )
    keys = (
        "clusterer", "channels", "mean_W_planted", "mean_W_control",
        "matlab", "cluster",
    )
    print(json.dumps({k: report[k] for k in keys}, indent=2, default=str))
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
