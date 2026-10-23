#!/usr/bin/env python3
"""Paper matcher display: 1 vs 8, 30 vs 5, minnbchan=0.

    python limo_sims/match_1v8_a5.py
"""
from __future__ import annotations

import json

from match_limo import HERE, match_mass_amps, plot_sig_maps

N1 = 1
N2 = 8
N_A = 5
N_B = 30
MINNBCHAN = 0
OUT = HERE / "results" / "match_1v8_a5"


def main():
    rec = match_mass_amps(n1=N1, n2=N2, n_a=N_A, n_b=N_B, minnbchan=MINNBCHAN)
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "match.json"
    path.write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps({
        "minnbchan": rec["minnbchan"],
        "n_A": rec["n_A"],
        "a1": rec["a1"],
        "a2": rec["a2"],
        "matched_mass": rec["matched_mass"],
        "n1": rec["region1"]["names"],
        "n2": rec["region2"]["names"],
        "fwer1": rec["region1"]["fwer_hit"],
        "fwer2": rec["region2"]["fwer_hit"],
        "floors": rec["floors"],
    }, indent=2))
    print(f"wrote {path}")
    plot_sig_maps(rec, OUT)


if __name__ == "__main__":
    main()
