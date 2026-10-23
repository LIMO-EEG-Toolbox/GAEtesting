#!/usr/bin/env python3
"""Joint isolation × rarity W grid (paper fig joint).

n1=1..5 vs 8; n planted in {5,6,7,8,11,14,17,20,25} vs 30.
LIMO-adjacent: match_limo 10/10 FWER then equal mass; Python t; cluster_test
(forming p<=0.05, mass t^2, boot max-mass T). minnbchan=0 if n1<3 else 2.
Reuses iso_bal match JSON when minnbchan matches.

Completed 10-MC grid is frozen in results/joint_grid/. This file is the
runner; do not bump N_MC without a new seed block.

    python limo_sims/joint_grid.py
"""
from __future__ import annotations

import json
import os
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

HERE = Path(__file__).resolve().parent
import _paths  # noqa: F401

from iso_bal import (
    BAL_MATCH,
    EVAL_SEED0,
    ISO_MATCH,
    load_mc,
    load_or_match,
    minnb_iso,
    one_mc,
    summarize_cell,
)
from pool import load_mixed_pad

N1_GRID = (1, 2, 3, 4, 5)
N2 = 8
NA_GRID = (5, 6, 7, 8, 11, 14, 17, 20, 25)
N_B = 30
N_MC = 10
OUT = HERE / "results" / "joint_grid"
MATCH = OUT / "match"
CSV_PATH = OUT / "joint_mc.csv"


def minnbchan_of(n1):
    return minnb_iso(n1)


def match_path(n1, n_a):
    name = f"match_n1{n1}_nA{n_a}.json"
    for d in (MATCH, BAL_MATCH, ISO_MATCH):
        p = d / name
        if p.is_file():
            rec = json.loads(p.read_text())
            if int(rec.get("minnbchan", -1)) == minnbchan_of(n1):
                return p
    return MATCH / name


def heatmap(ax, Z, title, vmin, vmax, cmap):
    im = ax.imshow(Z, origin="lower", aspect="auto", vmin=vmin, vmax=vmax, cmap=cmap)
    ax.set_xticks(np.arange(len(NA_GRID)), [str(n) for n in NA_GRID])
    ax.set_yticks(np.arange(len(N1_GRID)), [str(n) for n in N1_GRID])
    ax.set_xlabel("n planted")
    ax.set_ylabel("n1")
    ax.set_title(title)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    for i in range(len(N1_GRID)):
        for j in range(len(NA_GRID)):
            v = Z[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7,
                        color="white" if v < 0.45 * (vmin + vmax) else "black")
    return im


def plot_joint(summaries, path):
    w1 = np.full((len(N1_GRID), len(NA_GRID)), np.nan)
    sd = np.full_like(w1, np.nan)
    idx1 = {n: i for i, n in enumerate(N1_GRID)}
    idxa = {n: j for j, n in enumerate(NA_GRID)}
    for r in summaries:
        i, j = idx1[r["n1"]], idxa[r["n_A"]]
        w1[i, j] = r["mean_W_n1_planted"]
        sd[i, j] = r["mean_W_n1_planted_sd"]
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.4))
    heatmap(axes[0], w1, "mean W n1 (planted)", 0, 1.15, "viridis")
    heatmap(axes[1], sd, "SD W n1 (planted)", 0, max(0.05, float(np.nanmax(sd))), "magma")
    fig.suptitle(f"Joint isolation × rarity. Planted GAE W on n1  ({N_MC} MC)")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {path}")


def dump_rows(rows):
    if not rows:
        return
    import csv
    fields = list(rows[-1].keys())
    with CSV_PATH.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    MATCH.mkdir(parents=True, exist_ok=True)
    cohort = load_mixed_pad()
    rows, have = load_mc(CSV_PATH)
    summaries = []
    skipped = []
    for n1 in N1_GRID:
        k = minnbchan_of(n1)
        for n_a in NA_GRID:
            mpath = match_path(n1, n_a)
            try:
                rec = load_or_match(
                    mpath, n1=n1, n2=N2, n_a=n_a, n_b=N_B,
                    minnbchan=k, cohort=cohort,
                )
            except RuntimeError as e:
                print(f"SKIP n1={n1} n_A={n_a} minnbchan={k}: {e}", flush=True)
                skipped.append(dict(n1=n1, n_A=n_a, minnbchan=k, error=str(e)))
                continue
            n_new = 0
            for i in range(N_MC):
                seed = EVAL_SEED0 + i
                if (n1, n_a, seed) in have:
                    continue
                row = one_mc(cohort, rec, seed, k)
                rows.append(row)
                have.add((n1, n_a, seed))
                n_new += 1
                dump_rows(rows)
            block = [r for r in rows if r["n1"] == n1 and r["n_A"] == n_a]
            w1 = float(np.mean([r["mean_W_n1_planted"] for r in block]))
            print(f"joint n1={n1} n_A={n_a} minnb={k}  +{n_new}  W_n1={w1:.3f}", flush=True)
            summaries.append(summarize_cell(block, rec, k))
            (OUT / "joint.json").write_text(json.dumps({
                "n_mc": N_MC,
                "cells_done": len(summaries),
                "skipped": skipped,
                "summary": summaries,
                "seconds": time.time() - t0,
            }, indent=2) + "\n")
    report = dict(
        design=(
            f"Joint n1={list(N1_GRID)} vs {N2}, n planted {list(NA_GRID)} vs {N_B}. "
            f"match_limo 10/10 FWER then equal mass; Python t; LIMO cluster_test "
            f"(minnbchan=2 if n1>=3 else 0). {N_MC} held-out MCs. Mixed pad."
        ),
        n_mc=N_MC,
        skipped=skipped,
        summary=summaries,
        seconds=time.time() - t0,
    )
    (OUT / "joint.json").write_text(json.dumps(report, indent=2) + "\n")
    if summaries:
        plot_joint(summaries, OUT / "joint_W_n1.png")
    print(f"wrote {OUT / 'joint.json'}  cells={len(summaries)} skipped={len(skipped)}")


if __name__ == "__main__":
    main()
