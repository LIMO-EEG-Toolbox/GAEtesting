#!/usr/bin/env python3
"""Paper fig g1detect / g1W: whole-scalp 40 ms delay, LIMO FWER.

    python limo_sims/g1_latency.py
"""
from __future__ import annotations

import csv
import json
import os
import time
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

HERE = Path(__file__).resolve().parent
import _paths  # noqa: F401

from gae_weighting import compute_weights

from cluster_python import cluster_maps
from pool import load_mixed_pad
from tmaps import ESTIMATORS, all_estimator_boots, groups_from_arms

N_A_GRID = (20, 15, 10, 6)
N_B = 20
N_MC = 100
NBOOT = 101
SHIFT_MS = 40.0
DT_MS = 4.0
TAU = 3.0
MINNBCHAN = 2
YUEN_MIN_N = 6
OUT = HERE / "results" / "g1_latency"
CSV_PATH = OUT / "g1_latency.csv"
JSON_PATH = OUT / "g1_latency.json"

EST_LABEL = {"mean": "mean", "yuen": "Yuen 20%", "weighted": "weighted"}
EST_STYLE = {
    "mean": dict(color="#1f77b4", marker="o", ls="-"),
    "yuen": dict(color="#d55e00", marker="s", ls="--"),
    "weighted": dict(color="#009e73", marker="^", ls="-."),
}


def shift_Y(Y, n_samp, who):
    out = Y.copy()
    out[n_samp:, :, who] = Y[:-n_samp, :, who]
    out[:n_samp, :, who] = 0
    return out


def shift_datain(datain, n_samp, who):
    out = datain.copy()
    out[:, n_samp:, :, who] = datain[:, :-n_samp, :, who]
    out[:, :n_samp, :, who] = 0
    return out


def one_mc(datain, Y, A, n_a, mc):
    n_sub = datain.shape[3]
    n_samp = int(round(SHIFT_MS / DT_MS))
    rng = np.random.default_rng(2000 + 1000 * n_a + mc)
    pick = rng.choice(n_sub, n_a + N_B, replace=False)
    control, planted = pick[:N_B], pick[N_B:]
    keep = np.sort(np.concatenate([planted, control]))
    remap = {int(i): k for k, i in enumerate(keep)}
    planted_k = np.array([remap[int(i)] for i in planted])
    control_k = np.array([remap[int(i)] for i in control])
    d = shift_datain(datain[:, :, :, keep], n_samp, planted_k)
    y = shift_Y(Y[:, :, keep], n_samp, planted_k)
    t0 = time.time()
    W = compute_weights(d, A, beta="all", tau=TAU, seeds=(0,))
    w_sec = time.time() - t0
    groups = groups_from_arms(d.shape[3], planted_k, control_k)
    maps = all_estimator_boots(y, groups, W, NBOOT, 3000 + 100 * n_a + mc)
    names = [e for e in ESTIMATORS if e != "yuen" or n_a >= YUEN_MIN_N]
    all_ch = np.arange(y.shape[1], dtype=int)
    cluster = cluster_maps({k: maps[k] for k in names}, A, all_ch, minnbchan=MINNBCHAN)
    WA, WB = W[:, planted_k], W[:, control_k]
    rows = []
    for name in names:
        c = cluster[name]
        rows.append(dict(
            n_A=int(n_a), n_B=N_B, mc=int(mc), estimator=name,
            sig=int(c["planted_any_cluster"]),
            n_sig_clusters=int(c["n_sig_clusters"]),
            n_sig_cells=int(c["n_sig_cells"]),
            max_th=float(c["max_th"]),
            min_cluster_p=("" if c["min_cluster_p"] is None else float(c["min_cluster_p"])),
            max_abs_t=float(np.nanmax(np.abs(maps[name]["t"]))),
            mean_W_A=float(WA.mean()),
            mean_W_B=float(WB.mean()),
            gae_seconds=w_sec,
        ))
    return rows


def load_rows():
    if not CSV_PATH.is_file():
        return [], set()
    rows = []
    have = set()
    with CSV_PATH.open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            r["n_A"] = int(r["n_A"])
            r["n_B"] = int(r["n_B"])
            r["mc"] = int(r["mc"])
            r["sig"] = int(r["sig"])
            r["mean_W_A"] = float(r["mean_W_A"])
            r["mean_W_B"] = float(r["mean_W_B"])
            rows.append(r)
            have.add((r["n_A"], r["mc"]))
    return rows, have


def dump_rows(rows, fields):
    with CSV_PATH.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def plot_detect(buckets):
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    x = np.arange(len(N_A_GRID))
    for est in ESTIMATORS:
        xs, y = [], []
        for i, n in enumerate(N_A_GRID):
            block = buckets[n].get(est, [])
            if not block:
                continue
            xs.append(i)
            y.append(np.mean([int(r["sig"]) for r in block]))
        if xs:
            ax.plot(xs, y, label=EST_LABEL[est], **EST_STYLE[est], lw=1.8, ms=8)
    ax.set_xticks(x, [str(n) for n in N_A_GRID])
    ax.set_ylim(-0.05, 1.08)
    ax.set_xlabel("n delayed (Arm A)")
    ax.set_ylabel("LIMO FWER detection")
    ax.set_title("Power vs delayed-arm size")
    ax.axhline(0.05, color="0.6", lw=0.8, ls=":")
    ax.legend(frameon=False)
    path = OUT / "g1_detect.png"
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_weights(wrows):
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    x = np.arange(len(N_A_GRID))
    wa = [np.mean([r["mean_W_A"] for r in wrows[n]]) for n in N_A_GRID]
    wb = [np.mean([r["mean_W_B"] for r in wrows[n]]) for n in N_A_GRID]
    ax.plot(x, wa, "o-", color="#1f77b4", lw=1.8, ms=7, label="Arm A (delayed)")
    ax.plot(x, wb, "s--", color="0.45", lw=1.6, ms=7, label="Arm B (unshifted)")
    ax.axhline(1.0, color="0.6", lw=0.8, ls=":")
    ax.set_xticks(x, [str(n) for n in N_A_GRID])
    lo = min(0.985, min(wa + wb) - 0.01)
    hi = max(1.015, max(wa + wb) + 0.01)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("n delayed (Arm A)")
    ax.set_ylabel("mean W")
    ax.set_title("Mean GAE weight (channels × people)")
    ax.legend(frameon=False)
    path = OUT / "g1_weights.png"
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return path


def main():
    cohort = load_mixed_pad()
    datain, Y, A = cohort["datain"], cohort["Y"], cohort["adjacency"]
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    rows, have = load_rows()
    fields = [
        "n_A", "n_B", "mc", "estimator", "sig", "n_sig_clusters", "n_sig_cells",
        "max_th", "min_cluster_p", "max_abs_t", "mean_W_A", "mean_W_B", "gae_seconds",
    ]
    print(f"resume {len(have)} unique (n_A, mc)", flush=True)
    n_new = 0
    for n_a in N_A_GRID:
        for mc in range(N_MC):
            if (n_a, mc) in have:
                continue
            chunk = one_mc(datain, Y, A, n_a, mc)
            rows.extend(chunk)
            have.add((n_a, mc))
            n_new += 1
            hit = "  ".join(
                f"{r['estimator']} {r['sig']}" for r in chunk
            )
            w = chunk[0]["mean_W_A"]
            print(f"n_A={n_a} mc={mc}  {hit}  W_A={w:.3f}", flush=True)
            if n_new % 10 == 0:
                dump_rows(rows, fields)
    dump_rows(rows, fields)

    by = defaultdict(lambda: defaultdict(list))
    w_by = defaultdict(list)
    for r in rows:
        by[r["n_A"]][r["estimator"]].append(r)
        if r["estimator"] == "mean":
            w_by[r["n_A"]].append(r)
    summary = {}
    for n_a in N_A_GRID:
        summary[str(n_a)] = {
            "detect": {
                est: float(np.mean([r["sig"] for r in by[n_a][est]]))
                for est in ESTIMATORS if by[n_a][est]
            },
            "mean_W_A": float(np.mean([r["mean_W_A"] for r in w_by[n_a]])),
            "mean_W_B": float(np.mean([r["mean_W_B"] for r in w_by[n_a]])),
        }
    report = dict(
        design=(
            f"Unbalanced whole-scalp 40 ms delay on n_A vs {N_B} unshifted, mixed pad. "
            f"W refit each MC. Python t; LIMO cluster-mass FWER minnbchan={MINNBCHAN} "
            f"on all channels. {N_MC} MCs, nboot={NBOOT}."
        ),
        shift_ms=SHIFT_MS,
        n_mc=N_MC,
        nboot=NBOOT,
        minnbchan=MINNBCHAN,
        n_B=N_B,
        n_A_grid=list(N_A_GRID),
        summary=summary,
        seconds=time.time() - t0,
    )
    JSON_PATH.write_text(json.dumps(report, indent=2) + "\n")
    plot_detect(by)
    plot_weights(w_by)
    print(json.dumps(summary, indent=2))
    print(f"wrote {CSV_PATH}")
    print(f"wrote {JSON_PATH}")
    print(f"wrote {OUT / 'g1_detect.png'}")
    print(f"wrote {OUT / 'g1_weights.png'}")


if __name__ == "__main__":
    main()
