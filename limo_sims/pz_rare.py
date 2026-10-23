#!/usr/bin/env python3
"""Paper fig pzdet / pzW: Pz bump vs n planted.

    python limo_sims/pz_rare.py
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
from plant import add_bump
from pool import load_mixed_pad
from tmaps import ESTIMATORS, estimator_boot, groups_from_arms

N_A_GRID = tuple(range(30, 1, -2)) + (1,)
N_B = 30
N_MC = 100
NBOOT = 101
TAU = 3.0
MINNBCHAN = 0
YUEN_MIN_N = 6
CHAN_NAME = "Pz"
AMP_JSON = HERE / "results" / "pz_amp_pilot" / "pz_amp_pilot.json"
OUT = HERE / "results" / "pz_rare"
CSV_PATH = OUT / "pz_rare.csv"
JSON_PATH = OUT / "pz_rare.json"

EST_LABEL = {"mean": "mean", "yuen": "Yuen 20%", "weighted": "weighted"}
EST_STYLE = {
    "mean": dict(color="#1f77b4", marker="o", ls="-"),
    "yuen": dict(color="#d55e00", marker="s", ls="--"),
    "weighted": dict(color="#009e73", marker="^", ls="-."),
}


def load_amp():
    rec = json.loads(AMP_JSON.read_text())
    return float(rec["recommended_amp"])


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
            r["amp"] = float(r["amp"])
            r["n_sig_Pz"] = int(float(r["n_sig_Pz"]))
            r["t_Pz"] = float(r["t_Pz"])
            r["mean_W_Pz_A"] = float(r["mean_W_Pz_A"])
            r["mean_W_Pz_B"] = float(r["mean_W_Pz_B"])
            r["gae_seconds"] = float(r.get("gae_seconds") or 0.0)
            rows.append(r)
            have.add((r["n_A"], r["mc"]))
    return rows, have


def one_mc(datain, Y, A, pz, n_a, mc, amp):
    n_sub = datain.shape[3]
    rng = np.random.default_rng(6000 + 1000 * n_a + mc)
    pick = rng.choice(n_sub, n_a + N_B, replace=False)
    control, planted = pick[:N_B], pick[N_B:]
    keep = np.sort(np.concatenate([planted, control]))
    remap = {int(i): k for k, i in enumerate(keep)}
    planted_k = np.array([remap[int(i)] for i in planted])
    control_k = np.array([remap[int(i)] for i in control])
    d, y = add_bump(datain[:, :, :, keep], Y[:, :, keep], planted_k, [pz], amp)
    t0 = time.time()
    W = compute_weights(d, A, beta="all", tau=TAU, seeds=(0,))
    w_sec = time.time() - t0
    groups = groups_from_arms(d.shape[3], planted_k, control_k)
    names = [e for e in ESTIMATORS if e != "yuen" or n_a >= YUEN_MIN_N]
    maps = {
        name: estimator_boot(
            name, y, groups, W, NBOOT,
            np.random.default_rng(7000 + 100 * n_a + mc + i),
        )
        for i, name in enumerate(names)
    }
    cluster = cluster_maps(maps, A, [pz], minnbchan=MINNBCHAN)
    w_a = float(W[pz][planted_k].mean())
    w_b = float(W[pz][control_k].mean())
    rows = []
    for name in names:
        c = cluster[name]
        tsl = maps[name]["t"][pz]
        rows.append(dict(
            amp=float(amp), n_A=int(n_a), n_B=N_B, mc=int(mc), estimator=name,
            sig=int(c["planted_any_cluster"]),
            n_sig_Pz=int(c["planted_n_sig"]),
            t_Pz=float(np.nanmax(np.abs(tsl))),
            mean_W_Pz_A=w_a,
            mean_W_Pz_B=w_b,
            gae_seconds=w_sec,
        ))
    return rows


def plot_detect(buckets, amp):
    detect_n = [n for n in N_A_GRID if n >= YUEN_MIN_N]
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    for est in ESTIMATORS:
        xs, y = [], []
        for n in detect_n:
            block = buckets[n].get(est, [])
            if not block:
                continue
            xs.append(n)
            y.append(np.mean([int(r["sig"]) for r in block]))
        if xs:
            ax.plot(xs, y, label=EST_LABEL[est], **EST_STYLE[est], lw=1.8, ms=7)
    ax.set_xlim(31.5, 5.5)
    ax.set_ylim(-0.05, 1.08)
    ax.set_xticks(detect_n)
    ax.axhline(0.05, color="0.6", lw=0.8, ls=":")
    ax.set_xlabel("n planted (30 unplanted)")
    ax.set_ylabel("LIMO FWER detection")
    ax.set_title("Detection (Yuen valid n≥6)")
    ax.legend(frameon=False)
    path = OUT / "gae_nA_detect.png"
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_W(wrows, amp):
    ns = list(N_A_GRID)
    y_a = np.array([np.mean([r["mean_W_Pz_A"] for r in wrows[n]]) for n in ns])
    y_b = np.array([np.mean([r["mean_W_Pz_B"] for r in wrows[n]]) for n in ns])
    sd_a = np.array([
        np.std([r["mean_W_Pz_A"] for r in wrows[n]], ddof=1) if len(wrows[n]) > 1 else 0.0
        for n in ns
    ])
    sd_b = np.array([
        np.std([r["mean_W_Pz_B"] for r in wrows[n]], ddof=1) if len(wrows[n]) > 1 else 0.0
        for n in ns
    ])
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    ax.fill_between(ns, y_a - sd_a, y_a + sd_a, color="#1f77b4", alpha=0.15, lw=0)
    ax.fill_between(ns, y_b - sd_b, y_b + sd_b, color="#d55e00", alpha=0.15, lw=0)
    ax.plot(ns, y_a, "o-", color="#1f77b4", lw=1.8, ms=7, label="Pz planted")
    ax.plot(ns, y_b, "s--", color="#d55e00", lw=1.8, ms=7, label="Pz unplanted (n=30)")
    ax.axhline(1.0, color="0.6", lw=0.8, ls=":")
    ax.set_xlim(31.5, 0.5)
    hi = max(
        1.05,
        float(np.nanmax(y_a + sd_a)),
        float(np.nanmax(y_b + sd_b)),
    )
    ax.set_ylim(-0.05, hi + 0.08)
    ax.set_xticks(ns)
    ax.set_xlabel("n planted (30 unplanted)")
    ax.set_ylabel("mean W on Pz")
    nmc = max(len(wrows[n]) for n in ns)
    ax.set_title(f"Mean W on Pz (AMP={amp:.2f}, {nmc} MC, band = mean ± SD)")
    ax.legend(frameon=False)
    path = OUT / "gae_nA_W.png"
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return path


def main():
    amp = load_amp()
    cohort = load_mixed_pad()
    datain, Y, A = cohort["datain"], cohort["Y"], cohort["adjacency"]
    pz = cohort["labels"].index(CHAN_NAME)
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    fields = [
        "amp", "n_A", "n_B", "mc", "estimator", "sig", "n_sig_Pz", "t_Pz",
        "mean_W_Pz_A", "mean_W_Pz_B", "gae_seconds",
    ]
    rows, have = load_rows()
    if rows and abs(float(rows[0]["amp"]) - amp) > 0.05:
        raise SystemExit(
            f"CSV AMP={rows[0]['amp']} != floor AMP={amp:.4f}; "
            f"move {CSV_PATH.parent} aside and rerun"
        )
    print(f"resume {len(have)}/{len(N_A_GRID) * N_MC} MCs  AMP={amp:.2f}", flush=True)
    new_file = not CSV_PATH.is_file()
    with CSV_PATH.open("a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        if new_file:
            writer.writeheader()
        for n_a in N_A_GRID:
            n_new = 0
            for mc in range(N_MC):
                if (n_a, mc) in have:
                    continue
                chunk = one_mc(datain, Y, A, pz, n_a, mc, amp)
                for row in chunk:
                    writer.writerow(row)
                    rows.append(row)
                have.add((n_a, mc))
                n_new += 1
                f.flush()
            hit = defaultdict(list)
            for r in rows:
                if r["n_A"] == n_a:
                    hit[r["estimator"]].append(r["sig"])
            w = np.mean([r["mean_W_Pz_A"] for r in rows if r["n_A"] == n_a and r["estimator"] == "mean"])
            det = "  ".join(f"{k} {np.mean(v):.2f}" for k, v in hit.items())
            print(f"n_A={n_a:2d}  +{n_new}  {det}  W_Pz={w:.3f}", flush=True)

    by = defaultdict(lambda: defaultdict(list))
    w_by = defaultdict(list)
    for r in rows:
        by[r["n_A"]][r["estimator"]].append(r)
        if r["estimator"] == "mean":
            w_by[r["n_A"]].append(r)
    summary = {}
    for n_a in N_A_GRID:
        wa = np.array([x["mean_W_Pz_A"] for x in w_by[n_a]], dtype=float)
        wb = np.array([x["mean_W_Pz_B"] for x in w_by[n_a]], dtype=float)
        summary[str(n_a)] = {
            "detect": {
                est: float(np.mean([x["sig"] for x in by[n_a][est]]))
                for est in ESTIMATORS if by[n_a][est]
            },
            "mean_W_Pz_A": float(np.mean(wa)),
            "mean_W_Pz_A_sd": float(np.std(wa, ddof=1)) if len(wa) > 1 else 0.0,
            "mean_W_Pz_B": float(np.mean(wb)),
            "mean_W_Pz_B_sd": float(np.std(wb, ddof=1)) if len(wb) > 1 else 0.0,
        }
    report = dict(
        design=(
            f"Pz bump AMP={amp:.2f} (10/10 mean FWER bisection from pz_amp_pilot), "
            f"{N_B} unplanted, n planted {list(N_A_GRID)}. W refit each MC. "
            f"Python t; LIMO FWER minnbchan={MINNBCHAN} on Pz. "
            f"{N_MC} MCs, nboot={NBOOT}. Yuen omitted for n<{YUEN_MIN_N}."
        ),
        amp=amp,
        n_B=N_B,
        n_A_grid=list(N_A_GRID),
        n_mc=N_MC,
        nboot=NBOOT,
        minnbchan=MINNBCHAN,
        summary=summary,
        seconds=time.time() - t0,
    )
    JSON_PATH.write_text(json.dumps(report, indent=2) + "\n")
    plot_detect(by, amp)
    plot_W(w_by, amp)
    print(json.dumps(summary, indent=2))
    print(f"wrote {CSV_PATH}")
    print(f"wrote {JSON_PATH}")
    print(f"wrote {OUT / 'gae_nA_detect.png'}")
    print(f"wrote {OUT / 'gae_nA_W.png'}")


if __name__ == "__main__":
    main()
