#!/usr/bin/env python3
"""Isolation and balance W margins.

Isolation: n1=1..5 vs 8, 30 vs 5.
  minnbchan=0 if n1<3, else 2 (matcher and test).
Balance: 1 vs 8, n planted in the paper grid vs 30, minnbchan=0.

    python limo_sims/iso_bal.py
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
from load_cohort import load_cohort

from cluster_python import cluster_maps
from match_limo import match_mass_amps
from plant import add_bump
from pool import load_mixed_pad
from tmaps import ESTIMATORS, groups_from_arms

ISO_N1 = (1, 2, 3, 4, 5)
ISO_N2 = 8
ISO_NA = 5
ISO_NB = 30
def minnb_iso(n1):
    return 2 if int(n1) >= 3 else 0

BAL_N1 = 1
BAL_N2 = 8
BAL_NA = (5, 6, 7, 8, 10, 11, 14, 17, 20, 25)
BAL_NB = 30
BAL_MINNB = 0

N_MC = 100
EVAL_SEED0 = 10
NBOOT = 101
TAU = 3.0
YUEN_MIN_N = 6
OUT = HERE / "results" / "iso_bal"
ISO_MATCH = OUT / "iso_match"
BAL_MATCH = OUT / "bal_match"
ISO_CSV = OUT / "iso_mc.csv"
BAL_CSV = OUT / "bal_mc.csv"

W_KEYS = (
    ("mean_W_n1_planted", "n1 planted", dict(color="#1f77b4", marker="o", ls="-")),
    ("mean_W_n2_planted", "n2 planted", dict(color="#d55e00", marker="s", ls="-")),
    ("mean_W_n1_unplanted", "n1 unplanted", dict(color="#1f77b4", marker="o", ls="--")),
    ("mean_W_n2_unplanted", "n2 unplanted", dict(color="#d55e00", marker="s", ls="--")),
)


def mean_W(W, chans, who):
    sl = W[np.asarray(chans, dtype=int)][:, np.asarray(who, dtype=int)]
    return float(np.mean(sl)) if sl.size else float("nan")


def load_or_match(path, **kw):
    want = int(kw.get("minnbchan", -1))
    if path.is_file():
        rec = json.loads(path.read_text())
        if int(rec.get("minnbchan", -1)) == want:
            print(f"match cache {path.name}  a1={rec['a1']:.2f} a2={rec['a2']:.2f}", flush=True)
            return rec
        print(f"match cache minnbchan mismatch {path.name}, rematch", flush=True)
    print(f"matching {path.name}  n1={kw.get('n1')} n_A={kw.get('n_a')} minnbchan={kw.get('minnbchan')}", flush=True)
    rec = match_mass_amps(**kw)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rec, indent=2) + "\n")
    print(f"  a1={rec['a1']:.2f} a2={rec['a2']:.2f} mass={rec['matched_mass']:.1f}", flush=True)
    return rec


def one_mc(cohort, rec, seed, minnbchan):
    s1 = np.asarray(rec["region1"]["channels"], dtype=int)
    s2 = np.asarray(rec["region2"]["channels"], dtype=int)
    a1, a2 = float(rec["a1"]), float(rec["a2"])
    n_a, n_b = int(rec["n_A"]), int(rec["n_B"])
    datain, Y, A = cohort["datain"], cohort["Y"], cohort["adjacency"]
    n_sub = datain.shape[3]
    rng = np.random.default_rng(seed)
    pick = rng.choice(n_sub, n_a + n_b, replace=False)
    control, planted = pick[:n_b], pick[n_b:]
    keep = np.sort(np.concatenate([planted, control]))
    remap = {int(i): k for k, i in enumerate(keep)}
    planted_k = np.array([remap[int(i)] for i in planted])
    control_k = np.array([remap[int(i)] for i in control])
    d, y = add_bump(datain[:, :, :, keep], Y[:, :, keep], planted_k, s1, a1)
    d, y = add_bump(d, y, planted_k, s2, a2)
    W = compute_weights(d, A, beta="all", tau=TAU, seeds=(0,))
    groups = groups_from_arms(d.shape[3], planted_k, control_k)
    names = [e for e in ESTIMATORS if e != "yuen" or n_a >= YUEN_MIN_N]
    from tmaps import estimator_boot
    maps = {
        name: estimator_boot(name, y, groups, W, NBOOT, np.random.default_rng(seed + 1 + i))
        for i, name in enumerate(names)
    }
    both = np.concatenate([s1, s2])
    cluster = cluster_maps(maps, A, both, minnbchan=minnbchan, planted_ch1=s1, planted_ch2=s2)
    row = dict(
        seed=int(seed),
        n1=int(rec["n1"]),
        n2=int(rec["n2"]),
        n_A=n_a,
        n_B=n_b,
        minnbchan=int(minnbchan),
        a1=a1,
        a2=a2,
        mean_W_n1_planted=mean_W(W, s1, planted_k),
        mean_W_n2_planted=mean_W(W, s2, planted_k),
        mean_W_n1_unplanted=mean_W(W, s1, control_k),
        mean_W_n2_unplanted=mean_W(W, s2, control_k),
    )
    for est in ESTIMATORS:
        if est in cluster:
            row[f"{est}_n1"] = int(cluster[est]["n1_any"])
            row[f"{est}_n2"] = int(cluster[est]["n2_any"])
        else:
            row[f"{est}_n1"] = ""
            row[f"{est}_n2"] = ""
    return row


def load_mc(path):
    if not path.is_file():
        return [], set()
    rows = []
    have = set()
    with path.open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            r["seed"] = int(r["seed"])
            r["n1"] = int(r["n1"])
            r["n2"] = int(r["n2"])
            r["n_A"] = int(r["n_A"])
            r["n_B"] = int(r["n_B"])
            for k in r:
                if k.startswith("mean_W") or k in ("a1", "a2"):
                    r[k] = float(r[k])
                elif k.endswith("_n1") or k.endswith("_n2"):
                    if k not in ("n1", "n2"):
                        r[k] = int(float(r[k])) if str(r[k]).strip() != "" else ""
            rows.append(r)
            have.add((r["n1"], r["n_A"], r["seed"]))
    return rows, have


def run_cells(label, cells, cohort, minnbchan, csv_path, match_dir):
    match_dir.mkdir(parents=True, exist_ok=True)
    rows, have = load_mc(csv_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    summaries = []

    def dump():
        if not rows:
            return
        fields = list(rows[-1].keys())
        with csv_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            for r in rows:
                writer.writerow(r)

    for spec in cells:
        n1, n2, n_a, n_b = spec
        k = minnbchan(n1) if callable(minnbchan) else int(minnbchan)
        mpath = match_dir / f"match_n1{n1}_nA{n_a}.json"
        try:
            rec = load_or_match(
                mpath, n1=n1, n2=n2, n_a=n_a, n_b=n_b,
                minnbchan=k, cohort=cohort,
            )
        except RuntimeError as e:
            print(f"SKIP {label} n1={n1} n_A={n_a}: {e}", flush=True)
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
            if n_new % 10 == 0:
                dump()
        dump()
        block = [r for r in rows if r["n1"] == n1 and r["n_A"] == n_a]
        w1 = np.mean([r["mean_W_n1_planted"] for r in block])
        w2 = np.mean([r["mean_W_n2_planted"] for r in block])
        print(f"{label} n1={n1} n_A={n_a}  +{n_new}  W_n1={w1:.3f} W_n2={w2:.3f}", flush=True)
        summaries.append(summarize_cell(block, rec, k))
    return summaries


def summarize_cell(block, rec, minnbchan):
    out = dict(
        n1=int(rec["n1"]), n2=int(rec["n2"]),
        n_A=int(rec["n_A"]), n_B=int(rec["n_B"]),
        minnbchan=int(minnbchan),
        n_mc=len(block),
        a1=float(rec["a1"]), a2=float(rec["a2"]),
        n1_names=rec["region1"]["names"],
        n2_names=rec["region2"]["names"],
        region1_channels=rec["region1"]["channels"],
        region2_channels=rec["region2"]["channels"],
    )
    for key in (
        "mean_W_n1_planted", "mean_W_n2_planted",
        "mean_W_n1_unplanted", "mean_W_n2_unplanted",
    ):
        xs = np.array([r[key] for r in block], dtype=float)
        out[key] = float(xs.mean())
        out[key + "_sd"] = float(xs.std(ddof=1)) if xs.size > 1 else 0.0
    for est in ESTIMATORS:
        k1, k2 = f"{est}_n1", f"{est}_n2"
        if k1 in block[0] and str(block[0][k1]).strip() != "":
            out[f"{est}_hit_n1"] = float(np.mean([r[k1] for r in block]))
            out[f"{est}_hit_n2"] = float(np.mean([r[k2] for r in block]))
    return out


def plot_heads(iso_sum, cohort, path):
    coords = np.asarray(cohort["coords"])
    A = cohort["adjacency"]
    x, y = coords[:, 1], coords[:, 0]
    fig, axes = plt.subplots(1, len(iso_sum), figsize=(2.4 * len(iso_sum), 3.4))
    if len(iso_sum) == 1:
        axes = [axes]
    for ax, rec in zip(axes, iso_sum):
        s1 = set(int(i) for i in rec["region1_channels"])
        s2 = set(int(i) for i in rec["region2_channels"])
        n = A.shape[0]
        for i in range(n):
            for j in range(i + 1, n):
                if A[i, j] > 0:
                    ax.plot([x[i], x[j]], [y[i], y[j]], color="0.85", lw=0.5, zorder=1)
        for i in range(n):
            if i in s1:
                fc = "#1f77b4"
            elif i in s2:
                fc = "#d55e00"
            else:
                fc = "white"
            ax.scatter(x[i], y[i], s=28, c=fc, edgecolors="0.3", lw=0.6, zorder=4)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(f"n1={rec['n1']}", fontsize=9)
    fig.suptitle("Isolation neighbourhoods (blue n1, orange n2=8)", fontsize=10, y=0.98)
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {path}")


def plot_margins(iso_sum, bal_sum, path):
    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(10.8, 4.3))
    iso_sum = sorted(iso_sum, key=lambda r: r["n1"])
    bal_sum = sorted(bal_sum, key=lambda r: r["n_A"])
    x_l = np.array([r["n1"] for r in iso_sum], dtype=float)
    x_r = np.array([r["n_A"] for r in bal_sum], dtype=float)
    hi = 1.05
    for key, lab, st in W_KEYS:
        for ax, rows, x in ((ax_l, iso_sum, x_l), (ax_r, bal_sum, x_r)):
            if not rows:
                continue
            y = np.array([r[key] for r in rows], dtype=float)
            sd = np.array([r[key + "_sd"] for r in rows], dtype=float)
            ax.fill_between(x, y - sd, y + sd, color=st["color"], alpha=0.15, lw=0)
            ax.plot(x, y, label=lab, **st, lw=1.8, ms=7)
            hi = max(hi, float(np.nanmax(y + sd)))
    for ax in (ax_l, ax_r):
        ax.axhline(1.0, color="0.6", ls=":", lw=0.9)
        ax.set_ylim(-0.05, hi + 0.08)
        ax.set_ylabel("mean W")
        ax.legend(frameon=False, fontsize=8)
    ax_l.set_xticks(x_l)
    ax_l.set_xlabel("n1 (n2=8, 30 vs 5; minnbchan=0 if n1<3 else 2)")
    ax_l.set_title("Isolation")
    ax_r.set_xticks(x_r)
    ax_r.set_xlabel("n planted (1 vs 8, 30 unplanted, minnbchan=0)")
    ax_r.set_title("Balance")
    fig.suptitle(f"Mass-matched GAE W  ({N_MC} held-out MC, band = mean ± SD)")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {path}")


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    iso_cohort = load_cohort()
    bal_cohort = load_mixed_pad()
    iso_cells = [(n1, ISO_N2, ISO_NA, ISO_NB) for n1 in ISO_N1]
    bal_cells = [(BAL_N1, BAL_N2, n_a, BAL_NB) for n_a in BAL_NA]
    iso_sum = run_cells("iso", iso_cells, iso_cohort, minnb_iso, ISO_CSV, ISO_MATCH)
    bal_sum = run_cells("bal", bal_cells, bal_cohort, BAL_MINNB, BAL_CSV, BAL_MATCH)
    report = dict(
        design=(
            f"Isolation n1={list(ISO_N1)} vs {ISO_N2} on {ISO_NB} vs {ISO_NA}, "
            f"matcher and FWER minnbchan=0 if n1<3 else 2. "
            f"Balance n1={BAL_N1} vs {BAL_N2}, n planted {list(BAL_NA)} vs {BAL_NB}, "
            f"minnbchan={BAL_MINNB}. {N_MC} held-out MCs seeds {EVAL_SEED0}–{EVAL_SEED0 + N_MC - 1}."
        ),
        n_mc=N_MC,
        nboot=NBOOT,
        isolation=iso_sum,
        balance=bal_sum,
        seconds=time.time() - t0,
    )
    (OUT / "iso_bal.json").write_text(json.dumps(report, indent=2) + "\n")
    if iso_sum:
        plot_heads(iso_sum, iso_cohort, OUT / "robust_heads_n1.png")
    plot_margins(iso_sum, bal_sum, OUT / "robust_W_joint.png")
    print(json.dumps({
        "iso": [{k: r[k] for k in ("n1", "a1", "a2", "mean_W_n1_planted", "mean_W_n2_planted") if k in r} for r in iso_sum],
        "bal": [{k: r[k] for k in ("n_A", "a1", "a2", "mean_W_n1_planted", "mean_W_n2_planted") if k in r} for r in bal_sum],
        "seconds": report["seconds"],
    }, indent=2))
    print(f"wrote {OUT / 'iso_bal.json'}")


if __name__ == "__main__":
    main()
