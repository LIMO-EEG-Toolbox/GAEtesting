#!/usr/bin/env python3
"""Forming hit-rate, FWER hit-rate, and bootstrap max-mass H0 (2 vs 8, 30 vs 6).

    python limo_sims/hitrate.py
"""
from __future__ import annotations

import csv
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

from gae_weighting import compute_weights
from load_cohort import load_cohort

from limo_mass import cluster_test
from match_limo import match_mass_amps
from plant import PEAK_MS, add_bump, times_ms
from tmaps import ESTIMATORS, estimator_boot, groups_from_arms

N1 = 2
N2 = 8
N_A = 6
N_B = 30
MINNBCHAN = 0
N_MC = 200
EVAL_SEED0 = 10
NBOOT = 101
TAU = 3.0
OUT = HERE / "results" / "hitrate_2v8_a6"
MATCH_PATH = OUT / "match.json"
CSV_PATH = OUT / "mc.csv"
CELLS = OUT / "cells"

EST_LABEL = {"mean": "mean", "yuen": "Yuen 20%", "weighted": "weighted"}
EST_COLOR = {"mean": "#1f77b4", "yuen": "#d55e00", "weighted": "#009e73"}


def mean_W(W, chans, who):
    sl = W[np.asarray(chans, dtype=int)][:, np.asarray(who, dtype=int)]
    return float(np.mean(sl)) if sl.size else float("nan")


def load_or_match(cohort):
    OUT.mkdir(parents=True, exist_ok=True)
    if MATCH_PATH.is_file():
        rec = json.loads(MATCH_PATH.read_text())
        if int(rec.get("minnbchan", -1)) == MINNBCHAN and int(rec.get("n1", -1)) == N1:
            print(f"match cache {MATCH_PATH.name}  a1={rec['a1']:.2f} a2={rec['a2']:.2f}", flush=True)
            return rec
    print(f"matching n1={N1} n2={N2} {N_B} vs {N_A} minnbchan={MINNBCHAN}", flush=True)
    rec = match_mass_amps(n1=N1, n2=N2, n_a=N_A, n_b=N_B, minnbchan=MINNBCHAN, cohort=cohort)
    OUT.mkdir(parents=True, exist_ok=True)
    MATCH_PATH.write_text(json.dumps(rec, indent=2) + "\n")
    print(f"  a1={rec['a1']:.2f} a2={rec['a2']:.2f}", flush=True)
    return rec


def cell_path(mc):
    return CELLS / f"mc_{mc:03d}.npz"


def one_mc(cohort, rec, seed, mc):
    s1 = np.asarray(rec["region1"]["channels"], dtype=int)
    s2 = np.asarray(rec["region2"]["channels"], dtype=int)
    a1, a2 = float(rec["a1"]), float(rec["a2"])
    datain, Y, A = cohort["datain"], cohort["Y"], cohort["adjacency"]
    n_sub = datain.shape[3]
    rng = np.random.default_rng(seed)
    pick = rng.choice(n_sub, N_A + N_B, replace=False)
    control, planted = pick[:N_B], pick[N_B:]
    keep = np.sort(np.concatenate([planted, control]))
    remap = {int(i): k for k, i in enumerate(keep)}
    planted_k = np.array([remap[int(i)] for i in planted])
    control_k = np.array([remap[int(i)] for i in control])
    d, y = add_bump(datain[:, :, :, keep], Y[:, :, keep], planted_k, s1, a1)
    d, y = add_bump(d, y, planted_k, s2, a2)
    W = compute_weights(d, A, beta="all", tau=TAU, seeds=(0,))
    groups = groups_from_arms(d.shape[3], planted_k, control_k)
    maps = {
        name: estimator_boot(name, y, groups, W, NBOOT, np.random.default_rng(seed + 1 + i))
        for i, name in enumerate(ESTIMATORS)
    }
    both = np.concatenate([s1, s2])
    payload = dict(
        mc=np.int32(mc), seed=np.int32(seed),
        coords=np.asarray(cohort["coords"], dtype=np.float64),
        names=np.array([str(x).strip() for x in cohort["labels"]]),
        s1=s1, s2=s2,
        W_n1=mean_W(W, s1, planted_k),
        W_n2=mean_W(W, s2, planted_k),
    )
    row = dict(
        mc=int(mc), seed=int(seed),
        mean_W_n1_planted=payload["W_n1"],
        mean_W_n2_planted=payload["W_n2"],
    )
    for name, m in maps.items():
        ct = cluster_test(
            m["t"], m["p"], m["bootT"], m["bootP"], A,
            planted=both, planted1=s1, planted2=s2,
            minnbchan=MINNBCHAN, return_mask=True,
        )
        payload[f"form_n1_{name}"] = ct["form_n1"]
        payload[f"form_n2_{name}"] = ct["form_n2"]
        payload[f"form_chance_{name}"] = ct["form_chance"]
        payload[f"fwer_n1_{name}"] = ct["fwer_n1"]
        payload[f"fwer_n2_{name}"] = ct["fwer_n2"]
        payload[f"fwer_chance_{name}"] = ct["fwer_chance"]
        payload[f"boot_max_{name}"] = ct["boot_max"]
        row[f"max_th_{name}"] = float(ct["max_th"])
        row[f"obs_max_{name}"] = float(ct["obs_max_mass"])
        row[f"n1_form_{name}"] = int(ct["form_n1"].any())
        row[f"n2_form_{name}"] = int(ct["form_n2"].any())
        row[f"n1_fwer_{name}"] = int(ct["n1_any"])
        row[f"n2_fwer_{name}"] = int(ct["n2_any"])
    CELLS.mkdir(parents=True, exist_ok=True)
    tmp = CELLS / f"mc_{mc:03d}.tmp.npz"
    np.savez_compressed(tmp, **payload)
    tmp.replace(cell_path(mc))
    return row


def load_rows():
    if not CSV_PATH.is_file():
        return [], set()
    rows = []
    have = set()
    with CSV_PATH.open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            r["mc"] = int(r["mc"])
            r["seed"] = int(r["seed"])
            rows.append(r)
            have.add(r["mc"])
    return rows, have


def stack():
    files = sorted(p for p in CELLS.glob("mc_*.npz") if not p.name.endswith(".tmp.npz"))
    if not files:
        return None
    first = np.load(files[0], allow_pickle=True)
    acc = {}
    keys = []
    for est in ESTIMATORS:
        for kind in ("form_n1", "form_n2", "form_chance", "fwer_n1", "fwer_n2", "fwer_chance"):
            keys.append(f"{kind}_{est}")
    for k in keys:
        acc[k] = np.zeros(first[k].shape, dtype=np.float64)
    boot = {est: [] for est in ESTIMATORS}
    th = {est: [] for est in ESTIMATORS}
    for p in files:
        z = np.load(p, allow_pickle=True)
        for k in keys:
            acc[k] += z[k]
        for est in ESTIMATORS:
            boot[est].append(np.asarray(z[f"boot_max_{est}"], dtype=np.float64))
            th[est].append(float(np.percentile(z[f"boot_max_{est}"], 95)))
    n = len(files)
    for k in acc:
        acc[k] /= n
    return dict(
        n=n,
        names=[str(x) for x in first["names"].tolist()],
        coords=first["coords"],
        s1=first["s1"], s2=first["s2"],
        boot={est: np.concatenate(v) for est, v in boot.items()},
        max_th={est: np.array(v) for est, v in th.items()},
        **acc,
    )


def _panel_grid(st, kinds, title, path, nrows=2):
    order = np.argsort(-st["coords"][:, 0])
    t_ms = times_ms(st["form_n1_mean"].shape[1])
    ext = [float(t_ms[0]), float(t_ms[-1]), st["form_n1_mean"].shape[0] - 0.5, -0.5]
    fig, axes = plt.subplots(nrows, 3, figsize=(13.0, 3.6 * nrows), sharex=True, sharey=True)
    axes = np.atleast_2d(axes)
    im = None
    for j, est in enumerate(ESTIMATORS):
        for i, kind in enumerate(kinds):
            ax = axes[i, j]
            vis = st[f"{kind}_{est}"][order]
            im = ax.imshow(
                vis, aspect="auto", cmap="YlOrRd", vmin=0, vmax=1,
                origin="upper", extent=ext, interpolation="nearest",
            )
            ax.axvline(PEAK_MS, color="0.35", ls=":", lw=0.8)
            if i == 0:
                ax.set_title(EST_LABEL[est])
            if j == 0:
                ax.set_ylabel("channel (ant → post)" if i == 0 else "")
            if i == nrows - 1:
                ax.set_xlabel("time (ms)")
            labs = [st["names"][k] for k in order]
            ax.set_yticks(np.arange(0, len(order), 2))
            ax.set_yticklabels(labs[::2], fontsize=5)
    fig.subplots_adjust(right=0.90, wspace=0.12, hspace=0.18)
    cax = fig.add_axes([0.92, 0.18, 0.018, 0.64])
    fig.colorbar(im, cax=cax, label="hit rate")
    fig.suptitle(f"{title}  ·  {st['n']} MC  ·  2 vs 8, 30 vs 6", y=0.98)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_all(st):
    _panel_grid(st, ("form_n1", "form_n2"), "Forming overlap hit-rate", OUT / "hit_forming.png")
    _panel_grid(st, ("fwer_n1", "fwer_n2"), "FWER overlap hit-rate", OUT / "hit_fwer.png")
    _panel_grid(st, ("form_chance",), "Forming chance-blob hit-rate", OUT / "hit_form_chance.png", nrows=1)
    _panel_grid(st, ("fwer_chance",), "FWER chance-blob hit-rate", OUT / "hit_fwer_chance.png", nrows=1)

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    x = np.linspace(0, 1, 400)
    for est in ESTIMATORS:
        v = np.sort(st["boot"][est])
        y = np.arange(1, v.size + 1) / v.size
        ax.plot(v, y, color=EST_COLOR[est], lw=1.8, label=EST_LABEL[est])
    ax.set_xlabel("bootstrap max cluster mass (t²)")
    ax.set_ylabel("ECDF")
    ax.set_title("H0 max-mass  ·  pooled boots")
    ax.legend(frameon=False)
    fig.savefig(OUT / "boot_max_ecdf.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    data = [st["max_th"][est] for est in ESTIMATORS]
    bp = ax.boxplot(data, patch_artist=True)
    ax.set_xticklabels([EST_LABEL[e] for e in ESTIMATORS])
    for patch, est in zip(bp["boxes"], ESTIMATORS):
        patch.set_facecolor(EST_COLOR[est])
        patch.set_alpha(0.45)
    ax.set_ylabel("max_th (95% of boot max mass)")
    ax.set_title("FWER bar across Monte Carlos")
    fig.savefig(OUT / "max_th_box.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def main():
    t0 = time.time()
    cohort = load_cohort()
    rec = load_or_match(cohort)
    OUT.mkdir(parents=True, exist_ok=True)
    rows, have = load_rows()
    print(f"resume {len(have)}/{N_MC} MCs", flush=True)
    fields = None
    with CSV_PATH.open("a", newline="", encoding="utf-8") as f:
        writer = None
        for mc in range(N_MC):
            if mc in have:
                continue
            seed = EVAL_SEED0 + mc
            row = one_mc(cohort, rec, seed, mc)
            if writer is None:
                fields = list(row.keys())
                writer = csv.DictWriter(f, fieldnames=fields)
                if CSV_PATH.stat().st_size == 0:
                    writer.writeheader()
            writer.writerow(row)
            rows.append(row)
            have.add(mc)
            f.flush()
            if (mc + 1) % 10 == 0:
                print(f"  mc {mc + 1}/{N_MC}", flush=True)
    st = stack()
    if st is not None:
        plot_all(st)
    report = dict(
        design=(
            "2 vs 8 plants on 6 vs 30, matcher 10/10 LIMO FWER then equal t^2 mass, "
            "minnbchan=0. Each MC: GAE W, Python mean/Yuen/weighted t, limo_mass "
            "cluster_test. Forming occupancy, FWER occupancy, bootstrap max-mass H0."
        ),
        n_mc=N_MC,
        nboot=NBOOT,
        minnbchan=MINNBCHAN,
        a1=rec["a1"],
        a2=rec["a2"],
        n_done=len(have),
        seconds=time.time() - t0,
    )
    (OUT / "hitrate.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ("a1", "a2", "n_done", "seconds")}, indent=2))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
