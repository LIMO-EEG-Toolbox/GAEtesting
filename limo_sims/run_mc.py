"""One Monte Carlo: plant → LOCO W → Python t → clusterer from planted n.

n_planted < 3  → paper-style Python glue
n_planted >= 3 → LIMO spatial clustering (STUDY clique)
"""
from __future__ import annotations

import json
import os
import time
from pathlib import Path

import numpy as np

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

HERE = Path(__file__).resolve().parent
import _paths  # noqa: F401

from gae_weighting import compute_weights
from load_cohort import load_cohort

from cluster_limo import cluster_maps as limo_cluster
from cluster_python import cluster_maps as python_cluster
from plant import add_bump, draw_arms, planted_channels
from tmaps import all_estimator_boots, groups_from_arms
from write_limo_maps import write_bundle

TAU = 3.0


def clusterer_for(n_planted):
    n = int(n_planted)
    if n < 1:
        raise ValueError("n_planted must be >= 1")
    return "python" if n < 3 else "limo"


def design_text(n_planted, n_a, n_b, amp, clusterer):
    if clusterer == "python":
        cl = (
            "LIMO cluster-mass FWER in Python (t^2, forming p<0.05, minnbchan=0 so "
            "1- and 2-channel islands stay; bootstrap max-mass threshold)."
        )
    else:
        cl = (
            "LIMO spatial clustering on t^2 (MCC=2, minnbchan=2, forming p<0.05, "
            "bootstrap max-mass FWER)."
        )
    return (
        f"One Monte Carlo: {n_planted}-channel additive bump on {n_a} vs {n_b} maps, "
        f"AMP {amp}. Python two-sample t (mean / Yuen / LOCO-weighted Welch). {cl}"
    )


def run_one(*, n_planted, n_a, n_b, amp, nboot, split_seed, work, out, seed_chan="Pz"):
    clusterer = clusterer_for(n_planted)
    work = Path(work)
    out = Path(out)
    work.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)

    cohort = load_cohort()
    datain, Y, A = cohort["datain"], cohort["Y"], cohort["adjacency"]
    labels = [str(x).strip() for x in cohort["labels"]]
    n_sub = datain.shape[3]
    chans, ch_names = planted_channels(A, labels, n_planted, seed_name=seed_chan)

    rng = np.random.default_rng(split_seed)
    planted, control = draw_arms(n_sub, n_a, n_b, rng)
    keep = np.sort(np.concatenate([planted, control]))
    remap = {int(i): k for k, i in enumerate(keep)}
    planted_k = np.array([remap[int(i)] for i in planted])
    control_k = np.array([remap[int(i)] for i in control])

    d = datain[:, :, :, keep]
    y = Y[:, :, keep]
    d, y = add_bump(d, y, planted_k, chans, amp)

    t0 = time.time()
    W = compute_weights(d, A, beta="all", tau=TAU, seeds=(0,))
    w_sec = time.time() - t0

    groups = groups_from_arms(d.shape[3], planted_k, control_k)
    maps = all_estimator_boots(y, groups, W, nboot, split_seed + 1)

    matlab = None
    if clusterer == "limo":
        write_bundle(
            work / "bundle.mat",
            Y_planted=y[:, :, planted_k],
            Y_control=y[:, :, control_k],
            W=W,
            adjacency=A,
            labels=labels,
            planted_ch=chans,
            planted_idx=planted_k,
            control_idx=control_k,
            amp=amp,
            nboot=nboot,
            maps=maps,
        )
        matlab, cluster = limo_cluster(work)
    else:
        cluster = python_cluster(maps, A, chans)

    w_p = W[np.ix_(chans, planted_k)]
    w_c = W[np.ix_(chans, control_k)]
    report = {
        "design": design_text(n_planted, n_a, n_b, amp, clusterer),
        "clusterer": clusterer,
        "channels": ch_names,
        "channel_index": chans.tolist(),
        "n_A": int(n_a),
        "n_B": int(n_b),
        "amp": float(amp),
        "nboot": int(nboot),
        "mean_W_planted": float(w_p.mean()),
        "mean_W_control": float(w_c.mean()),
        "mean_W_planted_other": float(np.delete(W[:, planted_k], chans, axis=0).mean()),
        "gae_seconds": w_sec,
        "matlab": matlab,
        "cluster": cluster,
    }
    path = out / "check.json"
    path.write_text(json.dumps(report, indent=2) + "\n")
    return report, path
