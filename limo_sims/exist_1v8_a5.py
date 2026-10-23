#!/usr/bin/env python3
"""Paper existence proof: both plants on n=5, GAE, weighted FWER map.

    python limo_sims/exist_1v8_a5.py
"""
from __future__ import annotations

import json
import os
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

from match_limo import SPLIT_SEED, _sig_panel, frozen_split
from plant import add_bump
from tmaps import all_estimator_boots

MATCH = HERE / "results" / "match_1v8_a5" / "match.json"
OUT = HERE / "results" / "exist_1v8_a5"
TAU = 3.0
NBOOT = 101
ALPHA = 0.05


def mean_W(W, chans, who):
    sl = W[np.asarray(chans, dtype=int)][:, np.asarray(who, dtype=int)]
    return float(np.mean(sl)) if sl.size else float("nan")


def main():
    rec = json.loads(MATCH.read_text())
    cohort = load_cohort()
    ctx = frozen_split(
        cohort, rec["n_A"], rec["n_B"], rec["split_seed"], rec["minnbchan"],
    )
    s1 = np.asarray(rec["region1"]["channels"], dtype=int)
    s2 = np.asarray(rec["region2"]["channels"], dtype=int)
    a1, a2 = float(rec["a1"]), float(rec["a2"])
    d, y = ctx["datain"], ctx["Y"]
    d, y = add_bump(d, y, ctx["who_a"], s1, a1)
    d, y = add_bump(d, y, ctx["who_a"], s2, a2)
    W = compute_weights(d, ctx["A"], beta="all", tau=TAU, seeds=(0,))
    groups = ctx["groups"]
    maps = all_estimator_boots(y, groups, W, NBOOT, SPLIT_SEED + 11)
    both = np.concatenate([s1, s2])
    names = ctx["names"]
    order = np.argsort(-np.asarray(ctx["coords"])[:, 0])
    Wp = W[:, ctx["who_a"]]
    fig, ax = plt.subplots(figsize=(6.4, 5.6))
    im = ax.imshow(Wp[order], aspect="auto", cmap="viridis", vmin=0.0, vmax=1.3,
                   interpolation="nearest")
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels([names[i] for i in order], fontsize=6)
    planted = set(int(i) for i in both)
    for tick, ch in zip(ax.get_yticklabels(), order):
        if int(ch) in planted:
            tick.set_fontweight("bold")
    ax.set_xlabel("planted maps")
    ax.set_title(
        f"W on n={rec['n_A']}  Oz={mean_W(W, s1, ctx['who_a']):.2f}  "
        f"n2={mean_W(W, s2, ctx['who_a']):.2f}"
    )
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02, label="W")
    OUT.mkdir(parents=True, exist_ok=True)
    wpath = OUT / "exist_W.png"
    fig.savefig(wpath, dpi=140, bbox_inches="tight")
    plt.close(fig)

    def fwer_mask(name):
        from limo_mass import cluster_test
        m = maps[name]
        ct = cluster_test(
            m["t"], m["p"], m["bootT"], m["bootP"], ctx["A"],
            planted=both, planted1=s1, planted2=s2,
            minnbchan=int(rec["minnbchan"]), return_mask=True,
        )
        return m["t"], m["p"], ct

    t_m, p_m, c_m = fwer_mask("mean")
    t_w, p_w, c_w = fwer_mask("weighted")
    vmax = max(
        float(np.nanmax(np.abs(t_m[p_m <= ALPHA]))) if np.any(p_m <= ALPHA) else 1.0,
        float(np.nanmax(np.abs(t_w[p_w <= ALPHA]))) if np.any(p_w <= ALPHA) else 1.0,
        1.0,
    )
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 5.8), sharey=True)
    _sig_panel(
        axes[0], t_m, p_m, c_m["mask"], ctx["coords"], names, both,
        f"mean  n1={int(c_m['n1_any'])} n2={int(c_m['n2_any'])}\n"
        f"mass n1={c_m['fwer_overlap_mass']:.0f}",
        vmax,
    )
    _sig_panel(
        axes[1], t_w, p_w, c_w["mask"], ctx["coords"], names, both,
        f"weighted  n1={int(c_w['n1_any'])} n2={int(c_w['n2_any'])}\n"
        f"mass={c_w['fwer_overlap_mass']:.0f}",
        vmax,
    )
    axes[1].set_ylabel("")
    fig.subplots_adjust(left=0.08, right=0.90, top=0.86, bottom=0.12, wspace=0.12)
    cax = fig.add_axes([0.92, 0.12, 0.018, 0.72])
    fig.colorbar(axes[0].images[0], cax=cax, label="t (forming p<0.05)")
    tpath = OUT / "exist_tmap_weighted.png"
    fig.savefig(tpath, dpi=140)
    plt.close(fig)

    out = dict(
        design=(
            f"Existence on matcher split seed {rec['split_seed']}: both plants on n="
            f"{rec['n_A']} vs {rec['n_B']}, a1={a1:.2f} a2={a2:.2f}. GAE then "
            f"mean/weighted LIMO FWER minnbchan={rec['minnbchan']}."
        ),
        a1=a1, a2=a2,
        mean_W_n1=mean_W(W, s1, ctx["who_a"]),
        mean_W_n2=mean_W(W, s2, ctx["who_a"]),
        mean_n1=bool(c_m["n1_any"]),
        mean_n2=bool(c_m["n2_any"]),
        weighted_n1=bool(c_w["n1_any"]),
        weighted_n2=bool(c_w["n2_any"]),
        cluster={k: {kk: vv for kk, vv in c.items() if kk != "mask"}
                 for k, c in dict(mean=c_m, weighted=c_w).items()},
    )
    (OUT / "exist.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({k: out[k] for k in (
        "mean_W_n1", "mean_W_n2", "mean_n1", "mean_n2", "weighted_n1", "weighted_n2",
    )}, indent=2))
    print(f"wrote {wpath}")
    print(f"wrote {tpath}")


if __name__ == "__main__":
    main()
