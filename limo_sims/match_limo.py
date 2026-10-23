#!/usr/bin/env python3
"""Final matcher: 10-draw mean LIMO FWER floor, then equal t^2 mass.

    python limo_sims/match_limo.py          # minnbchan=2, 3 vs 8
    python limo_sims/match_1v8_a5.py        # minnbchan=0, 1 vs 8

Mean Welch t only. No GAE. Each neighbourhood planted alone.
Floor: smallest amp with 10/10 mean LIMO cluster-mass FWER hits on calibration
splits (plants alone). Then equalize median FWER-cluster t^2 mass. minnbchan
must match the later test. Confirmation MCs use new seeds.
"""
from __future__ import annotations

import json
import os
from collections import deque
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

HERE = Path(__file__).resolve().parent
import _paths  # noqa: F401

from gae_weighting import weighted_welch_t
from load_cohort import load_cohort

from limo_mass import cluster_test
from plant import add_bump, draw_arms, times_ms, PEAK_MS
from tmaps import estimator_boot, groups_from_arms, t_and_p

N1 = 3
N2 = 8
MINNBCHAN = 2
N_A = 6
N_B = 30
SPLIT_SEED = 0
NBOOT = 101
N_MATCH_MC = 10
N_HIT = 10
MATCH_SEEDS = tuple(range(N_MATCH_MC))
ALPHA = 0.05
A_LO = 0.0
A_HI = 256.0
MASS_REL = 0.02
N_SEARCH = 16
OUT = HERE / "results" / "match_limo"


def adj_list(A):
    return [np.flatnonzero(A[i] > 0).astype(int) for i in range(A.shape[0])]


def peak_index(n_time):
    t = times_ms(n_time)
    return int(np.argmin(np.abs(t - PEAK_MS)))


def share_edge(s1, s2, A):
    s1 = np.asarray(s1, dtype=int)
    s2 = np.asarray(s2, dtype=int)
    return bool(np.any(A[np.ix_(s1, s2)] > 0))


def bfs_grow(seed, n, neigh, blocked):
    vis = set(int(x) for x in blocked)
    if int(seed) in vis:
        return None
    q = deque([int(seed)])
    vis.add(int(seed))
    out = []
    while q and len(out) < n:
        u = q.popleft()
        out.append(u)
        for v in neigh[u]:
            v = int(v)
            if v not in vis:
                vis.add(v)
                q.append(v)
    if len(out) < n:
        return None
    return np.array(out[:n], dtype=int)


def _try_axis(coord, n1, n2, A, neigh):
    seed1 = int(np.argmin(coord))
    seed2 = int(np.argmax(coord))
    for order in ("s1", "s2"):
        if order == "s1":
            s1 = bfs_grow(seed1, n1, neigh, blocked=())
            s2 = None if s1 is None else bfs_grow(seed2, n2, neigh, blocked=s1)
        else:
            s2 = bfs_grow(seed2, n2, neigh, blocked=())
            s1 = None if s2 is None else bfs_grow(seed1, n1, neigh, blocked=s2)
        if s1 is None or s2 is None:
            continue
        if len(set(s1) & set(s2)):
            continue
        if share_edge(s1, s2, A):
            continue
        return s1, s2
    return None


def grow_neighbourhoods(n1, n2, A, coords, names):
    n1, n2 = int(n1), int(n2)
    if n1 < 1 or n2 < 1:
        raise ValueError("need n1,n2 >= 1")
    neigh = adj_list(A)
    ap = _try_axis(coords[:, 0], n1, n2, A, neigh)
    if ap is not None:
        s1, s2 = ap
        axis = "anterior-posterior"
    else:
        lr = _try_axis(coords[:, 1], n1, n2, A, neigh)
        if lr is None:
            raise ValueError(f"could not grow disconnected neighbourhoods of size {n1} and {n2}")
        s1, s2 = lr
        axis = "left-right"
    print(f"neighbourhoods on {axis}")
    print(f"  n1={n1}  {[names[i] for i in s1]}")
    print(f"  n2={n2}  {[names[i] for i in s2]}")
    return s1, s2, axis


def mean_tp(Y, groups):
    n_chan, n_sub = Y.shape[1], Y.shape[2]
    ones = np.ones((n_chan, n_sub), dtype=float)
    tmap, df, _ = weighted_welch_t(Y, groups, ones)
    t, p, d = t_and_p(tmap, df)
    return t, p, d, ones


def frozen_split(cohort, n_a, n_b, seed, minnbchan=MINNBCHAN, verbose=True):
    Y = cohort["Y"]
    datain = cohort["datain"]
    A = cohort["adjacency"]
    names = [str(x).strip() for x in cohort["labels"]]
    coords = np.asarray(cohort["coords"], dtype=np.float64)
    n_sub = Y.shape[2]
    rng = np.random.default_rng(seed)
    planted, control = draw_arms(n_sub, n_a, n_b, rng)
    keep = np.sort(np.concatenate([planted, control]))
    remap = {int(i): k for k, i in enumerate(keep)}
    planted_k = np.array([remap[int(i)] for i in planted])
    control_k = np.array([remap[int(i)] for i in control])
    y = Y[:, :, keep]
    d = datain[:, :, :, keep]
    groups = groups_from_arms(y.shape[2], planted_k, control_k)
    _, _, _, ones = mean_tp(y, groups)
    if verbose:
        print(f"frozen {n_b} vs {n_a}  seed={seed}  minnbchan={minnbchan}  FWER floor nboot={NBOOT}")
    return dict(
        Y=y, datain=d, groups=groups, who_a=planted_k, who_b=control_k,
        A=A, coords=coords, names=names, ones=ones, seed=int(seed),
        peak_i=peak_index(y.shape[0]), minnbchan=int(minnbchan),
    )


def eval_plant(ctx, chans, amp):
    y = ctx["Y"]
    d = ctx["datain"]
    _, yp = add_bump(d, y, ctx["who_a"], chans, amp)
    maps = estimator_boot(
        "mean", yp, ctx["groups"], ctx["ones"], NBOOT,
        np.random.default_rng(int(ctx.get("seed", SPLIT_SEED)) + 1),
    )
    ct = cluster_test(
        maps["t"], maps["p"], maps["bootT"], maps["bootP"], ctx["A"],
        planted=chans, minnbchan=ctx["minnbchan"],
    )
    sl = maps["t"][np.asarray(chans, dtype=int)]
    tmax = float(np.nanmax(np.abs(sl))) if sl.size else float("nan")
    tpeak = float(np.nanmax(np.abs(sl[:, ctx["peak_i"]]))) if sl.size else float("nan")
    return dict(
        Yp=yp, tmap=maps["t"], pmap=maps["p"],
        mass=float(ct["fwer_overlap_mass"]),
        detected=bool(ct["planted_any_cluster"]),
        max_th=float(ct["max_th"]),
        min_cluster_p=ct["min_cluster_p"],
        max_abs_t=tmax, t_peak=tpeak,
        forming_overlap_mass=float(ct["forming_overlap_mass"]),
    )


def _bisect(pred, lo, hi, n=N_SEARCH):
    if pred(lo):
        return lo
    if not pred(hi):
        return hi
    a, b = lo, hi
    for _ in range(n):
        m = 0.5 * (a + b)
        if pred(m):
            b = m
        else:
            a = m
    return b


def detect_floor(ctxs, chans):
    cache = {}
    names = ctxs[0]["names"]

    def mass_of(a):
        a = float(a)
        if a not in cache:
            recs = [eval_plant(ctx, chans, a) for ctx in ctxs]
            hits = [bool(r["detected"]) for r in recs]
            masses = [float(r["mass"]) for r in recs]
            n_hit = int(sum(hits))
            cache[a] = dict(
                detected=n_hit >= N_HIT,
                n_hit=n_hit,
                n_mc=len(ctxs),
                hit_rate=n_hit / len(ctxs),
                mass=float(np.median(masses)),
                recs=recs,
                max_th=float(np.median([r["max_th"] for r in recs])),
                min_cluster_p=recs[0]["min_cluster_p"],
                max_abs_t=float(np.median([r["max_abs_t"] for r in recs])),
                t_peak=float(np.median([r["t_peak"] for r in recs])),
                forming_overlap_mass=float(np.median([r["forming_overlap_mass"] for r in recs])),
            )
            print(
                f"  a={a:.2f}  FWER {n_hit}/{len(ctxs)}  "
                f"median mass={cache[a]['mass']:.1f}  "
                f"{[names[i] for i in chans]}",
                flush=True,
            )
        return cache[a]

    mass_of.cache = cache

    def detected(a):
        return mass_of(a)["detected"]

    if not detected(A_HI):
        rec = mass_of(A_HI)
        raise RuntimeError(
            f"amp {A_HI} FWER {rec['n_hit']}/{rec['n_mc']} on "
            f"{[names[i] for i in chans]} (need {N_HIT})"
        )
    a0 = _bisect(detected, A_LO, A_HI)
    return float(a0), mass_of(a0), mass_of


def match_to_mass(mass_of, a_lo, target):
    tgt = float(target)

    def ok(a):
        return mass_of(a)["mass"] >= tgt * (1.0 - MASS_REL)

    if ok(a_lo):
        return float(a_lo), mass_of(a_lo)
    if not ok(A_HI):
        return float(A_HI), mass_of(A_HI)
    a = _bisect(ok, a_lo, A_HI)
    return float(a), mass_of(a)


def region_out(chans, names, side, a, rec):
    return dict(
        n=int(len(chans)),
        channels=[int(i) for i in chans],
        names=[names[i] for i in chans],
        side=side,
        amp=float(a),
        cluster_mass=float(rec["mass"]),
        forming_overlap_mass=float(rec["forming_overlap_mass"]),
        max_th=float(rec["max_th"]),
        min_cluster_p=rec["min_cluster_p"],
        fwer_hit=bool(rec["detected"]),
        fwer_n_hit=int(rec.get("n_hit", int(rec["detected"]))),
        fwer_n_mc=int(rec.get("n_mc", 1)),
        max_abs_t=float(rec["max_abs_t"]),
        t_peak=float(rec["t_peak"]),
    )


def match_mass_amps(n1=N1, n2=N2, n_a=N_A, n_b=N_B, minnbchan=MINNBCHAN, cohort=None):
    cohort = cohort if cohort is not None else load_cohort()
    ctxs = [
        frozen_split(cohort, n_a, n_b, seed, minnbchan=minnbchan, verbose=(seed == MATCH_SEEDS[0]))
        for seed in MATCH_SEEDS
    ]
    ctx = ctxs[0]
    s1, s2, axis = grow_neighbourhoods(n1, n2, ctx["A"], ctx["coords"], ctx["names"])
    print(f"floor n1  ({N_HIT}/{N_MATCH_MC} mean FWER)", flush=True)
    a1_floor, rec1_floor, mass1 = detect_floor(ctxs, s1)
    print(f"floor n2  ({N_HIT}/{N_MATCH_MC} mean FWER)", flush=True)
    a2_floor, rec2_floor, mass2 = detect_floor(ctxs, s2)
    M = max(rec1_floor["mass"], rec2_floor["mass"])
    a1, rec1 = match_to_mass(mass1, a1_floor, M)
    a2, rec2 = match_to_mass(mass2, a2_floor, M)
    side1, side2 = ("posterior", "anterior") if axis == "anterior-posterior" else ("right", "left")
    return dict(
        design=(
            f"{n_b} vs {n_a} mean Welch t. Boxcar Gaussian (peak 380 ms, width 40 ms) "
            f"on n1={n1} vs n2={n2} disconnected {axis} neighbourhoods. "
            f"Floor is mean LIMO FWER (minnbchan={minnbchan}, nboot={NBOOT}) on "
            f"{N_MATCH_MC} calibration splits (seeds {list(MATCH_SEEDS)}), plants alone: "
            f"smallest amp with {N_HIT}/{N_MATCH_MC} hits. Then equalize median "
            f"FWER-cluster t^2 mass. Same minnbchan as the later test. "
            f"Confirmation MCs must use new seeds."
        ),
        clusterer="limo" if int(minnbchan) >= 2 else "limo_thin",
        minnbchan=int(minnbchan),
        nboot=NBOOT,
        match_seeds=list(MATCH_SEEDS),
        n_hit=N_HIT,
        n1=int(n1),
        n2=int(n2),
        n_A=int(n_a),
        n_B=int(n_b),
        split_seed=SPLIT_SEED,
        axis=axis,
        forming_alpha=ALPHA,
        matched_mass=float(M),
        a1=float(a1),
        a2=float(a2),
        region1=region_out(s1, ctx["names"], side1, a1, rec1),
        region2=region_out(s2, ctx["names"], side2, a2, rec2),
        floors=dict(
            a1=float(a1_floor), mass1=float(rec1_floor["mass"]),
            hit1=rec1_floor["n_hit"],
            a2=float(a2_floor), mass2=float(rec2_floor["mass"]),
            hit2=rec2_floor["n_hit"],
        ),
    )


def _sig_panel(ax, tmap, pmap, fwer_mask, coords, names, planted, title, vmax):
    t_ms = times_ms(tmap.shape[1])
    order = np.argsort(-np.asarray(coords)[:, 0])
    sig = np.asarray(pmap) <= ALPHA
    vis = np.asarray(tmap, dtype=np.float64)[order]
    vis[~sig[order]] = np.nan
    ext = [float(t_ms[0]), float(t_ms[-1]), tmap.shape[0] - 0.5, -0.5]
    im = ax.imshow(
        vis, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
        origin="upper", extent=ext, interpolation="nearest",
    )
    if fwer_mask is not None and np.any(fwer_mask):
        on = np.zeros(tmap.shape, dtype=float)
        on[np.asarray(fwer_mask) > 0] = 1.0
        ax.contour(
            on[order], levels=[0.5], colors=["#111111"], linewidths=1.1,
            origin="upper", extent=ext,
        )
    planted = set(int(i) for i in planted)
    yticks = np.arange(len(order))
    labs = [names[i] for i in order]
    ax.set_yticks(yticks)
    ax.set_yticklabels(labs, fontsize=6)
    for tick, ch in zip(ax.get_yticklabels(), order):
        if int(ch) in planted:
            tick.set_color("#d55e00")
            tick.set_fontweight("bold")
    ax.axvline(PEAK_MS, color="0.35", ls=":", lw=0.9)
    ax.set_xlabel("time (ms)")
    ax.set_ylabel("channel (ant → post)")
    ax.set_title(title, fontsize=9)
    return im


def plot_sig_maps(rec, out_dir):
    """Mean t, forming p<0.05 colour; black contour = FWER mask. Matcher split."""
    from load_cohort import load_cohort

    cohort = load_cohort()
    ctx = frozen_split(
        cohort, rec["n_A"], rec["n_B"], rec["split_seed"], rec["minnbchan"],
    )
    s1 = np.asarray(rec["region1"]["channels"], dtype=int)
    s2 = np.asarray(rec["region2"]["channels"], dtype=int)
    a1, a2 = float(rec["a1"]), float(rec["a2"])
    k = int(rec["minnbchan"])

    def one(chans, amp, extra=None):
        d, y = ctx["datain"], ctx["Y"]
        _, yp = add_bump(d, y, ctx["who_a"], chans, amp)
        if extra is not None:
            _, yp = add_bump(d, yp, ctx["who_a"], extra[0], extra[1])
        maps = estimator_boot(
            "mean", yp, ctx["groups"], ctx["ones"], NBOOT,
            np.random.default_rng(SPLIT_SEED + 1),
        )
        planted = chans if extra is None else np.concatenate([chans, extra[0]])
        p1, p2 = (None, None) if extra is None else (s1, s2)
        ct = cluster_test(
            maps["t"], maps["p"], maps["bootT"], maps["bootP"], ctx["A"],
            planted=planted, planted1=p1, planted2=p2, minnbchan=k,
            return_mask=True,
        )
        return maps["t"], maps["p"], ct

    t1, p1, c1 = one(s1, a1)
    t2, p2, c2 = one(s2, a2)
    tb, pb, cb = one(s1, a1, extra=(s2, a2))
    vmax = max(
        float(np.nanmax(np.abs(t1[p1 <= ALPHA]))) if np.any(p1 <= ALPHA) else 1.0,
        float(np.nanmax(np.abs(t2[p2 <= ALPHA]))) if np.any(p2 <= ALPHA) else 1.0,
        float(np.nanmax(np.abs(tb[pb <= ALPHA]))) if np.any(pb <= ALPHA) else 1.0,
        1.0,
    )
    fig, axes = plt.subplots(
        1, 3, figsize=(14.4, 5.8), sharey=True,
        gridspec_kw={"width_ratios": [1, 1, 1], "wspace": 0.12},
    )
    titles = [
        f"n1 alone  a={a1:.1f}\nFWER hit={int(c1['planted_any_cluster'])}  "
        f"mass={c1['fwer_overlap_mass']:.0f}  T={c1['max_th']:.0f}",
        f"n2 alone  a={a2:.1f}\nFWER hit={int(c2['planted_any_cluster'])}  "
        f"mass={c2['fwer_overlap_mass']:.0f}  T={c2['max_th']:.0f}",
        f"both (display)\nFWER n1={int(cb['n1_any'])} n2={int(cb['n2_any'])}  "
        f"minnbchan={k}",
    ]
    im = None
    for ax, tmap, pmap, ct, title, planted in (
        (axes[0], t1, p1, c1, titles[0], s1),
        (axes[1], t2, p2, c2, titles[1], s2),
        (axes[2], tb, pb, cb, titles[2], np.concatenate([s1, s2])),
    ):
        im = _sig_panel(
            ax, tmap, pmap, ct["mask"], ctx["coords"], ctx["names"],
            planted, title, vmax,
        )
    axes[1].set_ylabel("")
    axes[2].set_ylabel("")
    fig.subplots_adjust(left=0.07, right=0.90, top=0.86, bottom=0.12, wspace=0.14)
    cax = fig.add_axes([0.92, 0.12, 0.018, 0.72])
    fig.colorbar(im, cax=cax, label="mean t (forming p<0.05)")
    fig.suptitle(
        "Matcher split · colour = forming · black = FWER cluster",
        fontsize=11, y=0.98,
    )
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "sig_map.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    print(f"wrote {path}")
    return path


def main():
    rec = match_mass_amps()
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "match.json"
    path.write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps({
        "axis": rec["axis"],
        "a1": rec["a1"],
        "a2": rec["a2"],
        "matched_mass": rec["matched_mass"],
        "minnbchan": rec["minnbchan"],
        "n1": rec["region1"]["names"],
        "n2": rec["region2"]["names"],
        "fwer1": rec["region1"]["fwer_hit"],
        "fwer2": rec["region2"]["fwer_hit"],
    }, indent=2))
    print(f"wrote {path}")
    plot_sig_maps(rec, OUT)


if __name__ == "__main__":
    main()
