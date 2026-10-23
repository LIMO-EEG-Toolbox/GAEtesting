"""Python two-sample t (mean / Yuen / LOCO-weighted Welch) and label shuffles."""
from __future__ import annotations

import numpy as np
from scipy import stats

import _paths  # noqa: F401
from gae_weighting import weighted_welch_t, yuen_t

ESTIMATORS = ("mean", "yuen", "weighted")


def groups_from_arms(n_sub, planted, control):
    g = np.full(n_sub, -1, dtype=int)
    g[np.asarray(control, dtype=int)] = 0
    g[np.asarray(planted, dtype=int)] = 1
    return g


def t_and_p(tmap, df):
    t = np.asarray(tmap, dtype=np.float64)
    d = np.asarray(df, dtype=np.float64)
    p = np.ones(t.shape, dtype=np.float64)
    ok = np.isfinite(t) & np.isfinite(d) & (d > 1e-6)
    p[ok] = 2.0 * stats.t.sf(np.abs(t[ok]), d[ok])
    t = np.where(ok, t, 0.0)
    d = np.where(ok, d, np.nan)
    return t, p, d


def yuen_df_map(groups, shape, tr=0.2):
    n0 = int((groups == 0).sum())
    n1 = int((groups == 1).sum())
    h0 = n0 - 2 * int(tr * n0)
    h1 = n1 - 2 * int(tr * n1)
    return np.full(shape, float(h0 + h1 - 2))


def estimator_boot(name, Y, groups, W, nboot, rng):
    ones = np.ones_like(W)

    def draw(g):
        if name == "yuen":
            tmap, _ = yuen_t(Y, g)
            return t_and_p(tmap, yuen_df_map(g, tmap.shape))
        ww = ones if name == "mean" else W
        tmap, df, _ = weighted_welch_t(Y, g, ww)
        return t_and_p(tmap, df)

    t_obs, p_obs, df_obs = draw(groups)
    n_chan, n_time = t_obs.shape
    bootT = np.zeros((n_chan, n_time, nboot), np.float64)
    bootP = np.zeros((n_chan, n_time, nboot), np.float64)
    bootDF = np.zeros((n_chan, n_time, nboot), np.float64)
    idx = np.flatnonzero(groups >= 0)
    g0 = groups[idx]
    for b in range(nboot):
        gb = g0.copy()
        rng.shuffle(gb)
        g = np.full_like(groups, -1)
        g[idx] = gb
        tb, pb, dfb = draw(g)
        bootT[:, :, b] = tb
        bootP[:, :, b] = pb
        bootDF[:, :, b] = dfb
    return dict(t=t_obs, p=p_obs, df=df_obs, bootT=bootT, bootP=bootP, bootDF=bootDF)


def all_estimator_boots(Y, groups, W, nboot, rng0):
    maps = {}
    for i, name in enumerate(ESTIMATORS):
        maps[name] = estimator_boot(name, Y, groups, W, nboot, np.random.default_rng(int(rng0) + i))
    return maps
