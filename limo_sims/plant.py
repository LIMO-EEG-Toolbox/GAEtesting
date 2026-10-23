"""Additive Gaussian bump on 1–3 STUDY channels. Plant on datain and Y."""
from __future__ import annotations

import numpy as np

DT_MS = 4.0
T0_MS = -200.0
PEAK_MS = 380.0
WIDTH_MS = 40.0
SEED_CHAN = "Pz"


def times_ms(n_time):
    return T0_MS + np.arange(n_time, dtype=np.float64) * DT_MS


def gauss_bump(n_time, amp):
    t = times_ms(n_time)
    return amp * np.exp(-0.5 * ((t - PEAK_MS) / WIDTH_MS) ** 2).astype(np.float32)


def labels_list(labels):
    return [str(x).strip() for x in np.asarray(labels).tolist()]


def neighbour_pair(adjacency, labels, seed_name=SEED_CHAN):
    return planted_channels(adjacency, labels, 2, seed_name=seed_name)


def planted_channels(adjacency, labels, n, seed_name=SEED_CHAN):
    """n=1 seed; n=2 seed+neighbour; n=3 a STUDY triangle on the seed."""
    n = int(n)
    names = labels_list(labels)
    seed = names.index(seed_name)
    A = np.asarray(adjacency)
    A = np.maximum(A, A.T)
    np.fill_diagonal(A, 0.0)
    if n == 1:
        chans = [seed]
    elif n == 2:
        neigh = np.flatnonzero(A[seed] > 0)
        if neigh.size < 1:
            raise ValueError(f"{seed_name} has no STUDY neighbours")
        chans = [seed, int(neigh[0])]
    elif n == 3:
        neigh = np.flatnonzero(A[seed] > 0)
        chans = None
        for i, a in enumerate(neigh):
            for b in neigh[i + 1 :]:
                if A[int(a), int(b)] > 0:
                    chans = [seed, int(a), int(b)]
                    break
            if chans is not None:
                break
        if chans is None:
            raise ValueError(f"{seed_name} has no STUDY triangle")
    else:
        raise ValueError("n_planted is 1, 2, or 3")
    return np.array(chans, dtype=int), [names[c] for c in chans]


def add_bump(datain, Y, who, chans, amp):
    """who: subject indices. chans: channel indices. Same amp on each."""
    bump = gauss_bump(Y.shape[0], amp)
    d = datain.copy()
    y = Y.copy()
    who = np.asarray(who, dtype=int)
    for ch in np.asarray(chans, dtype=int):
        d[:, :, ch, who] += bump[None, :, None]
        y[:, ch, who] += bump[:, None]
    return d, y


def draw_arms(n_sub, n_a, n_b, rng):
    take = rng.permutation(n_sub)[: n_a + n_b]
    planted = np.sort(take[:n_a])
    control = np.sort(take[n_a:])
    return planted, control
