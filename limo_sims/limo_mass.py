"""LIMO spatial cluster geometry in Python (Time maps).

Forming p<=alpha, prune channels with fewer than minnbchan significant
neighbours, 4-connected time per channel, merge across STUDY edges.
Mass = sum of t^2 in a component. Inner loop for the LIMO matcher;
MATLAB limo_clustering is the FWER call, not this search.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import label as cc_label

MINNBCHAN = 2
ALPHA = 0.05


def findcluster(onoff, adjacency, minnbchan=MINNBCHAN):
    onoff = np.asarray(onoff, dtype=bool).copy()
    n_chan, n_time = onoff.shape
    select = np.maximum(np.asarray(adjacency), np.asarray(adjacency).T) > 0
    np.fill_diagonal(select, False)
    minnbchan = int(minnbchan)
    if minnbchan > 0:
        sel = select.astype(np.int32)
        nremoved = 1
        while nremoved:
            nsig = sel @ onoff.astype(np.int32)
            remove = onoff & (nsig < minnbchan)
            nremoved = int(remove.sum())
            onoff[remove] = False

    labels = np.zeros((n_chan, n_time), dtype=np.int32)
    total = 0
    struct = np.ones(3, dtype=int)
    for ch in range(n_chan):
        lab, num = cc_label(onoff[ch], structure=struct)
        if num:
            lab[lab > 0] += total
            labels[ch] = lab
            total += int(num)

    parent = np.arange(total + 1, dtype=np.int32)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return int(x)

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra < rb:
            parent[rb] = ra
        elif rb < ra:
            parent[ra] = rb

    for i in range(n_chan):
        for j in np.flatnonzero(select[i]):
            if int(j) <= i:
                continue
            both = (labels[i] > 0) & (labels[int(j)] > 0)
            for a, b in zip(labels[i, both], labels[int(j), both]):
                union(int(a), int(b))

    remap = {}
    out = np.zeros_like(labels)
    n = 0
    for i in range(n_chan):
        for t in range(n_time):
            v = int(labels[i, t])
            if v == 0:
                continue
            r = find(v)
            if r not in remap:
                n += 1
                remap[r] = n
            out[i, t] = remap[r]
    return out, n


def cluster_masses(tmap, labels, n):
    stat = np.asarray(tmap, dtype=np.float64) ** 2
    masses = np.zeros(n, dtype=np.float64)
    for c in range(1, n + 1):
        masses[c - 1] = stat[labels == c].sum()
    return masses


def max_t2_mass(tmap, pmap, adjacency, alpha=ALPHA, minnbchan=MINNBCHAN):
    labels, n = findcluster(np.asarray(pmap) <= alpha, adjacency, minnbchan)
    if n == 0:
        return 0.0
    return float(cluster_masses(tmap, labels, n).max())


def planted_overlap(tmap, pmap, adjacency, planted, alpha=ALPHA, minnbchan=MINNBCHAN):
    labels, n = findcluster(np.asarray(pmap) <= alpha, adjacency, minnbchan)
    empty = dict(mass=0.0, n_cells=0, n_chans=0, n_clusters=n)
    if n == 0:
        return empty
    planted = np.asarray(planted, dtype=int)
    stat = np.asarray(tmap, dtype=np.float64) ** 2
    best_c, best_m = 0, 0.0
    for c in range(1, n + 1):
        if np.any(labels[planted] == c):
            m = float(stat[labels == c].sum())
            if m > best_m:
                best_c, best_m = c, m
    if best_c == 0:
        return empty
    mask = labels == best_c
    return dict(
        mass=best_m,
        n_cells=int(mask.sum()),
        n_chans=int(mask.any(1).sum()),
        n_clusters=n,
    )


def planted_overlap_t2_mass(tmap, pmap, adjacency, planted, alpha=ALPHA, minnbchan=MINNBCHAN):
    return planted_overlap(tmap, pmap, adjacency, planted, alpha, minnbchan)["mass"]


def _sig_on(mask, chans):
    chans = np.asarray(chans, dtype=int)
    if mask.size == 0 or chans.size == 0:
        return False, 0
    sl = mask[chans]
    return bool(np.any(sl)), int(np.count_nonzero(sl))


def overlap_bool(labels, n, tmap, planted):
    """Largest forming (or FWER-labelled) component that touches planted."""
    shape = np.asarray(labels).shape
    empty = np.zeros(shape, dtype=bool)
    planted = np.asarray(planted, dtype=int)
    if n == 0 or planted.size == 0:
        return empty, 0.0
    stat = np.asarray(tmap, dtype=np.float64) ** 2
    best_c, best_m = 0, 0.0
    for c in range(1, n + 1):
        if np.any(labels[planted] == c):
            m = float(stat[labels == c].sum())
            if m > best_m:
                best_c, best_m = c, m
    if best_c == 0:
        return empty, 0.0
    return labels == best_c, best_m


def chance_bool(labels, n, planted):
    """Forming/FWER cells whose component never touches planted."""
    planted = np.asarray(planted, dtype=int)
    out = np.zeros(np.asarray(labels).shape, dtype=bool)
    if n == 0:
        return out
    touch = set()
    for c in range(1, n + 1):
        if planted.size and np.any(labels[planted] == c):
            touch.add(c)
    for c in range(1, n + 1):
        if c not in touch:
            out |= labels == c
    return out


def cluster_test(tmap, pmap, bootT, bootP, adjacency, planted=None,
                 planted1=None, planted2=None, alpha=ALPHA, minnbchan=MINNBCHAN,
                 return_mask=False):
    """LIMO limo_cluster_test: observed clusters with mass >= boot max-mass threshold."""
    nboot = int(bootT.shape[2])
    boot_max = np.array(
        [
            max_t2_mass(bootT[:, :, b], bootP[:, :, b], adjacency,
                        alpha=alpha, minnbchan=minnbchan)
            for b in range(nboot)
        ],
        dtype=np.float64,
    )
    order = np.sort(boot_max)
    u = int(np.round((1.0 - alpha) * nboot))
    u = min(max(u, 1), nboot)
    max_th = float(order[u - 1])

    labels, n = findcluster(np.asarray(pmap) <= alpha, adjacency, minnbchan)
    masses = cluster_masses(tmap, labels, n) if n else np.zeros(0)
    mask = np.zeros(np.asarray(tmap).shape, dtype=np.int32)
    cluster_p = []
    k = 0
    for c in range(1, n + 1):
        m = float(masses[c - 1])
        if m >= max_th:
            k += 1
            mask[labels == c] = k
            p = 1.0 - float(np.sum(m >= order)) / nboot
            cluster_p.append(max(p, 1.0 / nboot))

    planted = np.array([] if planted is None else planted, dtype=int)
    p1 = np.array([] if planted1 is None else planted1, dtype=int)
    p2 = np.array([] if planted2 is None else planted2, dtype=int)
    if planted.size == 0 and p1.size:
        planted = np.concatenate([p1, p2]) if p2.size else p1
    any_p, n_p = _sig_on(mask, planted)
    any1, n1 = _sig_on(mask, p1)
    any2, n2s = _sig_on(mask, p2)
    ov = planted_overlap(tmap, pmap, adjacency, planted, alpha, minnbchan) if planted.size else dict(
        mass=0.0, n_cells=0, n_chans=0, n_clusters=n
    )
    stat = np.asarray(tmap, dtype=np.float64) ** 2
    fwer_ov = 0.0
    if planted.size and k:
        for lab in range(1, k + 1):
            if np.any(mask[planted] == lab):
                fwer_ov = max(fwer_ov, float(stat[mask == lab].sum()))
    out = {
        "n_sig_cells": int(np.count_nonzero(mask)),
        "max_th": max_th,
        "obs_max_mass": float(masses.max()) if n else 0.0,
        "boot_max": boot_max,
        "min_cluster_p": None if not cluster_p else float(min(cluster_p)),
        "n_sig_clusters": k,
        "planted_any_cluster": any_p,
        "planted_n_sig": n_p,
        "n1_any": any1,
        "n1_n_sig": n1,
        "n2_any": any2,
        "n2_n_sig": n2s,
        "forming_n_clusters": int(n),
        "forming_overlap_mass": float(ov["mass"]),
        "fwer_overlap_mass": float(fwer_ov),
        "minnbchan": int(minnbchan),
    }
    if return_mask:
        both = planted
        form_n1, _ = overlap_bool(labels, n, tmap, p1)
        form_n2, _ = overlap_bool(labels, n, tmap, p2)
        fwer_n = int(k)
        fwer_n1, _ = overlap_bool(mask, fwer_n, tmap, p1)
        fwer_n2, _ = overlap_bool(mask, fwer_n, tmap, p2)
        out["mask"] = mask
        out["forming_labels"] = labels
        out["form_n1"] = form_n1
        out["form_n2"] = form_n2
        out["form_chance"] = chance_bool(labels, n, both)
        out["fwer_n1"] = fwer_n1
        out["fwer_n2"] = fwer_n2
        out["fwer_chance"] = chance_bool(mask, fwer_n, both)
    return out
