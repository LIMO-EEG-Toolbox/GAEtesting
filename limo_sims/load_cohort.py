"""P3 first-level maps -> tensors the weighter already understands.

datain  [nBeta, nTime, nChan, nSub]   all condition betas (GAE training)
Y       [nTime, nChan, nSub]          con_2 - con_1 (the t-test map)
A       [nChan, nChan]                STUDY neighbour graph
"""
from pathlib import Path

import h5py
import numpy as np
import scipy.io as sio

import _paths  # noqa: F401
from _paths import DATA as ROOT

BETA_LIST = ROOT / "Beta_files.txt"
CHANLOCS = ROOT / "limo_chanlocs.mat"


def _label(x):
    if isinstance(x, bytes):
        return x.decode()
    if isinstance(x, np.ndarray):
        return _label(x.item()) if x.size == 1 else "".join(_label(v) for v in x.flat)
    return str(x).strip()


def _subject(beta_path):
    for p in beta_path.parents:
        if p.name.startswith("sub-"):
            return p.name
    return beta_path.stem.split("_")[0]


def load_cohort(root=ROOT):
    files = []
    for ln in BETA_LIST.read_text().splitlines():
        if not ln.strip():
            continue
        path = Path(ln.strip())
        if not path.is_absolute():
            path = root / path
        files.append(path)
    maps = []
    for path in files:
        with h5py.File(path, "r") as f:
            maps.append(np.asarray(f["Betas"], np.float32))  # [nBeta, nTime, nChan]
    datain = np.stack(maps, axis=3)

    loc = sio.loadmat(CHANLOCS, squeeze_me=False, struct_as_record=False)
    A = np.asarray(loc["channeighbstructmat"], np.float32)
    np.fill_diagonal(A, 0.0)
    A = np.maximum(A, A.T)
    ch = loc["expected_chanlocs"]
    labels = [_label(ch[0, i].labels) for i in range(ch.shape[1])]
    coords = np.array(
        [[float(np.asarray(ch[0, i].X)),
          float(np.asarray(ch[0, i].Y)),
          float(np.asarray(ch[0, i].Z))] for i in range(ch.shape[1])],
        dtype=np.float64,
    )

    nTime, nChan, nSub = datain.shape[1], datain.shape[2], datain.shape[3]
    Y = np.empty((nTime, nChan, nSub), np.float32)
    subjects = []
    for s, path in enumerate(files):
        sub = _subject(path)
        subjects.append(sub)
        c2 = sio.loadmat(path.with_name(f"{sub}_desc-con_2.mat"))["con"][:, :, 0]
        c1 = sio.loadmat(path.with_name(f"{sub}_desc-con_1.mat"))["con"][:, :, 0]
        Y[:, :, s] = (c2 - c1).T.astype(np.float32)

    return dict(
        datain=datain, Y=Y, adjacency=A, coords=coords,
        labels=labels, subjects=subjects,
    )
