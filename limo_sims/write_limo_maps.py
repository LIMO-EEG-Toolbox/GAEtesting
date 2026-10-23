"""Write one Monte Carlo bundle for MATLAB LIMO (chan x time x subject)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.io import savemat


def to_limo_y(Y_time_chan_sub):
    """Python Y [nTime, nChan, nSub] -> LIMO [nChan, nTime, nSub]."""
    return np.transpose(np.asarray(Y_time_chan_sub, np.float64), (1, 0, 2))


def write_bundle(path, *, Y_planted, Y_control, W, adjacency, labels, planted_ch,
                 planted_idx, control_idx, amp, nboot, maps=None,
                 planted_ch1=None, planted_ch2=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "Y1": to_limo_y(Y_planted),
        "Y2": to_limo_y(Y_control),
        "W": np.asarray(W, np.float64),
        "adjacency": np.asarray(adjacency, np.float64),
        "labels": np.array(labels, dtype=object),
        "planted_ch": np.asarray(planted_ch, np.int32) + 1,
        "planted_idx": np.asarray(planted_idx, np.int32) + 1,
        "control_idx": np.asarray(control_idx, np.int32) + 1,
        "amp": float(amp),
        "nboot": int(nboot),
    }
    if planted_ch1 is not None:
        payload["planted_ch1"] = np.asarray(planted_ch1, np.int32) + 1
    if planted_ch2 is not None:
        payload["planted_ch2"] = np.asarray(planted_ch2, np.int32) + 1
    if maps:
        for name, m in maps.items():
            payload[f"{name}_t"] = np.asarray(m["t"], np.float64)
            payload[f"{name}_p"] = np.asarray(m["p"], np.float64)
            payload[f"{name}_bootT"] = np.asarray(m["bootT"], np.float64)
            payload[f"{name}_bootP"] = np.asarray(m["bootP"], np.float64)
    savemat(path, payload, do_compression=True)
    return path
