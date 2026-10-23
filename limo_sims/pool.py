"""Mixed-map pad used in the paper (37 real + 37 mix).

Build once after clone:

    python limo_sims/mix_people.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from _paths import POOL


def load_mixed_pad():
    path = Path(POOL)
    if not path.is_file():
        raise FileNotFoundError(
            f"mixed pad missing: {path}\n"
            "Run: python limo_sims/mix_people.py"
        )
    pool = np.load(path, allow_pickle=True)
    datain_all, Y_all, A = pool["datain"], pool["Y"], pool["adjacency"]
    n_real = len(pool["real_subjects"])
    n_mix = datain_all.shape[3] - n_real
    rng = np.random.default_rng(1)
    mix_keep = np.sort(rng.choice(np.arange(n_real, n_real + n_mix), n_real, replace=False))
    keep = np.concatenate([np.arange(n_real), mix_keep])
    labels = [str(x).strip() for x in np.asarray(pool["labels"]).tolist()]
    return dict(
        datain=datain_all[:, :, :, keep],
        Y=Y_all[:, :, keep],
        adjacency=A,
        labels=labels,
        coords=np.asarray(pool["coords"], dtype=np.float64),
    )
