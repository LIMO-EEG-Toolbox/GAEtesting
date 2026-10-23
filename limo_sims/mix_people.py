#!/usr/bin/env python3
"""Build mixed maps from the 37 reals (pool = 100). No plant, no GAE.

Each daughter is the equal average of three distinct parents.
Triples are a uniform random sample of unordered combinations.

    python limo_sims/mix_people.py
"""
from itertools import combinations
from pathlib import Path

import numpy as np

import _paths  # noqa: F401
from _paths import HERE
from load_cohort import load_cohort

N_MIX = 63
N_PARENTS = 3
SEED = 1
OUT = HERE / "cache" / "mixed_pool.npz"


def random_triples(n_real, n_mix, rng):
    all_t = np.array(list(combinations(range(n_real), N_PARENTS)), dtype=int)
    pick = rng.choice(len(all_t), size=n_mix, replace=False)
    return all_t[pick]


def mix_axis(arr, triples):
    a, b, c = triples[:, 0], triples[:, 1], triples[:, 2]
    return (arr[..., a] + arr[..., b] + arr[..., c]) / 3.0


def main():
    cohort = load_cohort()
    datain, Y = cohort["datain"], cohort["Y"]
    names = np.array(cohort["subjects"])
    n_real = datain.shape[3]
    rng = np.random.default_rng(SEED)
    triples = random_triples(n_real, N_MIX, rng)

    mixed_datain = mix_axis(datain, triples).astype(np.float32)
    mixed_Y = mix_axis(Y, triples).astype(np.float32)
    pool_datain = np.concatenate([datain, mixed_datain], axis=3)
    pool_Y = np.concatenate([Y, mixed_Y], axis=2)

    uses = np.bincount(triples.ravel(), minlength=n_real)
    parent_ids = names[triples]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        OUT,
        datain=pool_datain,
        Y=pool_Y,
        mixed_datain=mixed_datain,
        mixed_Y=mixed_Y,
        triples=triples,
        parent_ids=parent_ids,
        real_subjects=names,
        adjacency=cohort["adjacency"],
        coords=cohort["coords"],
        labels=np.array(cohort["labels"]),
    )

    print(f"real {n_real}  mixed {N_MIX}  pool {pool_datain.shape[3]}")
    print(f"wrote {OUT}")
    print(f"parent uses: min={uses.min()} max={uses.max()} mean={uses.mean():.2f}")


if __name__ == "__main__":
    main()
