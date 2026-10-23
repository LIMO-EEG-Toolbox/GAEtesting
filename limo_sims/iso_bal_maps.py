#!/usr/bin/env python3
"""Seed-0 forming/FWER maps for isolation (n1=1,2,5) and balance (n=5,10,20).

Reuses iso_bal match JSON. Frozen split seed 0. Plant both on datain.
Python t; LIMO cluster_test (minnbchan as in iso_bal). Mean then weighted.

    python limo_sims/iso_bal_maps.py
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

from iso_bal import (
    BAL_MATCH,
    BAL_MINNB,
    BAL_N1,
    BAL_N2,
    BAL_NB,
    ISO_MATCH,
    ISO_NA,
    ISO_NB,
    ISO_N2,
    TAU,
    load_or_match,
    mean_W,
    minnb_iso,
)
from limo_mass import cluster_test
from match_limo import SPLIT_SEED, _sig_panel, frozen_split
from plant import add_bump
from pool import load_mixed_pad
from tmaps import estimator_boot, groups_from_arms

ISO_MAP_N1 = (1, 2, 5)
BAL_MAP_NA = (5, 10, 20)
NBOOT = 101
OUT = HERE / "results" / "iso_bal"


def t2_mass(tmap, mask):
    on = np.asarray(mask) > 0
    if not np.any(on):
        return 0.0
    return float((np.asarray(tmap, dtype=np.float64) ** 2)[on].sum())


def masses_of(tmap, ct):
    return dict(
        n1=t2_mass(tmap, ct["fwer_n1"]),
        n2=t2_mass(tmap, ct["fwer_n2"]),
        both=float(ct["fwer_overlap_mass"]),
    )


def one_split(cohort, rec, minnbchan):
    s1 = np.asarray(rec["region1"]["channels"], dtype=int)
    s2 = np.asarray(rec["region2"]["channels"], dtype=int)
    a1, a2 = float(rec["a1"]), float(rec["a2"])
    ctx = frozen_split(
        cohort, rec["n_A"], rec["n_B"], SPLIT_SEED, minnbchan, verbose=True,
    )
    d, y = add_bump(ctx["datain"], ctx["Y"], ctx["who_a"], s1, a1)
    d, y = add_bump(d, y, ctx["who_a"], s2, a2)
    W = compute_weights(d, ctx["A"], beta="all", tau=TAU, seeds=(0,))
    groups = groups_from_arms(d.shape[3], ctx["who_a"], ctx["who_b"])
    both = np.concatenate([s1, s2])
    maps = {}
    for i, name in enumerate(("mean", "weighted")):
        ww = np.ones_like(W) if name == "mean" else W
        maps[name] = estimator_boot(
            name, y, groups, ww, NBOOT, np.random.default_rng(SPLIT_SEED + 1 + i),
        )
    cluster = {}
    mass = {}
    for name, m in maps.items():
        cluster[name] = cluster_test(
            m["t"], m["p"], m["bootT"], m["bootP"], ctx["A"],
            planted=both, planted1=s1, planted2=s2,
            minnbchan=minnbchan, return_mask=True,
        )
        mass[name] = masses_of(m["t"], cluster[name])
    return dict(
        rec=rec, ctx=ctx, s1=s1, s2=s2, maps=maps, cluster=cluster,
        W_n1=mean_W(W, s1, ctx["who_a"]),
        W_n2=mean_W(W, s2, ctx["who_a"]),
        mass=mass,
        minnbchan=int(minnbchan),
        a1=a1, a2=a2,
    )


def title_of(lab, pack, name):
    c = pack["cluster"][name]
    m = pack["mass"][name]
    wline = ""
    if name == "weighted":
        wline = f"W n1={pack['W_n1']:.2f} n2={pack['W_n2']:.2f}  "
    return (
        f"{lab}  {name}\n"
        f"{wline}"
        f"mass n1={m['n1']:.0f} n2={m['n2']:.0f}  "
        f"FWER {int(c['n1_any'])}/{int(c['n2_any'])}"
    )


def pair_maps(cases, path, suptitle):
    n = len(cases)
    fig, axes = plt.subplots(2, n, figsize=(4.6 * n, 9.4), squeeze=False)
    meta = []
    im = None
    for j, (lab, pack) in enumerate(cases):
        t_m, p_m = pack["maps"]["mean"]["t"], pack["maps"]["mean"]["p"]
        t_w, p_w = pack["maps"]["weighted"]["t"], pack["maps"]["weighted"]["p"]
        vmax = max(
            float(np.nanmax(np.abs(t_m[p_m <= 0.05]))) if np.any(p_m <= 0.05) else 1.0,
            float(np.nanmax(np.abs(t_w[p_w <= 0.05]))) if np.any(p_w <= 0.05) else 1.0,
            1.0,
        )
        planted = np.concatenate([pack["s1"], pack["s2"]])
        im = _sig_panel(
            axes[0, j], t_m, p_m, pack["cluster"]["mean"]["mask"],
            pack["ctx"]["coords"], pack["ctx"]["names"], planted,
            title_of(lab, pack, "mean"), vmax,
        )
        _sig_panel(
            axes[1, j], t_w, p_w, pack["cluster"]["weighted"]["mask"],
            pack["ctx"]["coords"], pack["ctx"]["names"], planted,
            title_of(lab, pack, "weighted"), vmax,
        )
        if j:
            axes[0, j].set_ylabel("")
            axes[1, j].set_ylabel("")
        cw = pack["cluster"]["weighted"]
        mm, mw = pack["mass"]["mean"], pack["mass"]["weighted"]
        meta.append(dict(
            label=lab, a1=pack["a1"], a2=pack["a2"],
            minnbchan=pack["minnbchan"],
            W_n1=pack["W_n1"], W_n2=pack["W_n2"],
            mean_n1=bool(pack["cluster"]["mean"]["n1_any"]),
            mean_n2=bool(pack["cluster"]["mean"]["n2_any"]),
            weighted_n1=bool(cw["n1_any"]),
            weighted_n2=bool(cw["n2_any"]),
            mean_mass_n1=mm["n1"], mean_mass_n2=mm["n2"],
            weighted_mass_n1=mw["n1"], weighted_mass_n2=mw["n2"],
        ))
    axes[0, 0].set_ylabel("mean  ·  channel (ant → post)")
    axes[1, 0].set_ylabel("weighted  ·  channel (ant → post)")
    fig.suptitle(suptitle, y=0.995)
    fig.tight_layout(rect=(0, 0, 0.96, 0.96))
    cax = fig.add_axes([0.97, 0.12, 0.012, 0.72])
    fig.colorbar(im, cax=cax, label="t (forming p<0.05)")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"wrote {path}", flush=True)
    return meta


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    iso_cohort = load_cohort()
    bal_cohort = load_mixed_pad()

    iso_cases = []
    for n1 in ISO_MAP_N1:
        k = minnb_iso(n1)
        rec = load_or_match(
            ISO_MATCH / f"match_n1{n1}_nA{ISO_NA}.json",
            n1=n1, n2=ISO_N2, n_a=ISO_NA, n_b=ISO_NB,
            minnbchan=k, cohort=iso_cohort,
        )
        pack = one_split(iso_cohort, rec, k)
        iso_cases.append((f"n1={n1}", pack))
        mm, mw = pack["mass"]["mean"], pack["mass"]["weighted"]
        print(
            f"iso n1={n1}  W={pack['W_n1']:.3f}/{pack['W_n2']:.3f}  "
            f"mass mean {mm['n1']:.0f}/{mm['n2']:.0f}  "
            f"weighted {mw['n1']:.0f}/{mw['n2']:.0f}",
            flush=True,
        )

    bal_cases = []
    for n_a in BAL_MAP_NA:
        rec = load_or_match(
            BAL_MATCH / f"match_n1{BAL_N1}_nA{n_a}.json",
            n1=BAL_N1, n2=BAL_N2, n_a=n_a, n_b=BAL_NB,
            minnbchan=BAL_MINNB, cohort=bal_cohort,
        )
        pack = one_split(bal_cohort, rec, BAL_MINNB)
        bal_cases.append((f"n={n_a}", pack))
        mm, mw = pack["mass"]["mean"], pack["mass"]["weighted"]
        print(
            f"bal n={n_a}  W={pack['W_n1']:.3f}/{pack['W_n2']:.3f}  "
            f"mass mean {mm['n1']:.0f}/{mm['n2']:.0f}  "
            f"weighted {mw['n1']:.0f}/{mw['n2']:.0f}",
            flush=True,
        )

    iso_meta = pair_maps(
        iso_cases,
        OUT / "sig_iso_n1_125.png",
        "Isolation  ·  30 vs 5, n2=8  ·  forming p<0.05  ·  frozen split seed 0",
    )
    bal_meta = pair_maps(
        bal_cases,
        OUT / "sig_bal_a5_10_20.png",
        "Balance  ·  1 vs 8, 30 unplanted  ·  forming p<0.05  ·  frozen split seed 0",
    )
    report = dict(
        design=(
            "Seed-0 isolation maps n1=1,2,5 vs 8, 30 vs 5; balance n=5,10,20 vs 30, "
            "n1=1 vs 8. Match cache from iso_bal. Plant both on datain; LOCO W; "
            "Python t; LIMO cluster_test FWER (minnbchan=0 if n1<3 else 2). "
            "Colour = forming p<=0.05; black = FWER cluster."
        ),
        split_seed=SPLIT_SEED,
        isolation=iso_meta,
        balance=bal_meta,
    )
    path = OUT / "sig_pairs.json"
    path.write_text(json.dumps(report, indent=2) + "\n")
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
