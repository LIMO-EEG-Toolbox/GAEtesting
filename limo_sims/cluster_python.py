"""LIMO cluster FWER with minnbchan=0 (1- and 2-channel islands kept).

Same test as limo_cluster_test: forming p<=0.05, mass = sum t^2, threshold =
bootstrap max cluster mass. Neighbour prune is off. Not MATLAB limo_clustering.
"""
from __future__ import annotations

from limo_mass import cluster_test

MINNBCHAN = 0


def cluster_maps(maps, adjacency, planted_ch, minnbchan=MINNBCHAN,
                 planted_ch1=None, planted_ch2=None):
    return {
        name: cluster_test(
            m["t"], m["p"], m["bootT"], m["bootP"], adjacency,
            planted=planted_ch, planted1=planted_ch1, planted2=planted_ch2,
            minnbchan=minnbchan,
        )
        for name, m in maps.items()
    }
