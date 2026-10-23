# LIMO GAE simulations

Audit of leave-one-channel-out graph-autoencoder weights at LIMO EEG’s
second level. The weighter is frozen. These scripts plant known hills on
preprocessed ERP CORE P3 maps and read out \(W\) and LIMO cluster-mass FWER.

## What to ship

Copy **this folder** (`limo_sims/`) as a unit (e.g. into
[LIMO-EEG-Toolbox/GAEtesting](https://github.com/LIMO-EEG-Toolbox/GAEtesting)).
It already contains the frozen weighter and the 37 P3 maps.

| Path | Role |
|---|---|
| `limo_sims/` | Python \(t\), matcher, paper runs, MATLAB glue, maps |
| `limo_sims/gae_weighting.py` | **frozen** LOCO weighter. Do not edit. |
| `limo_sims/data/P3_limo/` | 37 first-level maps + STUDY neighbourhood (~57 MB) |

Do **not** ship `gae_audit/` paper-era drivers, the historical `GAEtesting/`
xAI tree, or the Yhat MATLAB weighter (`limo_group_outliers.m` /
`NiPyAEoutliers.py`). MATLAB `limo_tools` is only needed for the 3-channel
smoke; set `LIMO_TOOLS` if that clone is not a sibling of this folder.

## Install

```bash
pip install -r limo_sims/requirements.txt
python limo_sims/mix_people.py            # once: mixed pad from the 37 reals
```

## Data (`load_cohort.py`)

No machine-specific paths. Maps sit **in this folder**:

```
limo_sims/data/P3_limo/
```

(`_paths.py` → `HERE / "data" / "P3_limo"`). You do **not** need the 1.5 GB
preprocessed EEG STUDY. `Beta_files.txt` lists relative paths
(`sub-001/sub-001_desc-Betas.mat`, …). Output: `datain` (26 betas, GAE
training), `Y` (`con_2 − con_1`), STUDY `adjacency` / `coords` / `labels`.

The mixed pad is gitignored (GitHub 100 MB cap). After clone, run
`mix_people.py`; it writes `limo_sims/cache/mixed_pool.npz`. Large-arm
paper runs (`pool.py`) need that file.

MATLAB is optional. Paper figures never open it. The one MATLAB smoke is
`check_threechan_n6.py`.

## LIMO clustering (why these files)

We do **not** call `limo_random_robust`. The two-sample tests themselves
are already in the frozen weighter (`yuen_t`, `weighted_welch_t`; mean is
the same Welch with \(W=1\)). `tmaps.py` only wraps those calls, adds
\(p\)-values, and shuffles labels for the bootstrap. What we still need
from LIMO is the **second-level cluster-mass FWER** (MCC = 2, Time maps):

1. Forming cells: uncorrected \(p \le 0.05\).
2. `limo_findcluster`: drop a channel unless it has `minnbchan` forming
   neighbours at that sample; 4-connected in time; merge across the STUDY
   graph.
3. Mass = \(\sum t^2\) in a component.
4. `limo_cluster_test`: same glue on each label shuffle; \(T\) is the
   \((1-\alpha)\) quantile of the **max** cluster mass. A cluster survives
   if its mass \(\ge T\).

`limo_clustering` in MATLAB **locks `minnbchan = 2`**. A 1- or 2-channel
plant is deleted before FWER. The matcher also needs thousands of inner
FWER calls. Opening MATLAB for each of those is too slow, and cannot turn
the prune off.

So:

- **`limo_mass.py`** — Python copy of that mechanism.
  `findcluster` = `limo_findcluster`; `cluster_test` = `limo_cluster_test`.
  `minnbchan` is an argument (0 or 2). Paper runs and the matcher use this
  so we never launch LIMO on an internal match.
- **`cluster_python.py`** — `cluster_test` with `minnbchan=0` (keep 1- and
  2-channel islands).
- **`cluster_limo.py` + `run_limo_ttest.m`** — the real toolbox on the same
  Python \(t\) maps (`write_limo_maps.py` → `limo_clustering(..., MCC=2)`).
  Default prune only. Kept so a 3-channel plant can still be checked in
  MATLAB (`run_mc.py` when \(n\ge 3\)).
- **`match_limo.py`** — not a clusterer. In the paper we need two plants
  of **different spatial extent** (broad moderate hill vs strong focal
  hole) that are **equally significant under a mean \(t\)-test**, so that
  later GAE weighting is not just comparing a weak effect to a strong
  one. The matcher never uses the GAE. It searches additive amplitudes on
  ordinary mean Welch \(t\) until, under this LIMO FWER, each plant is
  detectable alone (10/10 hits on seeds 0–9) and the two FWER-surviving
  overlap masses (\(\sum t^2\)) match. Peak \(t\) is not matched. Same
  `minnbchan` as the later test; confirmation MCs use new seeds. The GAE
  is applied afterwards, with those amplitudes frozen.

`minnbchan=0` if planted \(n<3\); `minnbchan=2` if \(n\ge 3\). Matcher and
test always share that choice.

We *could* have used the full MATLAB LIMO clusterer whenever the plant
has three or more channels (`minnbchan=2` is LIMO’s default) and kept
the Python copy only for \(n=1\) and \(n=2\), where MATLAB would prune
the island. Paper runs use the Python copy everywhere so matcher
bisections and Monte Carlos do not open MATLAB.

```bash
python limo_sims/check_threechan_n6.py    # MATLAB smoke, 3-channel triangle
python limo_sims/match_limo.py            # matcher CLI, default 3 vs 8, minnbchan=2
python limo_sims/match_1v8_a5.py          # matcher 1 vs 8, minnbchan=0
```

## Other core files

- `gae_weighting.py` — frozen LOCO weighter (do not edit)
- `load_cohort.py` / `pool.py` — 37 real P3 maps; mixed-map pad
- `plant.py` — Gaussian bump / latency on `datain` and \(Y\)
- `tmaps.py` — wrapper: weighter \(t\) + \(p\) + label shuffles
- `run_mc.py` — one MC; Python cluster if \(n<3\), MATLAB if \(n\ge 3\)

## Paper experiments

Constants at the top of each file. GAE never sees arm labels. Plants go on
`datain` and \(Y\). Confirmation MCs start at seed 10; the matcher uses 0–9.
Results: `limo_sims/results/<script>/`. Run in this order (later scripts
read earlier JSON).

A wiring smoke is just those constants: set `N_MC` to 2 or 3 in
`g1_latency.py`, `pz_rare.py`, `iso_bal.py` (paper: 100) and `hitrate.py`
(paper: 200). Cost is almost all one GAE fit per draw. Leave the matcher
alone (`N_MATCH_MC = 10`, `N_HIT = 10` in `match_limo.py`; the Pz floor
uses the same 10/10). Cutting that changes amplitudes. `NBOOT = 101` can
drop to ~21 for a wiring test (noisier FWER). `exist_1v8_a5.py` and
`iso_bal_maps.py` are already one split. Put the paper numbers back before
claiming results.

**1. Whole-scalp delay — is a minority enough?**
`g1_latency.py` · 40 ms delay on every channel of arm A vs 20 unshifted.
\(n_A \in \{20,15,10,6\}\), 100 MC, `minnbchan=2` (the whole scalp is a
clique). Neighbours can copy the delay, so LOCO should not see a hole.
Question: does a small delayed arm still get \(W<1\)? (Paper: no; \(W=1\);
the drop at \(n=6\) is power.) Needs the mixed pad. Figs `g1detect` / `g1W`.

**2. Focal Pz floor — amplitude so the mean can still see it**
`pz_amp_pilot.py` · Pz only, 20 vs 5, **no GAE**, `minnbchan=0`. Same
10/10 mean-FWER bisection as the matcher. Freezes AMP (paper: 65.58).
Must run before the rarity sweep.

**3. Focal Pz rarity — is a hole enough if many maps share it?**
`pz_rare.py` · same AMP, Pz bump, \(n_A\) from 30 down to 1 vs 30, 100 MC,
`minnbchan=0`. Reads `results/pz_amp_pilot/pz_amp_pilot.json`. Question:
when does planted \(W\) on Pz drop? (Paper: only when the hole is rare.)
Figs `pzdet` / `pzW`.

**4. Matcher display — two extents, same mean-test significance**
`match_1v8_a5.py` · 1 vs 8 channels, 5 vs 30, `minnbchan=0`. Writes
`results/match_1v8_a5/match.json` and the three-panel forming/FWER map.
No GAE.

**5. Existence split — both plants on the same five maps, GAE on**
`exist_1v8_a5.py` · frozen matcher split, both hills on the small arm,
weighted FWER contour. Needs `match.json` from (4). Shows the GAE can
drop the focal hole and keep the broad hill on one draw.

**6. Isolation and balance \(W\)**
`iso_bal.py` · 100 MC, rematch per cell.
- Isolation: \(n_1=1..5\) vs 8, 5 vs 30. Shrink the focal set; `minnbchan`
  switches 0 → 2 at \(n_1=3\).
- Balance: 1 vs 8, \(n\) planted vs 30, `minnbchan=0`. Grow how many maps
  carry both plants.
Readout is mean \(W\) on planted vs unplanted cells. Needs mixed pad.

`iso_bal_maps.py` · seed-0 forming/FWER maps for a few isolation and
balance cells (titles: FWER mass). Reuses `iso_bal` match JSON.

**7. Joint isolation × rarity**
`joint_grid.py` · \(n_1=1..5\) vs 8 and \(n_A\) in the balance grid vs 30.
Reuses iso/bal match files. **This runner is frozen at 10 MC**; the paper
figure is the older 100-MC grid, not a rerun of this file. Do not bump
`N_MC` without a new seed block.

**8. Occupancy / H0 mass (cluster-glue check)**
`hitrate.py` · 2 vs 8, 6 vs 30, `minnbchan=0`, 200 MC. Forming and FWER
occupancy on each plant, “chance” cells in clusters that never touch a
plant, and the bootstrap max-mass ECDF. Yuen fattens forming and raises
\(T\); GAE tracks the mean. Needs a match (writes its own if missing).

```bash
python limo_sims/g1_latency.py
python limo_sims/pz_amp_pilot.py
python limo_sims/pz_rare.py
python limo_sims/match_1v8_a5.py
python limo_sims/exist_1v8_a5.py
python limo_sims/iso_bal.py
python limo_sims/iso_bal_maps.py
python limo_sims/joint_grid.py          # 10 MC freeze; paper fig is the old 100-MC
python limo_sims/hitrate.py
```

Temp MATLAB maps: `limo_sims/work/` (gitignored).

## Rules

- Do not edit the LOCO weighter.
- GAE never sees arm labels. Typical = common and neighbour-predictable.
- Plants go on `datain`, not only on \(Y\).
- `minnbchan=0` if planted \(n<3\); `minnbchan=2` if \(n\ge 3\).
- Yuen omitted when planted \(n<6\).
- Confirmation Monte Carlos use seeds from 10; matcher uses 0–9.
