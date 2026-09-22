# Moment compression of visited-range distributions for extensive-lesion robustness in layered Markov chains

Reproducibility package for the preregistered study of extensive-lesion robustness
in feed-forward Markov chains. It contains the **frozen preregistered protocol**, the
**analysis code**, the **frozen data outputs**, and **SHA-256 manifests**, so that
every reported result and figure can be reproduced — either quickly from the
frozen outputs, or fully from scratch by rerunning the pipeline.

- Central result: the disorder-averaged dead-end robustness curve is an exact
  functional of the visited-range distribution P(V),
  `E_D[R_dead(η)]/R0 = E_ω[ C(n_b−V, r) / C(n_b, r) ]`,
  and across the preregistered 27-member × 5-resolution family the mean+variance
  closure reproduces the curve to 1% in 105 of 120 mean-insufficient cells.

## What is frozen

Nothing in `protocol/` or the listed frozen data files is regenerated when the
repository is assembled. The scientific record is fixed by:

- `protocol/PROTOCOL_v1.1.md` — the frozen preregistered protocol
  (SHA-256 `4a34a4e1dea12adbb94bb5c067e1691a34d4301632f4600974ab890d2660a9dd`).
- `protocol/PROTOCOL_v1.1_FREEZE_RECORD.md` — freeze record and review trail.
- `protocol/PROTOCOL_v1.0.md`, `protocol/PROTOCOL_v1.0_SUPERSEDED.md` — the earlier
  version, **superseded before any outcome was computed**, kept for the full
  preregistration trail.
- `MANIFEST.sha256` — SHA-256 of every file in the repository.

The seeds and sizes are registered in `protocol/PROTOCOL_v1.1.md` (the registered protocol lists them): n=256, T=64, λ=1−2^(−1/64), resolutions ℓ∈{1,2,4,8,16};
family random-permutation seeds, per-member trajectory seed `2026400000 + 1e4·index`
(M=10^6), and the validation removal seed.

## Requirements

Python 3.11, `numpy`, `scipy`, `matplotlib`, `networkx` (see `requirements.txt`):

```
pip install -r requirements.txt
```

All commands below are run from the repository root.

## Quick reproduction (from frozen outputs, seconds)

Rebuilds the endpoint table, the family verdict, and the figures **from the frozen
trajectory histograms** `data/tm_PV_hat.npz` and the frozen moments
`tm_candidates_v11.csv` — without the M=10^6 sampling.

```
python tm_outcome_analysis.py     # E2 / E3 / E6 / E1  ->  data/tm_E2_table.csv, data/tm_E2_summary.json
python tm_outcome_final.py        # E7, secondary (analytic), removal-set validation
python make_figs.py               # Fig. 2 (regime map) and Fig. 3 (representative curves)
python make_fig1.py               # Fig. 1 schematic
```

Expected: `data/tm_E2_table.csv` reproduces the frozen table (family verdict
`H_dist`; 15 mean-sufficient cells; 105/120 close at second order; 14 improve-not-
closed; 1 higher-order), and `fig2_regime_map.pdf`, `fig3_curves.pdf` reproduce the reported figures.

## Full reproduction (from scratch)

Regenerates the moments, the trajectory histograms, and all checks, then the
endpoints and figures. Deterministic under the registered seeds. Total runtime is a
few minutes on a laptop (the M=10^6 sampling dominates, ~2–3 min).

```
# 1. code checks (exact, no sampling)
python tm_T6.py                   # controls K-a, K-b        (T6)
python tm_T8.py                   # oracle curve validation  (T8)

# 2. exact moments  ->  tm_candidates_v11.csv
python tm_moments.py

# 3. trajectory histograms P̂(V), M=10^6, registered seeds  ->  data/tm_PV_hat.npz
python tm_sample_PV.py

# 4. endpoints and figures (as in Quick reproduction)
python tm_outcome_analysis.py
python tm_outcome_final.py
python make_figs.py
python make_fig1.py
```

`toy_verify.py`, `assess_analytic.py`, and `feasibility_PV.py` reproduce the
analytic identity check on the small enumerable chain, the diagonal/redirect
degeneration analysis, and the exact-vs-sampled P(V) feasibility check, respectively.

### Verifying against the frozen outputs

Quantities recomputed in a full run should match the frozen files to floating point
(the M=10^6 estimator agrees within its ≤5×10^(−4) pointwise standard error, so the
`H_dist` verdict, the 105/120 count, and the regime map are stable). Frozen
reference hashes:

| file | SHA-256 |
|---|---|
| `tm_candidates_v11.csv` (= `data/tm_moments_manifest.csv`) | `268b95ab1635fb452cd5636c081b250e0429fa6d35d7cf6d395c632f400e61c1` |
| `data/tm_PV_hat.npz` | `d31f823402b725c041cad6a664999c7cbb8dde31263a083a1f393b47547832ba` |
| `data/tm_E2_table.csv` | `5ae3d9eb2368c12f033c062b6f0e9f45da1f6e169ad6d10fdd2c9163e2e395b6` |
| `protocol/PROTOCOL_v1.1.md` | `4a34a4e1dea12adbb94bb5c067e1691a34d4301632f4600974ab890d2660a9dd` |

`MANIFEST.sha256` lists all files; verify with `sha256sum -c MANIFEST.sha256`.

## Layout

```
protocol/   frozen preregistered protocol (v1.1), freeze record, superseded v1.0
tm_*.py     analysis modules and pipeline entry points  (run from repo root)
data/       frozen outputs: P̂(V), moments manifest, E2 table, summary
tm_candidates_v11.csv, tm_checks_v11.json, tm_predictors_v11.csv   frozen Phase-1/2 records
```

The pipeline modules read their inputs by the relative paths used throughout the analysis
(`tm_candidates_v11.csv`, `data/…`), so they are run from the repository root.

## Family and endpoints in brief

27 doubly stochastic kernels (Loc, Perm, Mix, Block, ShiftMix), uniform reactive
marginals, uniform-leak embedding; persistent dead-end lesions of a random η-fraction
of coarse-grained blocks. Truth curves are exact for the 13 oracle members
(permutations by enumeration, Mix by a stable occupancy recursion) and, for the
other 14, the disorder-integrated trajectory estimator. Closures are maximum-entropy
laws under a mean, or a mean-and-variance, constraint; the decision statistic is the
sup-norm curve distance D_{k,j}(ℓ) with registered tolerance τ=10^(−2).

## License

Code: MIT (`LICENSE`). Data, protocol, and figures: CC-BY-4.0 (`LICENSE-DATA`).

## Citing

See `CITATION.cff`. Please cite the associated article and, for the exact frozen release, the
archived DOI (mint one via Zenodo on the tagged release; see `.zenodo.json`).
