#!/usr/bin/env python3
"""
Analytic protocol-failure assessment — unperturbed quantities only.

NO lesion run. NO registered removal seed is instantiated. Everything computed
here is a functional of the unperturbed reactive kernel, which §5 already
authorizes as predictor computation.

  (a) diagonal positivity  -> whether a dead row can exist at all
  (b) E[V(ell)]            -> expected number of distinct blocks visited
  (c) X(ell)               -> minimum trajectory crossing probability over
                              balanced contiguous block cuts
"""

import numpy as np, csv, time
from tm_core import build_family, N, T, RESOLUTIONS, block_of, MU_A


def stay_prob(P, keep):
    """P[x_0..x_T all inside `keep`] under the unperturbed reactive walk
    started from mu_A. One restricted-operator recursion, exact."""
    A = P * keep[None, :]
    f = keep.astype(float)
    for _ in range(T):
        f = A @ f
    return float(MU_A[keep] @ f[keep])


def expected_distinct_blocks(P, ell):
    """E[V] = sum_b (1 - a_b), a_b = P[block b never visited].
    Exact: one recursion per block, no lumpability required."""
    b = block_of(ell)
    nb = N // ell
    ev = 0.0
    for blk in range(nb):
        ev += 1.0 - stay_prob(P, b != blk)
    return ev


def min_cut_crossing(P, ell):
    """X = min over balanced contiguous block arcs U of
       1 - P[stay in U] / mu(U).
    A cut statistic of the reactive current over the whole T-layer horizon.
    Unlike the A->B cut of §5 P4 it is a path event, not a flux balance, so it
    is not pinned by conservation of reactive mass."""
    b = block_of(ell)
    nb = N // ell
    half = max(1, nb // 2)
    best = np.inf
    for start in range(nb):
        arc = set((start + j) % nb for j in range(half))
        keep = np.isin(b, list(arc))
        mu_U = keep.sum() / N
        best = min(best, 1.0 - stay_prob(P, keep) / mu_U)
    return float(best)


def main():
    fam = build_family()

    print("=" * 78)
    print("(a) dead-row feasibility: a row x in S dies iff supp(P^R(x,.)) is")
    print("    contained in the removed set. P^R(x,x) > 0 makes that impossible")
    print("    for every removal set, because x itself survives.")
    print("=" * 78)
    print(f"{'member':18s} {'min diag':>12s}  {'dead row possible?':>20s}")
    n_pos = 0
    for name, P in fam:
        d = float(np.min(np.diag(P)))
        possible = d <= 0.0
        n_pos += (not possible)
        print(f"{name:18s} {d:12.8f}  {'YES' if possible else 'no — never':>20s}")
    print(f"\n  {n_pos} of {len(fam)} members have a strictly positive diagonal")
    print(f"  -> R(eta)/R0 = |S|/n identically, hence eta_half = 1/2 exactly,")
    print(f"     at every resolution, for those {n_pos} members.")

    print()
    print("=" * 78)
    print("(b),(c) candidate predictors, unperturbed, exact for every member")
    print("=" * 78)
    rows = []
    t0 = time.time()
    for name, P in fam:
        line = f"{name:18s}"
        for ell in RESOLUTIONS:
            ev = expected_distinct_blocks(P, ell)
            rows.append(dict(member=name, resolution=ell,
                             predictor="N_blk", value=ev))
            if ell in (16, 4):
                line += f"  E[V](l={ell})={ev:9.4f}"
        for ell in RESOLUTIONS:
            xc = min_cut_crossing(P, ell)
            rows.append(dict(member=name, resolution=ell,
                             predictor="X_cut", value=xc))
            if ell in (16, 4):
                line += f"  X(l={ell})={xc:8.6f}"
        print(line)
    print(f"  ({time.time()-t0:.0f}s)")

    with open("tm_candidates_v11.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["member", "resolution",
                                          "predictor", "value"])
        w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
