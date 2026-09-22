#!/usr/bin/env python3
"""
Reproduce the exact visited-range moments E[V], Var(V) for every family member
and resolution, writing `tm_candidates_v11.csv` (also mirrored as
`data/tm_moments_manifest.csv`).

This is the moments stage of a full from-scratch reproduction (Sec. "Full
reproduction" in README.md). It recomputes the frozen moments from the family
specification alone, by exact single- and pairwise block-avoidance recursions —
no Monte Carlo, no lesion sampling. Running it reproduces the frozen CSV
bit-for-bit up to floating-point rounding.

Columns: member, resolution, predictor (N_blk = E[V], S_blk = Var(V)), value.
"""
import csv
import numpy as np
from tm_core import N, T, RESOLUTIONS, block_of, MU_A, build_family


def stay(P, keep):
    """P[x_0..x_T all inside the kept block-set] under the reactive walk from
    the uniform start marginal. Exact restricted-operator recursion."""
    A = P * keep[None, :]
    f = keep.astype(float)
    for _ in range(T):
        f = A @ f
    return float(MU_A[keep] @ f[keep])


def moments(P, ell):
    """Exact E[V] and Var(V) of the visited-block count at resolution ell."""
    b = block_of(ell)
    nb = N // ell
    a = np.array([stay(P, b != k) for k in range(nb)])   # P[block k never visited]
    p = 1.0 - a                                           # P[block k visited]
    EV = float(p.sum())
    var = float((p * (1.0 - p)).sum())
    for i in range(nb):
        for j in range(i + 1, nb):
            a_ij = stay(P, (b != i) & (b != j))
            both = 1.0 - a[i] - a[j] + a_ij               # P[i and j visited]
            var += 2.0 * (both - p[i] * p[j])
    return EV, var


def main():
    rows = []
    for name, P in build_family():
        for ell in RESOLUTIONS:
            EV, var = moments(P, ell)
            rows.append(dict(member=name, resolution=ell, predictor="N_blk", value=EV))
            rows.append(dict(member=name, resolution=ell, predictor="S_blk", value=var))
        print(f"{name:18s} done", flush=True)
    with open("tm_candidates_v11.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["member", "resolution", "predictor", "value"])
        w.writeheader(); w.writerows(rows)
    print("wrote tm_candidates_v11.csv")


if __name__ == "__main__":
    main()
