#!/usr/bin/env python3
"""
PROTOCOL v1.1 §7 T6 — controls K-a, K-b (§3.4).

NOT a lesion/outcome computation: no removal set, no eta, no curve, no
registered removal seed. Only the unperturbed visited-block moments E[V], Var(V)
are compared between the base Loc(1) chain and its two control augmentations,
via the exact block-avoidance recursion. Tolerance 1e-12.

K-a  Loc(1) with an independent coordinate z in {0..15} resampled uniformly
     every layer; blocks are on x only. V must be unchanged.
K-b  Loc(1) with a gate at layer T/2: with prob 1/2, eight identity waiting
     layers are inserted. Identity layers visit no new x-block, so V unchanged.

Both are exact-invariance controls: the augmented recursion is expected to
reproduce the base recursion to machine precision.
"""

import numpy as np
import time
from tm_core import N, T, MU_A, block_of, k_loc

NZ = 16          # auxiliary coordinate range for K-a, {0..15}  (§3.4)
GATE_WAIT = 8    # identity waiting layers for K-b                (§3.4)
TOL = 1e-12
P1 = k_loc(1)    # base member: Loc(1)


# ------------------------------------------------------------------ base
def stay_base(keep):
    """P[x_0..x_T all in keep-block-set] under Loc(1), start uniform. Exact."""
    A = P1 * keep[None, :]          # restrict columns to survivors
    f = keep.astype(float)
    for _ in range(T):
        f = A @ f
    return float(MU_A[keep] @ f[keep])


# ------------------------------------------------------------------ K-a
def stay_Ka(keep):
    """Same probability on the augmented chain (x ~ Loc(1), z ~ uniform resample
    each layer), state carried explicitly as an (N x NZ) array. Block depends on
    x only, so keep is broadcast over z. The recursion is run genuinely on the
    augmented state; if the control is correct it collapses to stay_base."""
    keepf = keep.astype(float)
    F = np.repeat(keepf[:, None], NZ, axis=1)      # F0[x,z] = keep(x)
    for _ in range(T):
        rowmean = F.mean(axis=1)                   # marginalize resampled z'
        h = P1 @ rowmean                           # x-move by Loc(1)
        F = (keepf * h)[:, None] * np.ones((1, NZ))  # restrict, broadcast over z
    mu_aug = (MU_A[:, None] / NZ)                  # uniform over (x,z)
    return float((mu_aug * F)[keep].sum())


# ------------------------------------------------------------------ K-b
def stay_Kb(keep):
    """Loc(1) with the bimodal timing gate at layer T/2. Backward recursion:
    second base segment, then the gate operator 1/2 (identity^8 restricted) +
    1/2 (skip), then first base segment. Identity restricted to keep is
    idempotent, so the gate is expected to be a no-op on V."""
    A = P1 * keep[None, :]
    keepf = keep.astype(float)
    g = keepf.copy()
    half = T // 2
    for _ in range(half):                          # second base segment
        g = A @ g
    wait = g.copy()                                # 8 identity waiting layers
    for _ in range(GATE_WAIT):
        wait = keepf * wait                        # identity restricted to keep
    g = 0.5 * wait + 0.5 * g                       # bimodal gate, equiprobable
    for _ in range(half):                          # first base segment
        g = A @ g
    return float(MU_A[keep] @ g[keep])


# ------------------------------------------------------------ moments
def moments(stay_fn, ell, with_var=True):
    b = block_of(ell)
    nb = N // ell
    a = np.array([stay_fn(b != k) for k in range(nb)])
    p = 1.0 - a
    EV = float(p.sum())
    if not with_var:
        return EV, None
    var = float((p * (1 - p)).sum())
    for i in range(nb):
        for j in range(i + 1, nb):
            a_ij = stay_fn((b != i) & (b != j))
            both = 1.0 - a[i] - a[j] + a_ij
            var += 2.0 * (both - p[i] * p[j])
    return EV, var


def main():
    # full E[V]+Var at coarse resolutions; E[V]-only at the two finest
    # (n_b = 128, 256 pairwise is not needed to establish the invariance).
    plan = {16: True, 8: True, 4: True, 2: False, 1: False}

    print("=" * 74)
    print("T6 controls K-a, K-b — visited-block moments vs base Loc(1)")
    print(f"tolerance {TOL:.0e}   NZ={NZ}  GATE_WAIT={GATE_WAIT}  T={T}")
    print("=" * 74)
    print(f"{'ell':>4} {'quantity':>9} {'base':>14} {'K-a':>14} {'K-b':>14}"
          f" {'|Ka-b|':>9} {'|Kb-b|':>9}")

    worst = 0.0
    for ell, wv in plan.items():
        t0 = time.time()
        eb, vb = moments(stay_base, ell, wv)
        ea, va = moments(stay_Ka, ell, wv)
        ek, vk = moments(stay_Kb, ell, wv)
        dEa, dEk = abs(ea - eb), abs(ek - eb)
        worst = max(worst, dEa, dEk)
        print(f"{ell:>4} {'E[V]':>9} {eb:14.10f} {ea:14.10f} {ek:14.10f}"
              f" {dEa:9.1e} {dEk:9.1e}")
        if wv:
            dVa, dVk = abs(va - vb), abs(vk - vb)
            worst = max(worst, dVa, dVk)
            print(f"{ell:>4} {'Var[V]':>9} {vb:14.10f} {va:14.10f} {vk:14.10f}"
                  f" {dVa:9.1e} {dVk:9.1e}")
        print(f"{'':>4} ({time.time()-t0:.0f}s)")

    print("=" * 74)
    verdict = "PASS" if worst <= TOL else "FAIL"
    print(f"T6 worst residual {worst:.2e}   tol {TOL:.0e}   -> {verdict}")
    return worst <= TOL


if __name__ == "__main__":
    ok = main()
    exit(0 if ok else 1)
