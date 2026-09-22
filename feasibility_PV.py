#!/usr/bin/env python3
"""
REPAIR SPEC v1.1 item 4 — algorithmic feasibility of exact P(V).

NO lesion computation. No R(eta), no eta_half, no removal set is formed.
Everything here is a functional of the UNPERTURBED reactive process:
V = number of distinct blocks visited by a trajectory over layers 0..T.

Three routes are compared against each other and against the exact E[V] that
the block-avoidance recursion already gives:

  route A  enumeration        deterministic kernels: n starts, one orbit each
  route B  occupancy closed   Mix(c): jumps ~ Binom(T,c), landings iid uniform
  route C  path sampling      general; V has support {1..n_b}, so its law is
                              a histogram on a small finite set
"""

import numpy as np
from math import comb
from scipy.stats import binom

from tm_core import build_family, N, T, block_of, MU_A, DETERMINISTIC
from assess_analytic import stay_prob


def exact_EV_recursion(P, ell):
    """E[V] by single-block avoidance. Exact, no sampling. O(n_b) recursions."""
    b = block_of(ell)
    return sum(1.0 - stay_prob(P, b != k) for k in range(N // ell))


# ---------------------------------------------------------------- route A
def PV_enumerate(P, ell):
    """Deterministic kernel: each start position gives one orbit."""
    b = block_of(ell)
    nxt = np.argmax(P, axis=1)
    counts = np.zeros(N // ell + 1)
    for x0 in range(N):
        x, seen = x0, {b[x0]}
        for _ in range(T):
            x = nxt[x]
            seen.add(b[x])
        counts[len(seen)] += 1
    return counts / N


# ---------------------------------------------------------------- route B
def PV_occupancy_mix(c, ell):
    """Mix(c) = (1-c) I + c U.  A step jumps to a uniform position w.p. c and
    stays w.p. 1-c, so the visited set is {x_0} union {K landing points} with
    K ~ Binom(T, c) and landings iid uniform over blocks.
    V is then the occupancy count of m = K+1 iid uniform draws into n_b cells."""
    nb = N // ell
    pv = np.zeros(nb + 1)
    for k in range(T + 1):
        w = binom.pmf(k, T, c)
        if w == 0.0:
            continue
        m = k + 1
        for v in range(1, min(nb, m) + 1):
            s = sum((-1) ** j * comb(v, j) * ((v - j) / nb) ** m
                    for j in range(v + 1))
            pv[v] += w * comb(nb, v) * s
    return pv


# ---------------------------------------------------------------- route C
def PV_sample(P, ell, M=200_000, seed=0):
    """General path sampling. Circulant kernels sample the step once from the
    common row law; otherwise per-row inverse CDF."""
    rng = np.random.Generator(np.random.PCG64(seed))
    b = block_of(ell)
    nb = N // ell
    idx = (np.arange(N)[None, :] - np.arange(N)[:, None]) % N
    circ = np.abs(P - P[0][idx]).max() <= 1e-14

    x = rng.integers(0, N, size=M)
    seen = np.zeros((M, nb), dtype=bool)
    seen[np.arange(M), b[x]] = True
    if circ:
        sup = np.flatnonzero(P[0]); pw = P[0][sup]
        for _ in range(T):
            x = (x + rng.choice(sup, size=M, p=pw)) % N
            seen[np.arange(M), b[x]] = True
    else:
        cdf = np.cumsum(P, axis=1); cdf[:, -1] = 1.0
        for _ in range(T):
            x = (cdf[x] < rng.random(M)[:, None]).sum(axis=1)
            seen[np.arange(M), b[x]] = True
    V = seen.sum(axis=1)
    pv = np.bincount(V, minlength=nb + 1).astype(float) / M
    return pv, float(V.std(ddof=1) / np.sqrt(M))


def EV(pv):
    return float(np.dot(np.arange(len(pv)), pv))


if __name__ == "__main__":
    fam = dict(build_family())
    ell = 8
    nb = N // ell
    print(f"resolution ell = {ell}, n_b = {nb}, T = {T}")
    print(f"{'member':16s} {'exact E[V]':>12s} {'route':>10s} "
          f"{'route E[V]':>12s} {'|diff|':>10s} {'supp':>5s} {'SE':>8s}")

    for name in ("Perm(identity)", "Perm(shift-4)", "Perm(random-0)"):
        P = fam[name]
        ex = exact_EV_recursion(P, ell)
        pv = PV_enumerate(P, ell)
        print(f"{name:16s} {ex:12.6f} {'enumerate':>10s} {EV(pv):12.6f} "
              f"{abs(EV(pv)-ex):10.2e} {int((pv>0).sum()):5d} {0.0:8.1e}")

    for c in (0.02, 0.1, 0.5):
        P = fam[f"Mix({c})"]
        ex = exact_EV_recursion(P, ell)
        pv = PV_occupancy_mix(c, ell)
        print(f"{'Mix(%s)'%c:16s} {ex:12.6f} {'occupancy':>10s} {EV(pv):12.6f} "
              f"{abs(EV(pv)-ex):10.2e} {int((pv>1e-15).sum()):5d} {0.0:8.1e}")

    for name in ("Loc(1)", "Loc(8)", "Loc(32)", "Mix(0.1)", "Block(16,0.02)",
                 "ShiftMix(0.05)", "Perm(random-0)"):
        P = fam[name]
        ex = exact_EV_recursion(P, ell)
        pv, se = PV_sample(P, ell, seed=4242)
        print(f"{name:16s} {ex:12.6f} {'sample':>10s} {EV(pv):12.6f} "
              f"{abs(EV(pv)-ex):10.2e} {int((pv>0).sum()):5d} {se:8.1e}")
