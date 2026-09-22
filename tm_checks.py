#!/usr/bin/env python3
"""
TRAJECTORY_MULTIPLICITY_PROTOCOL v1.0 — code checks T1-T5 (§7)

Every check computes a residual and compares it to the registered tolerance.
No check returns a hardcoded verdict, and no check compares an expression to
itself. Tolerances: T1, T2, T4 -> 1e-12; T3 -> 1e-10; T5 -> 1e-10.
"""

import numpy as np

from tm_core import (N, T, LAMBDA, R0, MU_A, committor_backward, doob_kernel,
                     entropy_pieces, block_of, RESOLUTIONS)

TOL_12 = 1e-12
TOL_10 = 1e-10


# ---------------------------------------------------------------------- T1

def T1(P):
    """Reactive marginals: forward propagation of mu_A under P^R reproduces
    uniform to 1e-12. Also reports the double-stochasticity residuals, which
    are what the propagation actually tests."""
    row = float(np.abs(P.sum(axis=1) - 1.0).max())
    col = float(np.abs(P.sum(axis=0) - 1.0).max())
    neg = float(max(0.0, -P.min()))
    mu = MU_A.copy()
    worst = 0.0
    for _ in range(T):
        mu = mu @ P
        worst = max(worst, float(np.abs(mu - 1.0 / N).max()))
    res = max(row, col, neg, worst)
    return res <= TOL_12, dict(row_sum=row, col_sum=col, negativity=neg,
                               marginal=worst, residual=res, tol=TOL_12)


# ---------------------------------------------------------------------- T2

def T2(P):
    """Embedding: (a) the committor of the embedded chain, obtained by backward
    recursion, equals (1-lam)^(T-t) to 1e-12; (b) the Doob transform of the
    embedded chain reproduces P^R to 1e-12; (c) R0 = 1/2."""
    q = committor_backward(P)
    t = np.arange(T + 1)
    closed = (1.0 - LAMBDA) ** (T - t)
    e_q = float(np.abs(q - closed[:, None]).max())

    e_doob = 0.0
    for layer in (0, T // 3, T // 2, T - 1):
        e_doob = max(e_doob, float(np.abs(doob_kernel(P, q, layer) - P).max()))

    e_r0 = abs(float(MU_A @ q[0]) - 0.5)
    res = max(e_q, e_doob, e_r0)
    return res <= TOL_12, dict(committor=e_q, doob=e_doob, R0=e_r0,
                               residual=res, tol=TOL_12)


# ---------------------------------------------------------------------- T3

def T3(P):
    """Entropy identity H(Gamma|S) = sum_t H_t - C to 1e-10, at every
    resolution. The three quantities are built by independent code paths
    (chain rule / propagated marginals / explicit joints) in tm_core."""
    worst = 0.0
    per = {}
    for ell in RESOLUTIONS:
        d = entropy_pieces(P, ell)
        per[ell] = d["identity_residual"]
        worst = max(worst, d["identity_residual"])
    return worst <= TOL_10, dict(per_resolution=per, residual=worst,
                                 tol=TOL_10)


# ---------------------------------------------------------------------- T4

def _reach_edge_lesion(P, t_les, x, y, mode):
    """Direct recomputation of R with a single edge (x->y) lesioned at layer
    t_les. redirect: row renormalized over survivors; dead-end: mass lost."""
    row = P[x].copy()
    row[y] = 0.0
    if mode == "redirect":
        s = row.sum()
        if s > 0:
            row = row / s
    q = np.ones(N)
    for t in range(T - 1, -1, -1):
        if t == t_les:
            qn = (1.0 - LAMBDA) * (P @ q)
            qn[x] = (1.0 - LAMBDA) * float(row @ q)
            q = qn
        else:
            q = (1.0 - LAMBDA) * (P @ q)
    return float(MU_A @ q)


def T4(P, n_lesions=100, seed=0):
    """Single-lesion exactness to 1e-12 for 100 random single-edge lesions per
    member. §7

    Feed-forward exactness (§2, rev.3 §4.3) gives closed forms, derived from the
    committor being flat within layers:

        redirect   dR = 0                     if P^R(x,y) < 1
                   dR = -R0 / N               if P^R(x,y) = 1  (row dies)
        dead-end   dR = -R0 * P^R(x,y) / N

    Both are compared against the direct backward recursion.
    """
    rng = np.random.Generator(np.random.PCG64(seed))
    nz = np.argwhere(P > 0.0)
    pick = rng.choice(len(nz), size=min(n_lesions, len(nz)), replace=False)
    layers = rng.integers(0, T, size=len(pick))

    worst_r = worst_d = 0.0
    for (idx, t_les) in zip(pick, layers):
        x, y = int(nz[idx][0]), int(nz[idx][1])
        pxy = float(P[x, y])

        pred_r = 0.0 if pxy < 1.0 else -R0 / N
        got_r = _reach_edge_lesion(P, t_les, x, y, "redirect") - R0
        worst_r = max(worst_r, abs(got_r - pred_r))

        pred_d = -R0 * pxy / N
        got_d = _reach_edge_lesion(P, t_les, x, y, "dead-end") - R0
        worst_d = max(worst_d, abs(got_d - pred_d))

    res = max(worst_r, worst_d)
    return res <= TOL_12, dict(redirect=worst_r, dead_end=worst_d,
                               n_lesions=len(pick), residual=res, tol=TOL_12)


# ---------------------------------------------------------------------- T5

def block_label_permutation(ell, rng):
    """Position permutation induced by a random relabelling of the blocks:
    block b -> sigma(b), within-block offset preserved."""
    nb = N // ell
    sigma = rng.permutation(nb)
    x = np.arange(N)
    b, r = x // ell, x % ell
    return sigma[b] * ell + r


def T5(P, predictor_fn, n_perm=10, seed=0, resolutions=RESOLUTIONS):
    """C1: every predictor identical under 10 random block-label
    permutations, to 1e-10. §7, §5

    A predictor that reads block indices as values rather than as labels fails
    here; per §8 K5 such a predictor is excluded and the study continues.
    """
    rng = np.random.Generator(np.random.PCG64(seed))
    worst = {}
    for ell in resolutions:
        base = predictor_fn(P, ell)
        for _ in range(n_perm):
            perm = block_label_permutation(ell, rng)
            Pp = P[np.ix_(perm, perm)]
            got = predictor_fn(Pp, ell)
            for k, v in base.items():
                if not isinstance(v, (int, float)) or isinstance(v, bool):
                    continue
                g = got.get(k)
                if not isinstance(g, (int, float)):
                    continue
                scale = max(1.0, abs(v))
                d = abs(g - v) / scale
                worst[k] = max(worst.get(k, 0.0), float(d))
    failed = {k: v for k, v in worst.items() if v > TOL_10}
    return len(failed) == 0, dict(worst_relative=worst, failed=failed,
                                  tol=TOL_10)
