#!/usr/bin/env python3
"""
PROTOCOL v1.1 §7 T8 — oracle curve validation.

For the 13 oracle members (8 permutations by direct enumeration, 5 Mix by the
occupancy closed form), build the EXACT robustness curve from the exact P(V),

    R_exact(η)/R0 = Σ_v P(v) w(v, r(η)),   w(v,r) = C(n_b−v, r)/C(n_b, r),

and compare it to the registered D-integrated trajectory estimator

    R̂(η)/R0 = (1/M) Σ_i w(V_i, r(η)),   M = 1e6,

on the registered η grid (§4.3), member trajectory seed 2026400000 + 1e4·index.
Agreement criterion: |R̂ − R_exact| ≤ 4 · SE(η) at every grid point
(pointwise implementation tolerance, §7.1 note).

This is the registered T8 validation. It computes curves only for the 13 oracle
members (truth exact there); it does not touch E2 or the 14 estimated members'
readout, and forms no removal set.
"""

import numpy as np
import time
from math import log2
from scipy.special import gammaln
from tm_core import N, T, block_of, build_family

M = 1_000_000
SEED_TRAJ_BASE = 2026400000
RESOLUTIONS = (1, 2, 4, 8, 16)
TOL_SE = 4.0
ABS_FLOOR = 1e-9      # SE→0 degenerate cells (deterministic V): the estimator
                      # has no sampling variability, so the meaningful check is
                      # machine-precision absolute agreement, not an SE ratio.
CHUNK = 100_000


def grid_eta(nb):
    """Registered η grid (§4.3): u_k = ln2·2^(-k/4), η_k = 1−exp(−u_k)."""
    K = int(round(4 * log2(nb)))
    u = np.log(2.0) * 2.0 ** (-np.arange(K + 1) / 4.0)
    return 1.0 - np.exp(-u)


def w_table(nb, r):
    """w(v,r) = C(n_b−v,r)/C(n_b,r) for v=1..nb, via log-gamma. 0 if n_b−v<r."""
    v = np.arange(1, nb + 1)
    top = nb - v
    logC = np.full(nb, -np.inf)
    ok = top >= r
    # log C(top,r) − log C(nb,r)
    lc = (gammaln(top[ok] + 1) - gammaln(r + 1) - gammaln(top[ok] - r + 1)
          - (gammaln(nb + 1) - gammaln(r + 1) - gammaln(nb - r + 1)))
    logC[ok] = lc
    w = np.zeros(nb)
    w[ok] = np.exp(logC[ok])
    return w                       # index 0 ↔ v=1


# ---------------------------------------------------------- exact P(V)
def pv_perm(P, ell):
    """Exact P(V) for a deterministic permutation: follow each of the N starts."""
    nb = N // ell
    b = block_of(ell)
    nxt = np.argmax(P, axis=1)
    Vs = np.empty(N, dtype=int)
    for x0 in range(N):
        x = x0
        seen = {b[x]}
        for _ in range(T):
            x = nxt[x]
            seen.add(b[x])
        Vs[x0] = len(seen)
    pv = np.bincount(Vs, minlength=nb + 1)[1:nb + 1].astype(float) / N
    return pv, Vs                  # pv index 0 ↔ v=1


def pv_mix(c, ell):
    """Exact P(V) for Mix(c): jumps K~Binom(T,c), K+1 iid uniform landings;
    V = distinct blocks among them (occupancy of m=K+1 balls in n_b cells).

    Occupancy distribution by the STABLE forward recursion (all-positive terms;
    the inclusion-exclusion / Stirling form cancels catastrophically for large m):
        q_m(v) = q_{m-1}(v)·v/n_b + q_{m-1}(v-1)·(n_b-v+1)/n_b,  q_0(0)=1.
    Then P(V) = Σ_k Binom(T,k;c) · q_{k+1}(·)."""
    from scipy.stats import binom
    nb = N // ell
    mmax = T + 1
    # q[m] over v=0..nb
    q = np.zeros((mmax + 1, nb + 1))
    q[0, 0] = 1.0
    v = np.arange(nb + 1)
    for m in range(1, mmax + 1):
        q[m, :] = q[m - 1, :] * (v / nb)
        q[m, 1:] += q[m - 1, :-1] * ((nb - v[1:] + 1) / nb)
    pv_full = np.zeros(nb + 1)
    for k in range(T + 1):
        wk = binom.pmf(k, T, c)
        if wk > 0.0:
            pv_full += wk * q[k + 1]
    return pv_full[1:nb + 1]        # index 0 ↔ v=1


# ---------------------------------------------------------- estimator samples
def sample_V_perm(Vs, seed):
    """Estimator V-samples for a permutation: M uniform starts, V = Vs[start]."""
    rng = np.random.Generator(np.random.PCG64(seed))
    out = np.empty(M, dtype=np.int32)
    for a in range(0, M, CHUNK):
        n = min(CHUNK, M - a)
        out[a:a + n] = Vs[rng.integers(0, N, size=n)]
    return out


def sample_V_mix(c, ell, seed):
    """Faithful Mix(c) V-sampler: K~Binom(T,c), K+1 uniform positions, count
    distinct blocks. Chunked; returns M visited-block counts."""
    rng = np.random.Generator(np.random.PCG64(seed))
    nb = N // ell
    out = np.empty(M, dtype=np.int32)
    for a in range(0, M, CHUNK):
        n = min(CHUNK, M - a)
        K = rng.binomial(T, c, size=n)
        mmax = int(K.max()) + 1
        pos = rng.integers(0, N, size=(n, mmax))
        blk = pos // ell
        # invalidate columns beyond K_i+1
        colidx = np.arange(mmax)[None, :]
        blk[colidx > K[:, None]] = -1
        blk.sort(axis=1)                     # -1s go first
        # distinct count = number of value-changes among non-(-1) entries
        diff = blk[:, 1:] != blk[:, :-1]
        distinct = diff.sum(axis=1) + 1      # +1 for first element of the row
        # but rows lead with some -1s; each row has ≥1 valid (start). The first
        # valid value is counted by the +1 only if row starts with a real block.
        # Correct: distinct valid = (#changes over full row) − (1 if any -1) ...
        # simpler: count distinct excluding sentinel.
        has_sentinel = (blk[:, 0] == -1)
        # changes include the -1→firstvalid transition; subtract it when present
        distinct = distinct - has_sentinel.astype(int)
        out[a:a + n] = distinct
    return out


def curve_from_pv(pv, nb, etas):
    r = np.ceil(etas * nb).astype(int)
    return np.array([pv @ w_table(nb, ri) for ri in r])


def curve_from_samples(Vs, nb, etas):
    r = np.ceil(etas * nb).astype(int)
    Rhat = np.empty(len(r)); SE = np.empty(len(r))
    for i, ri in enumerate(r):
        w = w_table(nb, ri)
        wi = w[Vs - 1]
        Rhat[i] = wi.mean()
        SE[i] = wi.std(ddof=1) / np.sqrt(len(wi))
    return Rhat, SE


def main():
    fam = dict((n, P) for n, P in build_family())
    order = [(6, "Perm(identity)"), (7, "Perm(shift-1)"), (8, "Perm(shift-4)"),
             (9, "Perm(shift-16)"), (10, "Perm(shift-64)"),
             (11, "Perm(random-0)"), (12, "Perm(random-1)"), (13, "Perm(random-2)"),
             (14, "Mix(0.02)"), (15, "Mix(0.05)"), (16, "Mix(0.1)"),
             (17, "Mix(0.2)"), (18, "Mix(0.5)")]

    print("=" * 78)
    print(f"T8 oracle curve validation   M={M}   tol {TOL_SE:.0f} pointwise SE")
    print("=" * 78)
    print(f"{'member':16s} {'ell':>3} {'maxΔ':>10} {'max|Δ|/SE':>10}"
          f" {'grid':>4} {'worstη':>8}")

    worst_ratio = 0.0
    fail = 0
    t0 = time.time()
    for idx, name in order:
        P = fam[name]
        seed = SEED_TRAJ_BASE + 10_000 * idx
        is_perm = name.startswith("Perm")
        if is_perm:
            per_ell_Vs = {}
        for ell in RESOLUTIONS:
            nb = N // ell
            etas = grid_eta(nb)
            if is_perm:
                pv, Vs_all = pv_perm(P, ell)
                Vsamp = sample_V_perm(Vs_all, seed + ell)
            else:
                c = float(name[4:-1])
                pv = pv_mix(c, ell)
                Vsamp = sample_V_mix(c, ell, seed + ell)
            R_exact = curve_from_pv(pv, nb, etas)
            Rhat, SE = curve_from_samples(Vsamp, nb, etas)
            d = np.abs(Rhat - R_exact)
            # a point passes if within 4 SE, OR (when SE≈0) within the absolute
            # machine-precision floor. ratio only meaningful where d > ABS_FLOOR.
            meaningful = d > ABS_FLOOR
            ratio = np.zeros_like(d)
            ratio[meaningful] = d[meaningful] / np.where(SE[meaningful] > 0,
                                                         SE[meaningful], np.inf)
            imax = int(np.argmax(ratio))
            worst_ratio = max(worst_ratio, float(ratio.max()))
            cell_fail = int(np.sum(ratio > TOL_SE))
            fail += cell_fail
            tag = "" if cell_fail == 0 else f"  <-- {cell_fail} FAIL"
            print(f"{name:16s} {ell:>3} {d.max():10.2e} {ratio.max():10.2f}"
                  f" {len(etas):>4} {etas[imax]:8.4f}{tag}")
    print("=" * 78)
    verdict = "PASS" if fail == 0 else f"FAIL ({fail} cells > {TOL_SE} SE)"
    print(f"T8 worst |Δ|/SE = {worst_ratio:.2f}   ({time.time()-t0:.0f}s)  -> {verdict}")
    return fail == 0


if __name__ == "__main__":
    ok = main()
    exit(0 if ok else 1)
