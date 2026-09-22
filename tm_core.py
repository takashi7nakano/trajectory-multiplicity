#!/usr/bin/env python3
"""
TRAJECTORY_MULTIPLICITY_PROTOCOL v1.0 — Pipeline A core
Exact kernel construction, embedding, committor, Doob transform.

Reimplementation (Opus 5, 2026-09-21) replacing the Haiku 4.5 draft whose
T2-T5 and P4/P6/P7/P8 were placeholders. No check returns a hardcoded verdict.

Protocol references in docstrings are to PROTOCOL_v1.0.md (FROZEN 2026-09-21).
"""

import numpy as np

# ---------------------------------------------------------------- constants §2
N = 256                              # positions per layer
T = 64                               # layers
LAMBDA = 1.0 - 0.5 ** (1.0 / T)      # uniform leak, R0 = (1-LAMBDA)^T = 1/2
R0 = (1.0 - LAMBDA) ** T
RESOLUTIONS = (1, 2, 4, 8, 16)
SEEDS_PERM = (2026101, 2026102, 2026103)      # §7
SEED_REMOVAL_BASE = 2026200000                # §7
SEED_SPLIT_E7 = 2026300                       # §6 E7

MU_A = np.full(N, 1.0 / N)           # §2 step 1: uniform reactive marginal


# ------------------------------------------------------------ kernel menus §3
# Every construction below is *exactly* doubly stochastic by symmetry.
# No Sinkhorn balancing is applied: iterating on an already-exact matrix only
# injects floating point noise, and T1 verifies exactness to 1e-12.

def k_loc(sigma):
    """Loc(sigma): symmetric nearest-range walk, periodic. §3"""
    P = np.zeros((N, N))
    idx = np.arange(N)
    w = 1.0 / (2 * sigma + 1)
    for d in range(-sigma, sigma + 1):
        P[idx, (idx + d) % N] = w
    return P


def k_perm(pi):
    """Perm(pi): deterministic permutation matrix. §3"""
    P = np.zeros((N, N))
    P[np.arange(N), pi] = 1.0
    return P


def k_shift(s):
    return k_perm((np.arange(N) + s) % N)


def k_mix(c):
    """Mix(c) = (1-c) I + c U. §3"""
    return (1.0 - c) * np.eye(N) + c * np.full((N, N), 1.0 / N)


def k_block(M, c):
    """Block(M,c): within-block uniform w.p. 1-c, uniform jump out w.p. c. §3"""
    bs = N // M
    blk = np.arange(N) // bs
    same = blk[:, None] == blk[None, :]
    return np.where(same, (1.0 - c) / bs, c / (N - bs))


def k_shiftmix(s, c):
    """ShiftMix = (1-c) shift(s) + c U. §3"""
    return (1.0 - c) * k_shift(s) + c * np.full((N, N), 1.0 / N)


def build_family():
    """The registered family, §3. Returns list of (name, kernel).

    Menu count: Loc 6 + Perm 8 + Mix 5 + Block 6 + ShiftMix 2 = 27.
    §3 states 'the menus below give >= 30'; the menus as written give 27.
    27 still satisfies the registered target N_fam >= 24, so the study is
    unaffected, but the discrepancy is a protocol-text error reported to B.
    """
    fam = []
    for sigma in (1, 2, 4, 8, 16, 32):
        fam.append((f"Loc({sigma})", k_loc(sigma)))
    fam.append(("Perm(identity)", k_perm(np.arange(N))))
    for s in (1, 4, 16, 64):
        fam.append((f"Perm(shift-{s})", k_shift(s)))
    for i, seed in enumerate(SEEDS_PERM):
        pi = np.random.Generator(np.random.PCG64(seed)).permutation(N)
        fam.append((f"Perm(random-{i})", k_perm(pi)))
    for c in (0.02, 0.05, 0.1, 0.2, 0.5):
        fam.append((f"Mix({c})", k_mix(c)))
    for M in (4, 16, 64):
        for c in (0.02, 0.1):
            fam.append((f"Block({M},{c})", k_block(M, c)))
    for c in (0.05, 0.2):
        fam.append((f"ShiftMix({c})", k_shiftmix(16, c)))
    return fam


DETERMINISTIC = lambda P: np.all((P == 0.0) | (P == 1.0))


# ------------------------------------------------------- embedding/committor §2

def committor_backward(P, lam=LAMBDA, n_layers=T):
    """Backward committor recursion on the *embedded* chain.

    Embedded chain (§2 step 3): P_t(x,y) = (1-lam) P^R(x,y), P_t(x,F) = lam.
    q_T = 1 on layer T (absorbed into B), q(F) = 0, so for t < T

        q_t(x) = (1-lam) sum_y P^R(x,y) q_{t+1}(y).

    Returns q[t] for t = 0..n_layers, shape (n_layers+1, N). Nothing about the
    closed form (1-lam)^(T-t) is assumed; T2 compares against it.
    """
    q = np.empty((n_layers + 1, N))
    q[n_layers] = 1.0
    for t in range(n_layers - 1, -1, -1):
        q[t] = (1.0 - lam) * (P @ q[t + 1])
    return q


def doob_kernel(P, q, t, lam=LAMBDA):
    """Doob h-transform of the embedded chain at layer t. §2 step 3.

    P^Doob_t(x,y) = P_t(x,y) q_{t+1}(y) / q_t(x),  P_t(x,y) = (1-lam) P^R(x,y).
    Should reproduce P^R exactly. Uses the numerically computed q, not the
    closed form, so T2 is a real test of the embedding rather than a tautology.
    """
    return (1.0 - lam) * P * q[t + 1][None, :] / q[t][:, None]


def reach(P, lam=LAMBDA, n_layers=T, mu=MU_A):
    """Overall reactive probability R = mu_A . q_0."""
    return float(mu @ committor_backward(P, lam, n_layers)[0])


# ------------------------------------------------------------ block machinery

def block_index(ell):
    """Nested dyadic partition, origin 0. §2"""
    return np.arange(N) // ell


def block_of(ell, shift=0):
    """Partition with origin shifted by `shift` positions. §2 C1 sensitivity."""
    return ((np.arange(N) - shift) % N) // ell


def lump_matrix(ell, shift=0):
    """(n_blocks x N) 0/1 aggregation matrix S with S[b,x] = 1 iff x in b."""
    b = block_of(ell, shift)
    nb = N // ell
    S = np.zeros((nb, N))
    S[b, np.arange(N)] = 1.0
    return S


def block_joint(P, ell, shift=0, mu=None):
    """Block-level joint J(b,b') = sum_{x in b, y in b'} mu(x) P(x,y)."""
    if mu is None:
        mu = MU_A
    S = lump_matrix(ell, shift)
    return S @ (mu[:, None] * P) @ S.T


def is_lumpable(P, ell, shift=0, tol=1e-12):
    """Strong lumpability: sum_{y in b'} P(x,y) equal for all x in the same b.

    Required for the block process to be Markov, which is what makes the exact
    itinerary DP (P7) available.
    """
    S = lump_matrix(ell, shift)
    Pb = P @ S.T                     # (N x nb): mass from x into each block
    b = block_of(ell, shift)
    nb = N // ell
    worst = 0.0
    for blk in range(nb):
        rows = Pb[b == blk]
        if rows.shape[0] > 1:
            worst = max(worst, float(np.abs(rows - rows[0]).max()))
    return worst <= tol, worst


# ------------------------------------------------------------------- entropies

def _h(p, axis=None):
    """Shannon entropy in nats, 0 log 0 = 0."""
    p = np.asarray(p, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        term = np.where(p > 0.0, p * np.log(p), 0.0)
    return -term.sum(axis=axis)


def entropy_pieces(P, ell=1, n_layers=T):
    """Pieces of the T3 identity H(Gamma|S) = sum_t H_t - C, computed by
    three independent code paths on the resolution-`ell` block process.

    H_path : chain rule, H(X_0) + sum_t H(X_{t+1}|X_t)
    sum_Ht : forward-propagated marginals, each entropy computed separately
    C      : sum_t I(X_t;X_{t+1}) from the explicitly formed joints  (= P6)
    """
    S = lump_matrix(ell)
    mu = MU_A.copy()

    h_marg = []          # H(X_t) on blocks
    h_cond = []          # H(X_{t+1} | X_t) on blocks
    mi = []              # I(X_t; X_{t+1}) on blocks

    for t in range(n_layers + 1):
        mb = S @ mu
        h_marg.append(_h(mb))
        if t < n_layers:
            J = S @ (mu[:, None] * P) @ S.T          # joint on blocks
            mb_next = J.sum(axis=0)
            h_cond.append(_h(J) - _h(mb))            # H(X,Y) - H(X)
            mi.append(_h(mb) + _h(mb_next) - _h(J))
            mu = mu @ P

    H_path = h_marg[0] + float(np.sum(h_cond))
    sum_Ht = float(np.sum(h_marg))
    C = float(np.sum(mi))
    return dict(H_path=H_path, sum_Ht=sum_Ht, C=C,
                identity_residual=abs(H_path - (sum_Ht - C)),
                h_marg=np.array(h_marg), mi=np.array(mi))
