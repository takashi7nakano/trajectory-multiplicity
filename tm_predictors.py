#!/usr/bin/env python3
"""
TRAJECTORY_MULTIPLICITY_PROTOCOL v1.0 — predictors P1-P8 (§5)

Replaces the Haiku 4.5 draft in which P4, P6, P7 and P8 were constants
computed from the wrong objects. Each predictor here is computed from the
object §5 names, and every value reports the method used to obtain it.
"""

import numpy as np
import networkx as nx

from tm_core import (N, T, LAMBDA, R0, MU_A, block_of, lump_matrix,
                     block_joint, is_lumpable, _h, DETERMINISTIC,
                     committor_backward, reach, SEED_REMOVAL_BASE)

EPS_OCC = 0.1      # §5 P2
EPS_FLUX = 1e-4    # §5 P4-int
MC_SAMPLES = 100_000


# ------------------------------------------------------------------ Q1: P1-P3

def p1_hill(P, ell):
    """P1: layerwise Hill numbers of mu^R_t on blocks, alpha in {0,1,2}. §5

    Identical across the uniform family (= n/ell) by construction; reported
    to show the matching (§5).
    """
    S = lump_matrix(ell)
    mu = MU_A.copy()
    n0, n1, n2 = [], [], []
    for _ in range(T + 1):
        mb = S @ mu
        n0.append(float(np.count_nonzero(mb > 0)))
        n1.append(float(np.exp(_h(mb))))
        n2.append(float(1.0 / np.sum(mb ** 2)))
        mu = mu @ P
    return {"P1(a=0)": float(np.mean(n0)),
            "P1(a=1)": float(np.mean(n1)),
            "P1(a=2)": float(np.mean(n2))}


def p2_occupancy(P, ell):
    """P2: per-layer occupancy width, blocks with mass > eps/(n/ell); min
    over layers (primary statistic) and mean. §5, §6 Q1 primary."""
    S = lump_matrix(ell)
    nb = N // ell
    thr = EPS_OCC / nb
    mu = MU_A.copy()
    w = []
    for _ in range(T + 1):
        w.append(float(np.count_nonzero((S @ mu) > thr)))
        mu = mu @ P
    return {"P2": float(np.min(w)), "P2_mean": float(np.mean(w))}


def p3_committor_margin(P, ell):
    """P3: committor-margin statistic, min over reactive paths of q.

    In the primary embedding q_t(x) = (1-lam)^(T-t) is flat within layers, so
    the minimum over a reactive path is attained at t = 0: min q = R0. Computed
    from the numerical committor rather than asserted.
    """
    q = committor_backward(P)
    return {"P3": float(q[0].min())}


# ------------------------------------------------------- Q2 structural: P4, P5

def p4_maxflow(P, ell, solver="preflow_push"):
    """P4 primary: fractional max-flow of the reactive current. §5

    Layered block graph: node (t,b) with capacity mu^R_t(b) = ell/N, edge
    (t,b)->(t+1,b') with capacity F^R(b,b') = sum_{x in b, y in b'} mu(x)P(x,y).
    Value normalized by total reactive mass (= 1). Node capacities are realized
    by splitting each node into in/out.

    Cost grows as T*(N/ell)^2 edges, so this is run at the coarse resolutions
    and the analytic value is used as the cross-check (see report).
    """
    nb = N // ell
    F = block_joint(P, ell)
    mu_b = ell / N

    G = nx.DiGraph()
    BIG = 10.0
    for b in range(nb):
        G.add_edge("S", (0, b, "i"), capacity=BIG)
        G.add_edge((T, b, "o"), "K", capacity=BIG)
    for t in range(T + 1):
        for b in range(nb):
            G.add_edge((t, b, "i"), (t, b, "o"), capacity=mu_b)
    nz = np.argwhere(F > 0.0)
    for t in range(T):
        for b, bp in nz:
            G.add_edge((t, int(b), "o"), (t + 1, int(bp), "i"),
                       capacity=float(F[b, bp]))
    val, _ = nx.maximum_flow(G, "S", "K", flow_func=getattr(nx.algorithms.flow,
                                                            solver))
    return {"P4": float(val)}


def p4_int_routes(P, ell):
    """P4-int secondary: number of node-disjoint A->B block routes carrying
    flux > eps_f of total. §5 (reported, not used in E2).

    On the layered graph with unit node capacities the number of node-disjoint
    routes is the min over layers of the number of blocks that carry flux, and
    the number of usable transitions bounds the routes that can be extended.
    """
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import maximum_bipartite_matching
    F = block_joint(P, ell)
    thr = EPS_FLUX * F.sum()
    A = csr_matrix(F > thr)
    m = maximum_bipartite_matching(A, perm_type="column")
    return {"P4-int": float(np.count_nonzero(m >= 0))}


def p5_edge_network(P, ell):
    """P5: adjacent-layer reactive edge network statistics on blocks:
    mean out-degree, and the persistence fraction (flux staying in the same
    block index across a layer). §5"""
    F = block_joint(P, ell)
    thr = EPS_FLUX * F.sum()
    deg = np.count_nonzero(F > thr, axis=1)
    persistence = float(np.trace(F) / F.sum())
    return {"P5_outdeg": float(deg.mean()), "P5_persist": persistence}


# ----------------------------------------------------------- Q2 scalar: P6, P7

def p6_coherence(P, ell):
    """P6: scalar coherence C = sum_t I(X_t; X_{t+1} | S) on blocks. §5

    Time-homogeneous kernel with uniform block marginals, so
    C = T * I(X_0; X_1) with I computed from the explicitly formed block joint.
    Check: for a permutation at ell = 1 this must give T ln N (§3 C6).
    """
    J = block_joint(P, ell)
    mb = J.sum(axis=1)
    mbp = J.sum(axis=0)
    I = float(_h(mb) + _h(mbp) - _h(J))
    return {"P6": T * I}


def _jump_chain(Q):
    """Jump (self-loop-removed) kernel of a block chain Q, and its diagonal."""
    s = np.diag(Q).copy()
    Qt = Q - np.diag(s)
    rows = Qt.sum(axis=1)
    safe = rows > 0
    Qt[safe] /= rows[safe][:, None]
    return Qt, s


def p7_regime(P, ell):
    """Which exact regime P7 falls into, without running any Monte Carlo."""
    if DETERMINISTIC(P):
        return "enumerate"
    lumpable, _ = is_lumpable(P, ell)
    if lumpable:
        Q = block_joint(P, ell) * (N // ell)
        if np.ptp(np.diag(Q)) <= 1e-12:
            return "closed"
    return "mc"


def _is_circulant(P, tol=1e-14):
    idx = (np.arange(N)[None, :] - np.arange(N)[:, None]) % N
    return bool(np.abs(P - P[0][idx]).max() <= tol)


def p7_itinerary(P, ell, seed=None):
    """P7: itinerary N_eff^(1), N_eff^(2) on blocks, consecutive repeats
    removed. §5, §6 Q2 scalar primary (N_eff^(1)).

    Three regimes, and the regime used is reported:

    'enumerate'  deterministic kernel: each of the N start positions gives one
                 deterministic block path, so the itinerary law has <= N atoms
                 and is obtained exactly by enumeration.
    'closed'     block process strongly lumpable with a constant self-loop
                 probability s (all registered kernels are circulant or
                 block-symmetric, so the lumped chain has constant diagonal).
                 Then the number of jumps K ~ Binomial(T, 1-s) is independent
                 of the jump path, the law factorizes, and
                     H = H(K) + ln nb + E[K] h~,
                     sum p^2 = (1/nb) sum_k P(K=k)^2 r^k
                 with h~, r the row entropy and row collision of the jump
                 kernel. Exact.
    'mc'         neither: the block process is a non-lumpable function of a
                 Markov chain, whose sequence entropy has no closed form.
                 Monte Carlo with a registered seed; sum p^2 by the unbiased
                 collision estimator, H by plug-in (downward biased, so
                 N_eff^(1) here is a LOWER bound). Reported with its SE.

    §5 asserts P7 is 'exact'; regime 'mc' shows that claim does not hold for
    the whole registered family. Reported to B.
    """
    nb = N // ell
    b_of = block_of(ell)

    if DETERMINISTIC(P):
        nxt = np.argmax(P, axis=1)
        counts = {}
        for x0 in range(N):
            x = x0
            itin = [b_of[x]]
            for _ in range(T):
                x = nxt[x]
                if b_of[x] != itin[-1]:
                    itin.append(b_of[x])
            k = tuple(itin)
            counts[k] = counts.get(k, 0) + 1
        p = np.array(list(counts.values()), dtype=float) / N
        return {"P7(1)": float(np.exp(_h(p))),
                "P7(2)": float(1.0 / np.sum(p ** 2)),
                "P7_method": "enumerate", "P7_se": 0.0}

    lumpable, _ = is_lumpable(P, ell)
    Q = block_joint(P, ell) * nb        # row-normalized block kernel
    if lumpable and np.ptp(np.diag(Q)) <= 1e-12:
        Qt, s = _jump_chain(Q)
        s0 = float(s[0])
        k = np.arange(T + 1)
        from scipy.stats import binom
        pk = binom.pmf(k, T, 1.0 - s0)
        h_jump = float(_h(Qt[0]))                     # constant across rows
        r = float(np.sum(Qt[0] ** 2))
        H = float(_h(pk) + np.log(nb) + np.sum(pk * k) * h_jump)
        sum_p2 = float(np.sum(pk ** 2 * r ** k) / nb)
        return {"P7(1)": float(np.exp(H)),
                "P7(2)": float(1.0 / sum_p2),
                "P7_method": "closed", "P7_se": 0.0}

    # Monte Carlo. All members reaching this branch are circulant (Loc(sigma)
    # at ell >= 2), so a step is sampled once from the common row law rather
    # than per-trajectory.
    rng = np.random.Generator(np.random.PCG64(seed if seed is not None
                                              else SEED_REMOVAL_BASE))
    M = MC_SAMPLES
    x = rng.integers(0, N, size=M)
    cur = b_of[x].copy()
    keys = [cur.copy()]
    if _is_circulant(P):
        row0 = P[0]
        sup = np.flatnonzero(row0)
        pw = row0[sup]
        for _ in range(T):
            x = (x + rng.choice(sup, size=M, p=pw)) % N
            b = b_of[x]
            moved = b != cur
            cur = np.where(moved, b, cur)
            keys.append(np.where(moved, b, -1))
    else:
        cdf = np.cumsum(P, axis=1); cdf[:, -1] = 1.0
        for _ in range(T):
            u = rng.random(M)
            x = (cdf[x] < u[:, None]).sum(axis=1)
            b = b_of[x]
            moved = b != cur
            cur = np.where(moved, b, cur)
            keys.append(np.where(moved, b, -1))

    K = np.stack(keys, axis=1)
    uniq = {}
    for row in K:
        k = tuple(row[row >= 0])
        uniq[k] = uniq.get(k, 0) + 1
    cnt = np.array(list(uniq.values()), dtype=float)
    coll = float(np.sum(cnt * (cnt - 1.0)) / (M * (M - 1.0)))   # unbiased
    H_plugin = float(_h(cnt / M))
    saturated = len(cnt) == M            # every draw distinct: MC exhausted
    return {"P7(1)": float(np.exp(H_plugin)),
            "P7(2)": float(1.0 / coll) if coll > 0 else float("inf"),
            "P7_method": "mc-saturated" if saturated else "mc",
            "P7_distinct": float(len(cnt)), "P7_samples": float(M)}


# ----------------------------------------------------------------------- P8

def lesion_block_redirect(P, removed_mask):
    """Row-renormalized kernel with all columns in `removed_mask` deleted.
    Rows with no surviving successor become dead ends (-> F). §4 redirect mode."""
    Pl = P.copy()
    Pl[:, removed_mask] = 0.0
    rs = Pl.sum(axis=1)
    alive = rs > 0
    Pl[alive] /= rs[alive][:, None]
    return Pl, alive


def p8_accessibility(P, ell):
    """P8: accessibility a(ell). For a single persistent block lesion in
    redirect mode, the fraction of the redirected mass that reaches B,
    averaged over blocks. §5

    Computed, not proxied: for each block b, lesion b at every layer, solve the
    committor of the lesioned chain, and take the reach of the mass that was
    redirected (i.e. of the surviving support) relative to the unlesioned reach.
    """
    nb = N // ell
    b_of = block_of(ell)
    vals = []
    for b in range(nb):
        rm = b_of == b
        Pl, alive = lesion_block_redirect(P, rm)
        q = committor_backward(Pl)[0]
        surv = ~rm
        # mass starting outside the lesion, carried by the redirected chain
        num = float(MU_A[surv] @ q[surv])
        den = float(MU_A[surv].sum() * R0)
        vals.append(num / den if den > 0 else 0.0)
    return {"P8": float(np.mean(vals)), "P8_min": float(np.min(vals))}


def all_predictors(P, ell, with_maxflow=True, seed=None):
    out = {}
    out.update(p1_hill(P, ell))
    out.update(p2_occupancy(P, ell))
    out.update(p3_committor_margin(P, ell))
    if with_maxflow:
        out.update(p4_maxflow(P, ell))
    out.update(p4_int_routes(P, ell))
    out.update(p5_edge_network(P, ell))
    out.update(p6_coherence(P, ell))
    out.update(p7_itinerary(P, ell, seed=seed))
    out.update(p8_accessibility(P, ell))
    return out
