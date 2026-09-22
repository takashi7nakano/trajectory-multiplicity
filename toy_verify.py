#!/usr/bin/env python3
"""
Analytic protocol-failure assessment — identity verification on a TOY chain.

NOT a registered family member and NOT a lesion run: n=8, T=4, exhaustive path
enumeration, exhaustive enumeration of block-removal sets. No registered seed is
touched. Purpose is only to confirm the algebra of the four identities used in
the assessment.
"""

import itertools
import numpy as np

n, T, ell = 8, 4, 2
nb = n // ell
lam = 1.0 - 0.5 ** (1.0 / T)
R0 = (1.0 - lam) ** T
mu = np.full(n, 1.0 / n)
b_of = np.arange(n) // ell


def circ(row):
    return np.array([[row[(y - x) % n] for y in range(n)] for x in range(n)])


# finite support (3 positions) -> dead rows possible
P_fin = circ(np.array([0.5, 0.3, 0, 0, 0, 0, 0, 0.2]))
# full support -> no dead rows ever
P_full = circ(np.array([0.3, 0.2, 0.1, 0.05, 0.05, 0.1, 0.1, 0.1]))

for P in (P_fin, P_full):
    assert abs(P.sum(1) - 1).max() < 1e-14 and abs(P.sum(0) - 1).max() < 1e-14


# ------------------------------------------------------------ brute force

def brute_reach(P, D, mode):
    """Exhaustive path sum for the lesioned chain."""
    S = [x for x in range(n) if x not in D]
    if mode == "redirect":
        Pl = P.copy()
        Pl[:, list(D)] = 0.0
        rs = Pl.sum(1)
        for x in range(n):
            if rs[x] > 0:
                Pl[x] /= rs[x]
    else:
        Pl = P.copy()
        Pl[:, list(D)] = 0.0
    tot = 0.0
    for path in itertools.product(range(n), repeat=T + 1):
        if path[0] in D:
            continue
        w = mu[path[0]]
        for t in range(T):
            w *= (1.0 - lam) * Pl[path[t], path[t + 1]]
            if w == 0.0:
                break
        tot += w
    return tot


def op_reach(P, D, mode):
    """Operator form: R = R0 * (1/n) 1_S^T A^T 1_S."""
    S = np.array([x not in D for x in range(n)])
    A = P.copy()
    A[:, ~S] = 0.0
    if mode == "redirect":
        rs = A.sum(1)
        alive = rs > 0
        A[alive] /= rs[alive][:, None]
    f = S.astype(float)
    for _ in range(T):
        f = A @ f
    return R0 * float(mu[S] @ f[S])


def dead_rows(P, D):
    A = P.copy(); A[:, list(D)] = 0.0
    return [x for x in range(n) if x not in D and A[x].sum() == 0.0]


def visited_count(P, path):
    v = set()
    for x in path:
        v.add(b_of[x])
    return len(v)


print("=" * 72)
print("identity 1/2 : lesioned reach, operator form vs exhaustive path sum")
print("identity 3   : redirect with no dead row  ->  R/R0 = |S|/n exactly")
print("=" * 72)
for label, P in (("finite-support", P_fin), ("full-support", P_full)):
    worst_r = worst_d = 0.0
    worst_deg = 0.0
    nZ = 0
    for r in range(1, nb):
        for Dblocks in itertools.combinations(range(nb), r):
            D = {x for x in range(n) if b_of[x] in Dblocks}
            for mode, w in (("redirect", "r"), ("dead-end", "d")):
                e = abs(brute_reach(P, D, mode) - op_reach(P, D, mode))
                if w == "r":
                    worst_r = max(worst_r, e)
                else:
                    worst_d = max(worst_d, e)
            Z = dead_rows(P, D)
            if Z:
                nZ += 1
            else:
                got = op_reach(P, D, "redirect") / R0
                worst_deg = max(worst_deg, abs(got - (n - len(D)) / n))
    print(f"  {label:16s} redirect {worst_r:.2e}  dead-end {worst_d:.2e}"
          f"  | no-dead-row cases: R/R0 = |S|/n to {worst_deg:.2e}"
          f"  ({nZ} of {2**nb-2} removal sets had dead rows)")

print()
print("=" * 72)
print("identity 4 : E_D[R_dead]/R0 = E_omega[ C(nb-V,r) / C(nb,r) ]")
print("             V = number of distinct blocks visited by the path")
print("=" * 72)
from math import comb
for label, P in (("finite-support", P_fin), ("full-support", P_full)):
    for r in range(1, nb):
        sets = list(itertools.combinations(range(nb), r))
        lhs = np.mean([op_reach(P, {x for x in range(n) if b_of[x] in Db},
                                "dead-end") for Db in sets]) / R0
        rhs = 0.0
        for path in itertools.product(range(n), repeat=T + 1):
            w = mu[path[0]]
            for t in range(T):
                w *= P[path[t], path[t + 1]]
                if w == 0.0:
                    break
            if w == 0.0:
                continue
            V = visited_count(P, path)
            rhs += w * (comb(nb - V, r) / comb(nb, r) if nb - V >= r else 0.0)
        print(f"  {label:16s} r={r}  LHS {lhs:.15f}  RHS {rhs:.15f}"
              f"  diff {abs(lhs-rhs):.2e}")
