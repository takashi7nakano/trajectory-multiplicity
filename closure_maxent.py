#!/usr/bin/env python3
"""
PROTOCOL v1.1 — maximum-entropy closures on the support v = 1..n_b.

B REVIEW of REPAIR SPEC v1.1, item 1: closure-1 and closure-2 must be defined
symmetrically, as maximum-entropy laws under moment constraints, so that the
only difference between the rungs of the hierarchy is how much information is
supplied and not how the closure is built.

    closure-1   p1(v) ∝ exp(l1 v)              mean matched
    closure-2   p2(v) ∝ exp(l1 v + l2 v^2)     mean and variance matched

NO lesion computation. Only the already-computed unperturbed moments E[V] and
Var(V) are used, to check that the two closures are well posed on the family.
"""

import numpy as np
from scipy.optimize import brentq, root


def _law(nb, lams):
    v = np.arange(1, nb + 1, dtype=float)
    z = lams[0] * v + (lams[1] * v * v if len(lams) > 1 else 0.0)
    z = z - z.max()
    p = np.exp(z)
    return v, p / p.sum()


def closure1(nb, mean, tol=1e-12):
    """Max-ent law on {1..n_b} with the mean matched. Unique for mean in (1,n_b);
    the endpoints are the degenerate point masses."""
    if mean <= 1.0 + tol:
        p = np.zeros(nb); p[0] = 1.0; return p
    if mean >= nb - tol:
        p = np.zeros(nb); p[-1] = 1.0; return p
    f = lambda l: _law(nb, [l])[0] @ _law(nb, [l])[1] - mean
    lo, hi = -50.0, 50.0
    l1 = brentq(f, lo, hi, xtol=1e-14, rtol=1e-15)
    return _law(nb, [l1])[1]


def var_max(nb, mean):
    """Largest variance attainable on {1..n_b} at the given mean: the two-point
    law on the endpoints. Var = (mean-1)(n_b-mean)."""
    return (mean - 1.0) * (nb - mean)


def closure2(nb, mean, var, tol=1e-10):
    """Max-ent law on {1..n_b} with mean and variance matched. Well posed when
    0 < var < var_max(n_b, mean); var -> 0 degenerates to a point mass."""
    if var <= tol:
        p = np.zeros(nb)
        k = int(round(mean))
        k = min(max(k, 1), nb)
        p[k - 1] = 1.0
        return p, "degenerate(var=0)"
    vmax = var_max(nb, mean)
    if var >= vmax * (1 - 1e-9):
        return None, f"outside moment space (var {var:.4f} >= max {vmax:.4f})"

    def res(l):
        v, p = _law(nb, l)
        m1 = v @ p
        m2 = (v * v) @ p
        return [m1 - mean, (m2 - m1 * m1) - var]

    for guess in ([0.0, 0.0], [0.1, -0.01], [-0.1, 0.01], [1.0, -0.1]):
        s = root(res, guess, method="hybr", tol=1e-13)
        if s.success:
            v, p = _law(nb, s.x)
            m1 = v @ p
            m2 = (v * v) @ p
            if abs(m1 - mean) < 1e-8 and abs(m2 - m1 * m1 - var) < 1e-6:
                return p, "ok"
    return None, "no convergence"


if __name__ == "__main__":
    import csv, collections
    rows = list(csv.DictReader(open("tm_candidates_v11.csv")))
    d = collections.defaultdict(dict)
    for r in rows:
        d[(r["member"], int(r["resolution"]))][r["predictor"]] = float(r["value"])

    print("well-posedness of the two closures on the registered family")
    print("(moments are the already-computed unperturbed E[V], Var(V))\n")
    bad = []
    for ell in (16, 8, 4, 2, 1):
        nb = 256 // ell
        n_ok = n_deg = n_bad = 0
        for (m, e), vals in d.items():
            if e != ell:
                continue
            mean, var = vals["N_blk"], vals["S_blk"]
            p1 = closure1(nb, mean)
            m1 = np.arange(1, nb + 1) @ p1
            ok1 = abs(m1 - mean) < 1e-7 or mean <= 1 + 1e-9 or mean >= nb - 1e-9
            p2, status = closure2(nb, mean, var)
            if status == "ok":
                n_ok += 1
            elif status.startswith("degenerate"):
                n_deg += 1
            else:
                n_bad += 1
                bad.append((m, ell, mean, var, var_max(nb, mean), status))
            if not ok1:
                bad.append((m, ell, mean, var, None, "closure1 failed"))
        print(f"  ell={ell:2d} n_b={nb:3d}   closure-2: {n_ok} solved, "
              f"{n_deg} degenerate(var=0), {n_bad} outside/failed")

    if bad:
        print("\n  cases needing a registered rule:")
        for m, e, mean, var, vmax, s in bad:
            print(f"    {m:18s} ell={e:2d} mean={mean:9.4f} var={var:9.4f} "
                  f"var_max={vmax if vmax is None else round(vmax,4)}  {s}")
    else:
        print("\n  no failures")
