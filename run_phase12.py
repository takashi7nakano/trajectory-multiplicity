#!/usr/bin/env python3
"""Phase 1-2 runner: T1-T5 on the whole family, then all predictors at all
resolutions. Writes tm_checks_v11.json and tm_predictors_v11.csv."""

import json, time, csv, hashlib
import numpy as np

from tm_core import build_family, RESOLUTIONS, N, T, LAMBDA, R0
from tm_checks import T1, T2, T3, T4, T5
import tm_predictors as pr

MAXFLOW_RESOLUTIONS = (16, 8, 4)     # exact solve; coarser is O(T (N/ell)^2)


def cheap_predictors(P, ell):
    """Predictor set used inside T5 (C1 invariance): everything except the
    max-flow solve, which is verified separately, and with P7 in its exact
    regimes only (the MC regime carries sampling noise that is not a C1
    violation and is reported separately)."""
    out = {}
    out.update(pr.p1_hill(P, ell))
    out.update(pr.p2_occupancy(P, ell))
    out.update(pr.p3_committor_margin(P, ell))
    out.update(pr.p4_int_routes(P, ell))
    out.update(pr.p5_edge_network(P, ell))
    out.update(pr.p6_coherence(P, ell))
    out.update(pr.p8_accessibility(P, ell))
    if pr.p7_regime(P, ell) != "mc":
        d7 = pr.p7_itinerary(P, ell, seed=99)
        out["P7(1)"] = d7["P7(1)"]
        out["P7(2)"] = d7["P7(2)"]
    return out


def main():
    fam = build_family()
    print(f"family: {len(fam)} members   N={N} T={T} lambda={LAMBDA:.10f} R0={R0:.12f}")

    checks = {}
    print("\n=== T1-T5 ===")
    t_all = time.time()
    for name, P in fam:
        rec = {}
        ok1, d1 = T1(P); rec["T1"] = dict(passed=bool(ok1), **d1)
        ok2, d2 = T2(P); rec["T2"] = dict(passed=bool(ok2), **d2)
        ok3, d3 = T3(P); rec["T3"] = dict(passed=bool(ok3),
                                          residual=d3["residual"], tol=d3["tol"])
        ok4, d4 = T4(P, n_lesions=100, seed=7)
        rec["T4"] = dict(passed=bool(ok4), **d4)
        ok5, d5 = T5(P, cheap_predictors, n_perm=10, seed=11)
        rec["T5"] = dict(passed=bool(ok5), failed=d5["failed"], tol=d5["tol"],
                         worst=max(d5["worst_relative"].values()))
        rec["all_passed"] = all(rec[k]["passed"] for k in ("T1", "T2", "T3", "T4", "T5"))
        checks[name] = rec
        flag = "PASS" if rec["all_passed"] else "FAIL"
        print(f"  {name:18s} {flag}  T1 {d1['residual']:.1e}  T2 {d2['residual']:.1e}"
              f"  T3 {d3['residual']:.1e}  T4 {d4['residual']:.1e}"
              f"  T5 {rec['T5']['worst']:.1e}")
    print(f"  ({time.time()-t_all:.0f}s)")

    n_fail = sum(1 for r in checks.values() if not r["all_passed"])
    print(f"\nT1-T5: {len(fam)-n_fail}/{len(fam)} members pass")

    print("\n=== predictors ===")
    rows = []
    t_all = time.time()
    for mi, (name, P) in enumerate(fam):
        for li, ell in enumerate(RESOLUTIONS):
            seed = 2026200000 + 10000 * mi + 100 * li
            d = {}
            d.update(pr.p1_hill(P, ell))
            d.update(pr.p2_occupancy(P, ell))
            d.update(pr.p3_committor_margin(P, ell))
            if ell in MAXFLOW_RESOLUTIONS:
                d.update(pr.p4_maxflow(P, ell))
                d["P4_method"] = "solved"
            else:
                d["P4"] = float("nan")
                d["P4_method"] = "not-solved"
            d.update(pr.p4_int_routes(P, ell))
            d.update(pr.p5_edge_network(P, ell))
            d.update(pr.p6_coherence(P, ell))
            d.update(pr.p7_itinerary(P, ell, seed=seed))
            d.update(pr.p8_accessibility(P, ell))
            for k, v in d.items():
                rows.append(dict(member=name, member_index=mi, resolution=ell,
                                 predictor=k, value=v))
        print(f"  {name:18s} done ({time.time()-t_all:.0f}s)")

    with open("tm_predictors_v11.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["member", "member_index",
                                          "resolution", "predictor", "value"])
        w.writeheader(); w.writerows(rows)
    with open("tm_checks_v11.json", "w") as f:
        json.dump(checks, f, indent=1, default=float)

    for fn in ("tm_predictors_v11.csv", "tm_checks_v11.json"):
        h = hashlib.sha256(open(fn, "rb").read()).hexdigest()
        print(f"  {fn}  SHA-256 {h}")


if __name__ == "__main__":
    main()
