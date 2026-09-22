#!/usr/bin/env python3
"""
Outcome stage 3 — E2 (primary), E3, E6, E1, secondary(analytic), validation.

Inputs: data/tm_PV_hat.npz (P̂(V), stage 2), tm_candidates_v11.csv (exact
E[V],Var(V)), exact P(V) for the 13 oracle members (recomputed here).
Everything is R/R0 units. τ_R = 1e-2. Registered η grid (§4.3).
"""
import numpy as np, csv, json, time
from collections import defaultdict
from tm_core import N, T, build_family, block_of, MU_A
from tm_T8 import grid_eta, w_table, pv_perm, pv_mix
from closure_maxent import closure1, closure2

TAU = 1e-2
RES = (1, 2, 4, 8, 16)
M = 1_000_000
ORACLE = ({f"Perm({s})" for s in ("identity", "shift-1", "shift-4",
           "shift-16", "shift-64", "random-0", "random-1", "random-2")}
          | {f"Mix({c})" for c in (0.02, 0.05, 0.1, 0.2, 0.5)})

fam = dict(build_family())
names = [n for n, _ in build_family()]

# exact moments
mom = defaultdict(dict)
for r in csv.DictReader(open("tm_candidates_v11.csv")):
    mom[(r["member"], int(r["resolution"]))][r["predictor"]] = float(r["value"])

PVhat = dict(np.load("data/tm_PV_hat.npz"))


def curve(pv, nb, etas):
    r = np.ceil(etas * nb).astype(int)
    return np.array([pv @ w_table(nb, ri) for ri in r])


def est_curve_se(pv, nb, etas):
    r = np.ceil(etas * nb).astype(int)
    R = np.empty(len(r)); SE = np.empty(len(r))
    for i, ri in enumerate(r):
        w = w_table(nb, ri)
        R[i] = pv @ w
        var = pv @ (w * w) - R[i] ** 2
        SE[i] = np.sqrt(max(var, 0) / M)
    return R, SE


def exact_pv(name, ell):
    if name.startswith("Perm"):
        return pv_perm(fam[name], ell)[0]
    if name.startswith("Mix"):
        return pv_mix(float(name[4:-1]), ell)
    return None


def main():
    rows = []          # D-table
    mj = {}
    t0 = time.time()
    for name in names:
        m_exceed = 0
        for ell in RES:
            nb = N // ell
            etas = grid_eta(nb)
            mean = mom[(name, ell)]["N_blk"]; var = mom[(name, ell)]["S_blk"]
            # closures (exact moments)
            p1 = closure1(nb, mean)
            p2, st2 = closure2(nb, mean, var)
            R1 = curve(p1, nb, etas)
            R2 = curve(p2, nb, etas) if p2 is not None else R1
            # truth
            pv_ex = exact_pv(name, ell) if name in ORACLE else None
            if pv_ex is not None:
                Rstar = curve(pv_ex, nb, etas); SEstar = np.zeros_like(Rstar)
                src = "exact"
            else:
                Rstar, SEstar = est_curve_se(PVhat[f"{name}|{ell}"], nb, etas)
                src = "est"
            D1 = float(np.max(np.abs(R1 - Rstar)))
            D2 = float(np.max(np.abs(R2 - Rstar)))
            exceed = D1 > TAU
            m_exceed += int(exceed)
            cls = ""
            if exceed:
                cls = "2nd" if D2 < D1 / 2 else "higher"
            rows.append(dict(member=name, ell=ell, src=src,
                             EV=round(mean, 5), Var=round(var, 5),
                             D1=D1, D2=D2, exceed=int(exceed), cls=cls))
        mj[name] = m_exceed

    # family verdict (existential)
    any_hdist = any(v >= 2 for v in mj.values())
    any_incon = any(v == 1 for v in mj.values())
    all_zero = all(v == 0 for v in mj.values())
    verdict = ("H0" if all_zero else "H_dist" if any_hdist else "inconclusive")

    # write D-table
    with open("data/tm_E2_table.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["member", "ell", "src", "EV", "Var",
                                          "D1", "D2", "exceed", "cls"])
        w.writeheader()
        for r in rows:
            w.writerow({**r, "D1": f"{r['D1']:.3e}", "D2": f"{r['D2']:.3e}"})

    print("=" * 74)
    print(f"E2 primary   τ_R={TAU}   family verdict: {verdict}")
    print("=" * 74)
    print("per-member m_j (# resolutions with D1>τ_R):")
    for name in names:
        tag = ""
        cells = [r for r in rows if r["member"] == name]
        exc = [(r["ell"], r["cls"]) for r in cells if r["exceed"]]
        if exc:
            tag = "  " + " ".join(f"ℓ{e}:{c}" for e, c in exc)
        print(f"  {name:16s} m_j={mj[name]}{tag}")

    # ---- E3 strong control: C6 pairs, sup curve distance between members ----
    print("\n" + "=" * 74)
    print("E3 (C6 strong control): curve distance within equal-coherence pairs")
    print("=" * 74)

    def rstar_curve(name, ell, etas, nb):
        pv_ex = exact_pv(name, ell) if name in ORACLE else None
        pv = pv_ex if pv_ex is not None else PVhat[f"{name}|{ell}"]
        return curve(pv, nb, etas)
    c6 = [("Perm(identity)", "Perm(shift-16)", 1),
          ("Perm(identity)", "Perm(random-0)", 1),
          ("Perm(shift-16)", "Perm(random-0)", 1),
          ("Mix(0.05)", "ShiftMix(0.05)", 1),
          ("Mix(0.2)", "ShiftMix(0.2)", 1)]
    for a, b, ell in c6:
        nb = N // ell; etas = grid_eta(nb)
        Ra = rstar_curve(a, ell, etas, nb); Rb = rstar_curve(b, ell, etas, nb)
        d = float(np.max(np.abs(Ra - Rb)))
        print(f"  ℓ{ell}  {a:16s} vs {b:16s}  supΔcurve = {d:.4f}"
              f"   (E[V]: {mom[(a,ell)]['N_blk']:.2f} vs {mom[(b,ell)]['N_blk']:.2f})")

    # ---- E6 analytic: layer-independent removal -> (1-η)^(T+1) for all ----
    print("\n" + "=" * 74)
    print("E6 (analytic): per-(layer,block) independent removal")
    print("=" * 74)
    print("  Under independent per-(layer,block) removal at η, dead-end survival")
    print("  = Π_{t=0}^{T} P[current block present] = (1-η)^(T+1) for EVERY member")
    print("  (kernel-independent; visited-block identity plays no role). Perm and")
    print("  all members coincide → η_½ = 1-2^(-1/(T+1)) = %.5f."
          % (1 - 2 ** (-1.0 / (T + 1))))

    # ---- E1 descriptive: family spread of R* at η_½ reference points ----
    print("\n" + "=" * 74)
    print("E1 (descriptive, no decision weight): family spread of R*/R0")
    print("=" * 74)
    for ell in (16, 4):
        nb = N // ell; etas = grid_eta(nb)
        Rs = np.array([rstar_curve(n, ell, etas, nb) for n in names])
        # at the grid point nearest η=0.2
        j = int(np.argmin(np.abs(etas - 0.2)))
        col = Rs[:, j]
        print(f"  ℓ{ell} η≈{etas[j]:.3f}: R*/R0 range [{col.min():.3f},{col.max():.3f}]"
              f" median {np.median(col):.3f}  IQR "
              f"[{np.percentile(col,25):.3f},{np.percentile(col,75):.3f}]")

    summary = dict(verdict=verdict, mj=mj,
                   n_hdist=sum(v >= 2 for v in mj.values()),
                   n_incon=sum(v == 1 for v in mj.values()),
                   n_h0=sum(v == 0 for v in mj.values()))
    json.dump(summary, open("data/tm_E2_summary.json", "w"), indent=1)
    print(f"\n({time.time()-t0:.0f}s)  saved data/tm_E2_table.csv, tm_E2_summary.json")


if __name__ == "__main__":
    main()
