#!/usr/bin/env python3
"""Figures for the paper, from frozen Pipeline A outputs. No new analysis:
categories and curves are read/recomputed from the frozen E2 table, moments,
and P(V). Okabe-Ito CVD-safe categorical palette; secondary glyph encoding."""
import numpy as np, csv
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from tm_core import N, T
from tm_T8 import grid_eta, w_table, pv_perm, pv_mix
from closure_maxent import closure1, closure2

RES = (1, 2, 4, 8, 16)
TAU = 1e-2
mpl.rcParams.update({"font.size": 8, "font.family": "serif",
                     "axes.linewidth": 0.6, "figure.dpi": 200})

rows = list(csv.DictReader(open("data/tm_E2_table.csv")))
members = []
for r in rows:
    if r["member"] not in members:
        members.append(r["member"])
cell = {(r["member"], int(r["ell"])): r for r in rows}

# categories: 0 H0, 1 2nd-closes, 2 2nd-improves, 3 higher
OI = {"H0": "#009E73", "closes": "#56B4E9", "improves": "#E69F00",
      "higher": "#D55E00"}
GLYPH = {"H0": "·", "closes": "o", "improves": "^", "higher": "x"}
order = ["H0", "closes", "improves", "higher"]


def cat(r):
    if r["exceed"] == "0":
        return "H0"
    if r["cls"] == "higher":
        return "higher"
    return "closes" if float(r["D2"]) <= TAU else "improves"


# ---------------- Fig 2: regime map + D1,D2 heatmaps ----------------
def fig2():
    ncat = np.zeros((len(members), len(RES)), dtype=int)
    D1 = np.zeros_like(ncat, dtype=float); D2 = np.zeros_like(ncat, dtype=float)
    for i, m in enumerate(members):
        for j, l in enumerate(RES):
            r = cell[(m, l)]
            ncat[i, j] = order.index(cat(r))
            D1[i, j] = float(r["D1"]); D2[i, j] = float(r["D2"])

    fig, axs = plt.subplots(1, 3, figsize=(7.0, 6.2),
                            gridspec_kw={"width_ratios": [1.25, 1, 1]})
    cmap = mpl.colors.ListedColormap([OI[c] for c in order])
    ax = axs[0]
    ax.imshow(ncat, aspect="auto", cmap=cmap, vmin=-0.5, vmax=3.5)
    for i in range(len(members)):
        for j in range(len(RES)):
            c = order[ncat[i, j]]
            ax.text(j, i, GLYPH[c], ha="center", va="center", fontsize=5,
                    color="black" if c in ("closes", "H0") else "white")
    ax.set_xticks(range(len(RES))); ax.set_xticklabels([f"$\\ell$={l}" for l in RES])
    ax.set_yticks(range(len(members)))
    ax.set_yticklabels(members, fontsize=5)
    ax.set_title("(a) regime map", fontsize=8)
    leg = [Patch(facecolor=OI[c], label=lab) for c, lab in
           [("H0", "mean sufficient"), ("closes", r"2nd order closes ($D_2\leq\tau$)"),
            ("improves", "2nd improves, not closed"), ("higher", "higher-order")]]
    ax.legend(handles=leg, loc="upper center", bbox_to_anchor=(1.9, -0.06),
              ncol=2, fontsize=6, frameon=False)

    for ax, Dm, ttl in ((axs[1], D1, "(b) $D_1$ (mean-only)"),
                        (axs[2], D2, "(c) $D_2$ (mean+var)")):
        im = ax.imshow(np.clip(Dm, 1e-4, None), aspect="auto", cmap="magma_r",
                       norm=mpl.colors.LogNorm(vmin=1e-3, vmax=1.0))
        ax.set_xticks(range(len(RES))); ax.set_xticklabels([f"{l}" for l in RES])
        ax.set_yticks([]); ax.set_title(ttl, fontsize=8)
        ax.set_xlabel("$\\ell$")
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.ax.tick_params(labelsize=5)
        # mark tau on colorbar
        cb.ax.axhline(TAU, color="#009E73", lw=1.2)
    fig.tight_layout()
    fig.savefig("fig2_regime_map.pdf", bbox_inches="tight")
    fig.savefig("fig2_regime_map.png", bbox_inches="tight")
    print("fig2 written")


# ---------------- Fig 3: representative curves ----------------
def curves(name, ell):
    nb = N // ell; etas = grid_eta(nb)
    mom = {(r["member"], int(r["ell"])): r for r in rows}
    mean = float(mom[(name, ell)]["EV"]); var = float(mom[(name, ell)]["Var"])
    p1 = closure1(nb, mean); p2, _ = closure2(nb, mean, var)
    pv = (pv_perm_cached(name, ell) if name.startswith("Perm")
          else pv_mix(float(name[4:-1]), ell) if name.startswith("Mix")
          else np.load("data/tm_PV_hat.npz")[f"{name}|{ell}"])

    def cv(pv):
        return np.array([pv @ w_table(nb, int(np.ceil(e * nb))) for e in etas])
    return etas, cv(pv), cv(p1), cv(p2 if p2 is not None else p1)


_pvc = {}
def pv_perm_cached(name, ell):
    if (name, ell) not in _pvc:
        from tm_core import build_family
        _pvc[(name, ell)] = pv_perm(dict(build_family())[name], ell)[0]
    return _pvc[(name, ell)]


def fig3():
    mpl.rcParams.update({"font.size": 10})
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.3))
    panels = [("Loc(8)", 8, "(a) Loc(8), $\\ell$=8 — 2nd order closes"),
              ("Perm(random-0)", 1, "(b) Perm(random-0), $\\ell$=1 — not closed")]
    for ax, (name, ell, ttl) in zip(axs, panels):
        etas, Rt, R1, R2 = curves(name, ell)
        ax.plot(etas, Rt, "-", color="black", lw=1.8, label="truth $R^*$")
        ax.plot(etas, R1, "--", color=OI["improves"], lw=1.6, label="closure-1 (mean)")
        ax.plot(etas, R2, ":", color=OI["closes"], lw=2.2, label="closure-2 (mean+var)")
        ax.set_xlabel("$\\eta$"); ax.set_title(ttl, fontsize=10)
        ax.set_ylabel("$R(\\eta)/R_0$")
        ax.legend(fontsize=9, frameon=False)
    fig.tight_layout()
    fig.savefig("fig3_curves.pdf", bbox_inches="tight")
    fig.savefig("fig3_curves.png", bbox_inches="tight")
    print("fig3 written")


if __name__ == "__main__":
    fig2(); fig3()
