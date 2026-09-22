#!/usr/bin/env python3
"""
Outcome stage 2 — freeze P̂(V) from M=1e6 trajectory samples, all 27 members.

ONE trajectory set per member (seed 2026400000 + 1e4·index); V(ℓ) for every ℓ
is derived from the SAME trajectories (§5.4). Class-appropriate exact samplers:
  Loc/Mix/ShiftMix  circulant → iid offset increments
  Perm*             deterministic orbit from the start position
  Block(M,c)        stepwise stay-in-block / jump-to-other-block
Saves data/tm_PV_hat.npz : arrays "<member>|<ell>" = P̂(V), v=1..nb.
"""
import numpy as np, time
from tm_core import N, T, build_family, block_of

M = 1_000_000
SEED = 2026400000
RES = (1, 2, 4, 8, 16)
CH = 50_000


def positions_circulant(P, rng, n):
    off = P[0].copy(); off /= off.sum()
    cdf = np.cumsum(off); cdf[-1] = 1.0
    inc = np.searchsorted(cdf, rng.random((n, T)))
    x0 = rng.integers(0, N, size=n)
    pos = np.empty((n, T + 1), dtype=np.int32)
    pos[:, 0] = x0
    pos[:, 1:] = (x0[:, None] + np.cumsum(inc, axis=1)) % N
    return pos


def orbit(perm):
    orb = np.empty((T + 1, N), dtype=np.int32)
    orb[0] = np.arange(N)
    for t in range(1, T + 1):
        orb[t] = perm[orb[t - 1]]
    return orb


def positions_perm(orb, rng, n):
    x0 = rng.integers(0, N, size=n)
    return orb[:, x0].T.copy()


def positions_block(M_blk, c, rng, n):
    bs = N // M_blk
    x = rng.integers(0, N, size=n)
    pos = np.empty((n, T + 1), dtype=np.int32)
    pos[:, 0] = x
    for t in range(T):
        cur = x // bs
        within = rng.random(n) < (1.0 - c)
        tgt = np.where(within, cur,
                       (cur + 1 + rng.integers(0, M_blk - 1, size=n)) % M_blk)
        x = tgt * bs + rng.integers(0, bs, size=n)
        pos[:, t + 1] = x
    return pos


def member_kind(name):
    if name.startswith("Perm"):
        return "perm"
    if name.startswith("Block"):
        return "block"
    return "circ"     # Loc, Mix, ShiftMix


def main():
    fam = build_family()
    out = {}
    for idx, (name, P) in enumerate(fam):
        t0 = time.time()
        rng = np.random.Generator(np.random.PCG64(SEED + 10_000 * idx))
        kind = member_kind(name)
        if kind == "perm":
            nxt = np.argmax(P, axis=1)
            orb = orbit(nxt)
        elif kind == "block":
            Mb = int(name[6:name.index(",")])
            c = float(name[name.index(",") + 1:-1])
        hist = {ell: np.zeros(N // ell + 1, dtype=np.int64) for ell in RES}
        done = 0
        while done < M:
            n = min(CH, M - done)
            if kind == "circ":
                pos = positions_circulant(P, rng, n)
            elif kind == "perm":
                pos = positions_perm(orb, rng, n)
            else:
                pos = positions_block(Mb, c, rng, n)
            for ell in RES:
                blk = pos // ell
                blk.sort(axis=1)
                cnt = 1 + (blk[:, 1:] != blk[:, :-1]).sum(axis=1)
                h = np.bincount(cnt, minlength=N // ell + 1)
                hist[ell][:len(h)] += h
            done += n
        for ell in RES:
            pv = hist[ell][1:N // ell + 1].astype(float) / M
            out[f"{name}|{ell}"] = pv
        print(f"[{idx:2d}] {name:16s} {kind:5s} done ({time.time()-t0:.0f}s)  "
              f"E[V](l16)={(np.arange(1,17)*out[name+'|16']).sum():.4f}", flush=True)
    np.savez_compressed("data/tm_PV_hat.npz", **out)
    print("saved data/tm_PV_hat.npz")


if __name__ == "__main__":
    main()
