"""Second-order coefficient d^2 a_m/dt^2|_0 of the non-harmonic multi-valley gate.

Continues branching_nonsplittable_general.py (which fixed the harmonic anchor
a_m(0) = 2 sin(pi m/k) and the first-order slope a_m'(0)). Here we obtain the
curvature a_m''(0), the last elementary-to-evaluate piece of the perturbation
expansion, via a 2D resolvent solve.

Setup.  E(t) = S_m a_m(t) = min_h sum_edges c_pq(t)(h_p-h_q)^2, c_pq=0.5(w_p+w_q),
w_p=exp(-F_t(p)), F_t = s(-Re w^k + t r^k), S_m = sum_j cos^2(2 pi m j/k).
Let L(t) be the weighted graph Laplacian (H^T L H = sum_edges c (dh)^2).  With
'=d/dt at fixed h  (c' = 0.5(w_p'+w_q'), w_p' = -s r_p^k w_p;
                    c''= 0.5(w_p''+w_q''), w_p''= (s r_p^k)^2 w_p):

    E'(0)  = [h]^T L' [h],
    E''(0) = [h]^T L'' [h] + 2 [hdot]^T L' [h],   hdot = 0 on reservoirs,
    L_II hdot_I = -(L'[h])_I               <-- the 2D resolvent solve.

a_m''(0) = E''(0)/S_m.  Validated against the central second difference of a_m(t);
the same LU factorization of L_II serves the committor and the resolvent.

Taylor:  a_m(t) = 2 sin(pi m/k) + a_m'(0) t + 1/2 a_m''(0) t^2 + O(t^3).
For k=4 the coupling map a=(6+2t)/(1-t) (t=(a-6)/(a+2), da/dt|_0=8,
d^2a/dt^2|_0 = ... ; dt/da|_6=1/8, d^2t/da^2|_6=-1/32) gives, by the chain rule,
    d^2 a_m/da^2|_6 = a_m''(0)/64 - a_m'(0)/32.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix, csc_matrix
from scipy.sparse.linalg import splu, spsolve

FIGURE_PATH = Path("figures/branching_nonsplittable_second_order.png")


def _grid(k, t, m, s, L, n, rad):
    xr = L - rad - 0.1
    wells = [(xr * math.cos(2 * math.pi * j / k), xr * math.sin(2 * math.pi * j / k))
             for j in range(k)]
    vals = np.array([math.cos(2 * math.pi * m * j / k) for j in range(k)])
    xs = np.linspace(-L, L, n)
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    rk = (X ** 2 + Y ** 2) ** (k / 2)
    F = s * (-np.real((X + 1j * Y) ** k) + t * rk)
    w = np.exp(-F)
    dead = w < 1e-250
    wid = np.full((n, n), -1)
    for j, (cx, cy) in enumerate(wells):
        wid[(X - cx) ** 2 + (Y - cy) ** 2 <= rad ** 2] = j
    return X, Y, rk, w, dead, wid, vals


def admittance(k, t, m, s=1.0, L=1.9, n=541, rad=0.20):
    """Plain committor admittance a_m(t) (for finite-difference cross-checks)."""
    _, _, _, w, dead, wid, vals = _grid(k, t, m, s, L, n, rad)
    isf = wid >= 0
    idx = -np.ones((n, n), int)
    unk = np.argwhere(~isf & ~dead)
    for kk, (i, j) in enumerate(unk):
        idx[i, j] = kk
    rows, cols, vv, rhs = [], [], [], np.zeros(len(unk))
    for p, (i, j) in enumerate(unk):
        dg = 0.0
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            ii, jj = i + di, j + dj
            if not (0 <= ii < n and 0 <= jj < n) or dead[ii, jj]:
                continue
            ae = 0.5 * (w[i, j] + w[ii, jj]); dg -= ae
            if isf[ii, jj]:
                rhs[p] -= ae * vals[wid[ii, jj]]
            else:
                rows.append(p); cols.append(idx[ii, jj]); vv.append(ae)
        rows.append(p); cols.append(p); vv.append(dg)
    sol = spsolve(csr_matrix((vv, (rows, cols)), shape=(len(unk),) * 2), rhs)
    H = np.where(isf, vals[wid], 0.0).astype(float)
    for p, (i, j) in enumerate(unk):
        H[i, j] = sol[p]
    I0 = 0.0
    for (i, j) in np.argwhere(wid == 0):
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            ii, jj = i + di, j + dj
            if not (0 <= ii < n and 0 <= jj < n) or dead[ii, jj] or wid[ii, jj] == 0:
                continue
            I0 += 0.5 * (w[i, j] + w[ii, jj]) * (H[i, j] - H[ii, jj])
    return I0


def perturbation(k, m, s=1.0, L=1.9, n=541, rad=0.20):
    """Return (a0, a1, a2) = a_m(0), a_m'(0), a_m''(0) via the graph-Laplacian resolvent."""
    X, Y, rk, w, dead, wid, vals = _grid(k, 0.0, m, s, L, n, rad)
    live = ~dead
    isf = wid >= 0
    gid = -np.ones((n, n), int)
    live_nodes = np.argwhere(live)
    gid[live] = np.arange(len(live_nodes))
    N = len(live_nodes)

    wp = w                       # w
    wp1 = -s * rk * w            # dw/dt
    wp2 = (s * rk) ** 2 * w      # d^2w/dt^2

    def edges_axis(axis):
        a_id = gid[:-1, :] if axis == 0 else gid[:, :-1]
        b_id = gid[1:, :] if axis == 0 else gid[:, 1:]
        mask = (a_id >= 0) & (b_id >= 0)
        return a_id[mask], b_id[mask], mask

    gh_a, gh_b, mh = edges_axis(0)
    gv_a, gv_b, mv = edges_axis(1)
    ea = np.concatenate([gh_a, gv_a]); eb = np.concatenate([gh_b, gv_b])

    def assemble(node_w):
        ch = 0.5 * (node_w[:-1, :] + node_w[1:, :])[mh]
        cv = 0.5 * (node_w[:, :-1] + node_w[:, 1:])[mv]
        c = np.concatenate([ch, cv])
        R = np.concatenate([ea, eb, ea, eb]); C = np.concatenate([ea, eb, eb, ea])
        V = np.concatenate([c, c, -c, -c])
        return csr_matrix((V, (R, C)), shape=(N, N))

    Lmat, Lp, Lpp = assemble(wp), assemble(wp1), assemble(wp2)

    safe_wid = np.where(wid >= 0, wid, 0)
    Hgrid = np.where(isf, vals[safe_wid], 0.0).astype(float)
    inter = gid[live & ~isf]
    bnd = gid[live & isf]
    Hvec = np.zeros(N)
    Hvec[gid[live]] = Hgrid[live]

    Lii = csc_matrix(Lmat[inter][:, inter])
    Lib = Lmat[inter][:, bnd]
    lu = splu(Lii)
    Hvec[inter] = lu.solve(-(Lib @ Hvec[bnd]))     # harmonic committor

    Sm = float(np.sum(vals ** 2))
    u = Lp @ Hvec                                  # L'[h]
    Ep = float(Hvec @ u)
    hdot_I = lu.solve(-u[inter])                   # 2D resolvent
    Epp = float(Hvec @ (Lpp @ Hvec) + 2.0 * (hdot_I @ u[inter]))
    E = float(Hvec @ (Lmat @ Hvec))
    return E / Sm, Ep / Sm, Epp / Sm


def main() -> None:
    print("=== Perturbation of a_m at the harmonic point t=0 (resolvent vs FD) ===")
    print(f"{'k':>2} {'m':>2} {'a0':>9} {'a1(res)':>9} {'a1(FD)':>9} "
          f"{'a2(res)':>9} {'a2(FD)':>9}")
    hstep = 0.14
    results = {}
    for k, s in [(4, 1.0), (6, 0.20)]:
        for m in range(1, k // 2 + 1):
            a0, a1, a2 = perturbation(k, m, s)
            ap = admittance(k, hstep, m, s)
            am = admittance(k, -hstep, m, s)
            fd1 = (ap - am) / (2 * hstep)
            fd2 = (ap - 2 * a0 + am) / hstep ** 2
            results[(k, m)] = (a0, a1, a2)
            print(f"{k:>2} {m:>2} {a0:>9.4f} {a1:>9.4f} {fd1:>9.4f} "
                  f"{a2:>9.4f} {fd2:>9.4f}")

    print("\n=== Taylor a_m(t) = 2sin(pi m/k) + a1 t + 1/2 a2 t^2 ===")
    for (k, m), (a0, a1, a2) in results.items():
        anchor = 2 * math.sin(math.pi * m / k)
        print(f" k={k} m={m}: a_m(t) ~ {anchor:.4f} + ({a1:+.4f}) t "
              f"+ 1/2 ({a2:+.4f}) t^2")

    print("\n=== k=4: convert to coupling a via a=(6+2t)/(1-t) (chain rule at a=6) ===")
    for m in [1, 2]:
        a0, a1, a2 = results[(4, m)]
        da_da = a1 / 8.0                       # dt/da|_6 = 1/8
        d2a_da2 = a2 / 64.0 - a1 / 32.0        # a2*(1/8)^2 + a1*(-1/32)
        anchor = 2 * math.sin(math.pi * m / 4)
        print(f" m={m}: a_m(a) ~ {anchor:.4f} + ({da_da:+.4f})(a-6) "
              f"+ 1/2 ({d2a_da2:+.4f})(a-6)^2")

    # figure: k=4 data with the osculating parabola from (a0, a1, a2)
    ts = np.linspace(-0.5, 0.5, 11)
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    tt = np.linspace(-0.55, 0.55, 200)
    colors = {1: "#0072B2", 2: "#D55E00"}
    for m in [1, 2]:
        a0, a1, a2 = results[(4, m)]
        data = [admittance(4, t, m) for t in ts]
        ax.plot(ts, data, "o", color=colors[m], ms=6, label=fr"$a_{m}(t)$ (solver)")
        par = a0 + a1 * tt + 0.5 * a2 * tt ** 2
        ax.plot(tt, par, "-", color=colors[m], lw=1.6, alpha=0.9,
                label=fr"2nd-order Taylor $m={m}$")
        ax.plot(0, 2 * math.sin(math.pi * m / 4), "*", color=colors[m], ms=15, zorder=5)
    ax.axvline(0, color="#888888", ls="--", lw=1.0)
    ax.annotate("harmonic $t=0$", (0.0, 0.82), ha="center", fontsize=9, color="#555")
    ax.set_xlabel(r"anharmonicity  $t$  in  $-\mathrm{Re}\,w^4 + t\,(x^2+y^2)^2$")
    ax.set_ylabel(r"channel admittance  $a_m$")
    ax.set_title("Second-order expansion: osculating parabola at the harmonic point",
                 fontweight="bold", fontsize=11)
    ax.grid(True, color="#eaeaea")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=8, ncol=2)
    fig.tight_layout()
    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_PATH, dpi=150)
    print(f"\nSaved {FIGURE_PATH}")


if __name__ == "__main__":
    main()
