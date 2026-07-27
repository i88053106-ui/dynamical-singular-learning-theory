"""Large-k asymptotics of the branching-gate channel admittances a_m.

The junction of a degree-k homogeneous gate is a circulant k-terminal network;
its channel admittances are lambda_m = a_m T (scale invariant a_m). We measure
a_m directly: fix ALL k wells to the Fourier mode v_j = cos(2 pi m j / k) and read
the net current out of well 0, which equals the eigenvalue lambda_m.

Large-k prediction. Zooming near the unit circle (u = k(rho-1), psi = k phi) sends
the gate to the universal cylinder weight e^{e^u cos psi}; the circulant spectrum
becomes a dispersion relation Lambda(q) = sum_s beta_s (1 - cos q s) at wavenumber
q = 2 pi m / k. If the chain couplings decay as beta_s ~ (2/pi)/s^2 then
Lambda(q) ~ |q| for small q, i.e.

    a_m(k)  ~  2 pi m / k    (fixed m, k -> infinity).

This script tests that law (a_m k / (2 pi m) -> 1) and the coupling decay
beta_s ~ (2/pi)/s^2.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

LDOM, RADIUS = 2.0, 0.24


def make_F(k):
    def F(x, y):
        return (x * x + y * y) ** k / (2.0 * k) - np.real((x + 1j * y) ** k)
    return F


def well_centers(k):
    r0 = k ** (1.0 / k)
    return [(r0 * math.cos(2 * math.pi * j / k), r0 * math.sin(2 * math.pi * j / k))
            for j in range(k)]


def channel_admittance(k, m, T, n=361):
    """Return a_m = lambda_m / T via the all-terminal Fourier-mode response."""
    F = make_F(k)
    centers = well_centers(k)
    values = np.array([math.cos(2 * math.pi * m * j / k) for j in range(k)])

    xs = np.linspace(-LDOM, LDOM, n)
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    a = np.exp(-F(X, Y) / T)
    dead = a < 1e-250

    wid = np.full((n, n), -1, dtype=int)
    r2 = RADIUS ** 2
    for j, (cx, cy) in enumerate(centers):
        wid[(X - cx) ** 2 + (Y - cy) ** 2 <= r2] = j
    is_fixed = wid >= 0

    idx = -np.ones((n, n), dtype=int)
    unknown = np.argwhere(~is_fixed & ~dead)
    for kk, (i, j) in enumerate(unknown):
        idx[i, j] = kk
    n_unk = len(unknown)

    rows, cols, vals = [], [], []
    rhs = np.zeros(n_unk)
    for p, (i, j) in enumerate(unknown):
        diag = 0.0
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            ii, jj = i + di, j + dj
            if not (0 <= ii < n and 0 <= jj < n) or dead[ii, jj]:
                continue
            ae = 0.5 * (a[i, j] + a[ii, jj])
            diag -= ae
            if is_fixed[ii, jj]:
                rhs[p] -= ae * values[wid[ii, jj]]
            else:
                rows.append(p); cols.append(idx[ii, jj]); vals.append(ae)
        rows.append(p); cols.append(p); vals.append(diag)

    A = csr_matrix((vals, (rows, cols)), shape=(n_unk, n_unk))
    sol = spsolve(A, rhs)

    H = np.where(is_fixed, values[wid], 0.0).astype(float)
    for p, (i, j) in enumerate(unknown):
        H[i, j] = sol[p]

    # net current out of well 0 = sum over edges (well0 node -> outside) a_edge (v0 - h)
    I0 = 0.0
    for (i, j) in np.argwhere(wid == 0):
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            ii, jj = i + di, j + dj
            if not (0 <= ii < n and 0 <= jj < n) or dead[ii, jj]:
                continue
            if wid[ii, jj] == 0:
                continue
            ae = 0.5 * (a[i, j] + a[ii, jj])
            I0 += ae * (H[i, j] - H[ii, jj])
    # cap = T * int |grad h|^2 e^{-F/T} (Dirichlet form carries the T), so the
    # scale-invariant admittance a_m = lambda_m/T equals the raw current I0.
    return I0


def main() -> None:
    import matplotlib.pyplot as plt

    T = 0.12
    ks = [3, 4, 5, 6, 8, 10, 12, 14]
    data = {1: {}, 2: {}, 3: {}}
    for k in ks:
        for m in (1, 2, 3):
            if 2 * m < k:  # genuine low-wavenumber channel (q = 2 pi m/k < pi)
                data[m][k] = channel_admittance(k, m, T)

    print(f"=== channel admittances a_m and the law a_m ~ 2 pi m / k  (T={T}) ===")
    print(f"{'k':>3} {'a_1':>8} {'a_2':>8} {'a_3':>8} "
          f"{'a1·k/2π':>9} {'a2·k/4π':>9} {'a3·k/6π':>9}")
    for k in ks:
        a = {m: data[m].get(k, float('nan')) for m in (1, 2, 3)}
        c = {m: a[m] * k / (2 * math.pi * m) for m in (1, 2, 3)}
        print(f"{k:>3} {a[1]:>8.4f} {a[2]:>8.4f} {a[3]:>8.4f} "
              f"{c[1]:>9.4f} {c[2]:>9.4f} {c[3]:>9.4f}")

    print("\nIf a_m ~ 2 pi m / k then each c_m = a_m k /(2 pi m) -> 1 as k grows.")
    print("Leading slope 1 <=> dispersion Lambda(q) ~ |q|, i.e. a chain-coupling")
    print("tail beta_s ~ (2/pi)/s^2 (an asymptotic tail; individual finite-k b_s differ).")

    from pathlib import Path
    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    colors = {1: "#0072B2", 2: "#D55E00", 3: "#009E73"}
    ax.axhline(1.0, color="#888888", lw=1.2, ls="--", label="limit 1")
    for m in (1, 2, 3):
        kk = sorted(data[m])
        cc = [data[m][k] * k / (2 * math.pi * m) for k in kk]
        ax.plot(kk, cc, "o-", color=colors[m], ms=7, label=f"m={m}")
    ax.set_xlabel("valley number  k")
    ax.set_ylabel(r"$a_m \, k / (2\pi m)$")
    ax.set_title(r"Large-$k$ admittance law:  $a_m \sim 2\pi m/k$", fontweight="bold")
    ax.grid(True, color="#eaeaea")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    fig.tight_layout()
    out = Path("figures/branching_admittance_largek.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
