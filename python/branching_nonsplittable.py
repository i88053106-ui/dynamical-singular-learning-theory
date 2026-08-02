"""Branching non-splittable gate: coefficient-dependent circulant admittances.

Case II (branching, harmonic Re w^k) gave universal admittances a_m = 2 sin(pi m/k),
independent of any coefficient. Case III (non-splittable) gave coefficient-dependent
capacity. Here we merge them with the homogeneous degree-4 family

    P_a(x,y) = -x^4 + a x^2 y^2 - y^4      (a > 2: four valleys on +-x, +-y),

which is 4-fold symmetric (a circulant 4-terminal junction, Case II structure) and
harmonic ONLY at a = 6, where P_6 = -(x^4 - 6 x^2 y^2 + y^4) = -Re(w^4). So:

  * at a = 6 the gate is the harmonic k=4 branching gate and the scale-invariant
    channel admittances are the universal Steklov values a_m = 2 sin(pi m/4);
  * for a != 6 the gate is non-harmonic, and the admittances a_m(a) DEPEND on the
    coefficient a -- the branching network survives, but its admittances become
    PDE-algebraic (Case III) rather than universal.

a_m(a) is scale invariant (P_a homogeneous), so we solve the pure germ at T=1 with
four reservoirs deep in the valleys and read the current into valley 0 under the
Fourier boundary data cos(2 pi m j / 4).
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

FIGURE_PATH = Path("figures/branching_nonsplittable.png")


def channel_admittance(a, m, L=1.9, n=541, rad=0.20):
    """Scale-invariant a_m for P_a = -x^4 + a x^2 y^2 - y^4 (pure germ, T=1)."""
    xr = L - rad - 0.1
    wells = [(xr, 0.0), (0.0, xr), (-xr, 0.0), (0.0, -xr)]
    F = lambda x, y: -x ** 4 + a * x ** 2 * y ** 2 - y ** 4
    vals = np.array([math.cos(2 * math.pi * m * j / 4) for j in range(4)])

    xs = np.linspace(-L, L, n)
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    w = np.exp(-F(X, Y))
    dead = w < 1e-250

    wid = np.full((n, n), -1)
    for j, (cx, cy) in enumerate(wells):
        wid[(X - cx) ** 2 + (Y - cy) ** 2 <= rad ** 2] = j
    isf = wid >= 0

    idx = -np.ones((n, n), int)
    unk = np.argwhere(~isf & ~dead)
    for k, (i, j) in enumerate(unk):
        idx[i, j] = k
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


def main() -> None:
    s1, s2 = 2 * math.sin(math.pi / 4), 2 * math.sin(math.pi / 2)  # 2sin(pi m/4)
    print(f"Harmonic Steklov values: a_1=2sin(pi/4)={s1:.4f}, a_2=2sin(pi/2)={s2:.4f}")

    print("\nConvergence at a=6 (=-Re w^4), reservoirs deepening -> a_1 -> 2sin(pi/4):")
    for L, n, rad in [(1.3, 421, 0.16), (1.6, 481, 0.18), (1.9, 541, 0.20)]:
        print(f"  L={L}, n={n}: a_1 = {channel_admittance(6, 1, L, n, rad):.4f}")

    print("\nAdmittances vs coupling a (L=1.9, n=541):")
    print(f"{'a':>4} {'a_1':>8} {'a_2':>8}   note")
    a_list = [3, 4, 5, 6, 7, 8]
    A1 = [channel_admittance(a, 1) for a in a_list]
    A2 = [channel_admittance(a, 2) for a in a_list]
    for a, a1, a2 in zip(a_list, A1, A2):
        note = "<- harmonic: 2sin = (1.414, 2.000)" if a == 6 else ""
        print(f"{a:>4} {a1:>8.4f} {a2:>8.4f}   {note}")
    print("\nThe admittances depend on a: the branching network survives, but its")
    print("admittances are coefficient-dependent (Case III) away from the harmonic")
    print("point a=6, where they equal the universal Steklov values 2 sin(pi m/4).")

    fig, ax = plt.subplots(figsize=(6.8, 5.0))
    ax.axhline(s1, color="#0072B2", ls=":", lw=1.2)
    ax.axhline(s2, color="#D55E00", ls=":", lw=1.2)
    ax.plot(a_list, A1, "o-", color="#0072B2", ms=7, label=r"$a_1(a)$")
    ax.plot(a_list, A2, "s-", color="#D55E00", ms=7, label=r"$a_2(a)$")
    ax.axvline(6, color="#888888", ls="--", lw=1.0)
    ax.annotate("harmonic\n$a=6$ ($-\\mathrm{Re}\\,w^4$)", (6, 3.6), ha="center",
                fontsize=9, color="#555555")
    ax.annotate(r"$2\sin(\pi/2)$", (8.1, s2), fontsize=9, color="#D55E00", va="center")
    ax.annotate(r"$2\sin(\pi/4)$", (8.1, s1), fontsize=9, color="#0072B2", va="center")
    ax.set_xlabel(r"coupling  $a$  in  $-x^4 + a\,x^2y^2 - y^4$")
    ax.set_ylabel(r"channel admittance  $a_m$")
    ax.set_title("Branching non-splittable: admittances depend on the coefficient",
                 fontweight="bold")
    ax.grid(True, color="#eaeaea")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    fig.tight_layout()
    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_PATH, dpi=150)
    print(f"Saved {FIGURE_PATH}")


if __name__ == "__main__":
    main()
