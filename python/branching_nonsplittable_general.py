"""General-k non-harmonic multi-valley branching gate: anchor, slope, scale reduction.

Reparametrize the homogeneous branching family through the harmonic point:

    P_{k,t}(x,y) = -Re(w^k) + t * (x^2 + y^2)^{k/2},   w = x + i y,   k even,

so that t = 0 is the harmonic gate -Re(w^k) (Case II) and t sweeps the
anharmonicity. The added term t r^k is isotropic, so it does NOT move the
valleys (theta_j = 2 pi j / k): the circulant k-terminal network is preserved
for every t, only the channel admittances a_m change.

For k = 4 this is exactly the family of `branching_nonsplittable.py`:
    -Re(w^4) + t r^4 = (t-1) x^4 + (6+2t) x^2 y^2 + (t-1) y^4
                     proportional to  -x^4 + a x^2 y^2 - y^4,   a = (6+2t)/(1-t),
i.e.  t = (a-6)/(a+2)   (a=6 <-> t=0 harmonic;  a=2 <-> t=-1;  a->inf <-> t->1).
Since a_m is scale invariant for a homogeneous germ (P -> cP), a_m depends on the
single variable t; the two-parameter look of "a_m(a)" collapses to one curve.

Three results (see docs/branching-nonsplittable-general.md):

  (a) Closed form of a_m(a):  there is NO elementary closed form for t != 0.
      a_m(t) is a PDE-algebraic invariant (Case III). The only exact anchor is
      the harmonic value a_m(0) = 2 sin(pi m / k).

  (b) Harmonic-point perturbation da_m/da|_{a=6}:  the envelope theorem gives
          da_m/dt|_0 = -(1/S_m) * sum_edges 0.5(w_p r_p^k + w_q r_q^k)(dh_m)^2,
      S_m = sum_j cos^2(2 pi m j/k)  ( = k/2 generically, = k at the Nyquist mode
      m=k/2 ), h_m = harmonic committor of mode m. This finite-difference-free
      slope matches the central difference of a_m(t). For k=4, da_m/da = (1/8) da_m/dt.

  (c) General k:  the harmonic anchor a_m(0) = 2 sin(pi m/k) is universal for all
      k, but a_m(t) is coefficient-dependent for t != 0 and every k. Verified k=4,6.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

FIGURE_PATH = Path("figures/branching_nonsplittable_general.png")


def _solve(k, t, m, s, L, n, rad):
    """Weighted committor for P_{k,t} scaled by fixed s; return grid quantities."""
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
    return X, Y, w, dead, wid, H, rk, vals


def channel_admittance(k, t, m, s=1.0, L=1.9, n=541, rad=0.20):
    """Scale-invariant channel admittance a_m of P_{k,t} = -Re w^k + t r^k (T=1).

    a_m = current into valley 0 under Fourier boundary data cos(2 pi m j/k).
    s rescales the germ (P -> sP, a_m invariant) for numerical range control.
    """
    _, _, w, dead, wid, H, _, _ = _solve(k, t, m, s, L, n, rad)
    n = w.shape[0]
    I0 = 0.0
    for (i, j) in np.argwhere(wid == 0):
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            ii, jj = i + di, j + dj
            if not (0 <= ii < n and 0 <= jj < n) or dead[ii, jj] or wid[ii, jj] == 0:
                continue
            I0 += 0.5 * (w[i, j] + w[ii, jj]) * (H[i, j] - H[ii, jj])
    return I0


def slope_envelope(k, m, s=1.0, L=1.9, n=541, rad=0.20):
    """da_m/dt at t=0 by the discrete envelope theorem (finite-difference-free).

    E(t) = min_h sum_edges 0.5(w_p+w_q)(dh)^2 = S_m * a_m(t),  S_m = sum_j vals_j^2.
    dE/dt|_0 = sum_edges 0.5(wdot_p+wdot_q)(dh)^2,  wdot = -w s r^k.
    => da_m/dt|_0 = dE/dt|_0 / S_m.
    """
    X, Y, w, dead, wid, H, rk, vals = _solve(k, 0.0, m, s, L, n, rad)
    n = w.shape[0]
    wdot = -w * s * rk
    dEdt = 0.0
    for (i, j) in np.argwhere(~dead):
        for di, dj in ((1, 0), (0, 1)):   # each edge once
            ii, jj = i + di, j + dj
            if not (0 <= ii < n and 0 <= jj < n) or dead[ii, jj]:
                continue
            dh = H[i, j] - H[ii, jj]
            dEdt += 0.5 * (wdot[i, j] + wdot[ii, jj]) * dh * dh
    Sm = float(np.sum(vals ** 2))
    return dEdt / Sm


def a_of_t(t):
    return (6 + 2 * t) / (1 - t)


def main() -> None:
    # (a) k=4 reparametrization: a_m(t) is a single-variable curve; a maps to t.
    print("=== (a) k=4: a_m as a single-variable curve a_m(t), a=(6+2t)/(1-t) ===")
    print(f"{'t':>6} {'a':>7} {'a_1':>8} {'a_2':>8}")
    for t in [-0.5, -0.25, 0.0, 0.25, 0.5]:
        a1 = channel_admittance(4, t, 1)
        a2 = channel_admittance(4, t, 2)
        print(f"{t:>6.2f} {a_of_t(t):>7.3f} {a1:>8.4f} {a2:>8.4f}")
    print("Monotone decreasing in t (increasing coupling a) -> coefficient dependent.")

    # (b) harmonic anchor a_m(0) = 2 sin(pi m/k)
    print("\n=== (b) harmonic anchor a_m(0) = 2 sin(pi m/k) ===")
    for k, s in [(4, 1.0), (6, 0.20)]:
        for m in range(1, k // 2 + 1):
            am = channel_admittance(k, 0.0, m, s)
            print(f" k={k} m={m}: a_m(0)={am:.4f}   2sin(pi m/k)={2*math.sin(math.pi*m/k):.4f}")

    # (b/c) slope da_m/dt|_0: envelope theorem vs central difference
    print("\n=== (b) slope da_m/dt|_0: envelope theorem vs central difference ===")
    hstep = 0.12
    for k, s in [(4, 1.0), (6, 0.20)]:
        for m in range(1, k // 2 + 1):
            env = slope_envelope(k, m, s)
            ap = channel_admittance(k, hstep, m, s)
            am = channel_admittance(k, -hstep, m, s)
            fd = (ap - am) / (2 * hstep)
            extra = ""
            if k == 4:
                extra = f"   da_{m}/da|_6 = (1/8)*slope = {env/8:+.4f}"
            print(f" k={k} m={m}: envelope={env:+.4f}  central-diff={fd:+.4f}{extra}")
    print("\nFor k=4, da/dt|_0 = 8, so da_m/da|_{a=6} = (1/8) da_m/dt|_0.")

    # figure: k=4 curve a_m(t) with harmonic anchors and slope tangents
    ts = np.linspace(-0.55, 0.55, 11)
    A1 = [channel_admittance(4, t, 1) for t in ts]
    A2 = [channel_admittance(4, t, 2) for t in ts]
    s1, s2 = 2 * math.sin(math.pi / 4), 2 * math.sin(math.pi / 2)
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.plot(ts, A1, "o-", color="#0072B2", ms=6, label=r"$a_1(t)$")
    ax.plot(ts, A2, "s-", color="#D55E00", ms=6, label=r"$a_2(t)$")
    ax.plot(0, s1, "*", color="#0072B2", ms=16, zorder=5)
    ax.plot(0, s2, "*", color="#D55E00", ms=16, zorder=5)
    ax.annotate(r"$2\sin(\pi/4)$", (0.06, s1), color="#0072B2", fontsize=9, va="center")
    ax.annotate(r"$2\sin(\pi/2)$", (0.06, s2), color="#D55E00", fontsize=9, va="center")
    ax.axvline(0, color="#888888", ls="--", lw=1.0)
    ax.annotate("harmonic\n$t=0$ ($-\\mathrm{Re}\\,w^4$)", (0.0, 0.85), ha="center",
                fontsize=9, color="#555555")
    ax.set_xlabel(r"anharmonicity  $t$  in  $-\mathrm{Re}\,w^4 + t\,(x^2+y^2)^2$")
    ax.set_ylabel(r"channel admittance  $a_m$")
    ax.set_title("Non-harmonic multi-valley: single-variable $a_m(t)$, "
                 "universal only at $t=0$", fontweight="bold", fontsize=11)
    ax.grid(True, color="#eaeaea")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    fig.tight_layout()
    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_PATH, dpi=150)
    print(f"\nSaved {FIGURE_PATH}")


if __name__ == "__main__":
    main()
