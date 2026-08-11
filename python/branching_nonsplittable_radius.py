"""Convergence radius and large-coupling regime of the non-harmonic multi-valley gate.

The single-variable admittance a_m(t) of P_{4,t} = -Re(w^4) + t (x^2+y^2)^2
(equivalently -x^4 + a x^2 y^2 - y^4 with a=(6+2t)/(1-t)) is a genuine isolated
4-valley gate only for t in (-1,1) <-> a in (2, inf).  At both ends the germ,
though still degree-4 homogeneous, loses its isolated-well structure:

    t -> +1 :  P -> 8 x^2 y^2        (the axes become flat valleys; Newton-degenerate),
    t -> -1 :  P -> -(x^2 - y^2)^2   (the diagonal ridges collapse; wells coalesce).

So a_m(t) is real-analytic on (-1,1) with algebraic branch points at t = +-1.
Because both singularities sit at |t| = 1, the Taylor series of a_m about the
harmonic point t=0 (built in branching-nonsplittable-second-order.md) has radius
of convergence exactly R = 1.

Measured endpoint laws (k=4):
    t -> +1 :  a_m ~ C_+ (1-t)^{gamma_+},  gamma_+ small (slow decay);
               in the coupling,  a_m ~ const * a^{-gamma_+}   as a -> inf.
    t -> -1 :  a_m ~ C_- (1+t)^{-1/2}      (square-root divergence);
               since 1+t ~ (a-2)/2,        a_m ~ C' (a-2)^{-1/2} as a -> 2+.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from branching_nonsplittable_general import channel_admittance as ca

FIGURE_PATH = Path("figures/branching_nonsplittable_radius.png")


def endpoint_fit(ts, vals):
    """Fit log(vals) = gamma*log(dist) + c; return (gamma, C)."""
    lx, ly = np.log(ts), np.log(vals)
    g, c = np.polyfit(lx, ly, 1)
    return g, math.exp(c)


def main() -> None:
    # ---- full-range sweep --------------------------------------------------
    ts = [-0.97, -0.95, -0.9, -0.8, -0.7, -0.6, -0.4, -0.2, 0.0,
          0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95, 0.97]
    A = {1: {}, 2: {}}
    print("=== a_m(t), k=4 (L=1.9, n=541) ===")
    print(f"{'t':>7} {'a':>10} {'a_1':>9} {'a_2':>9}")
    for t in ts:
        for m in (1, 2):
            A[m][t] = ca(4, t, m)
        a = (6 + 2 * t) / (1 - t)
        print(f"{t:>7.2f} {a:>10.3f} {A[1][t]:>9.4f} {A[2][t]:>9.4f}")

    # ---- endpoint t -> +1 (large coupling) --------------------------------
    print("\n=== t -> +1 (a -> inf):  a_m ~ C_+ (1-t)^{gamma_+} ===")
    tp = [0.7, 0.8, 0.9, 0.95, 0.97]
    gpluses = {}
    for m in (1, 2):
        gp, Cp = endpoint_fit([1 - t for t in tp], [A[m][t] for t in tp])
        gpluses[m] = (gp, Cp)
        print(f" m={m}: gamma_+ = {gp:.3f},  C_+ = {Cp:.3f}   "
              f"=>  a_m ~ {Cp*8**gp:.3f} * a^(-{gp:.3f})  (a->inf)")

    # ---- endpoint t -> -1 (a -> 2): square-root divergence ----------------
    print("\n=== t -> -1 (a -> 2):  a_m ~ C_- (1+t)^{gamma_-} ===")
    tm = [-0.7, -0.8, -0.9, -0.95, -0.97]
    gminuses = {}
    for m in (1, 2):
        gm, Cm = endpoint_fit([1 + t for t in tm], [A[m][t] for t in tm])
        gminuses[m] = (gm, Cm)
        print(f" m={m}: gamma_- = {gm:.3f}  (cf. -1/2),  C_- = {Cm:.3f}")
    # local two-point exponents near -1 (grid-convergence sanity)
    print(" local exponents near t=-1 (should approach -1/2):")
    for m in (1, 2):
        pairs = [(-0.9, -0.95), (-0.95, -0.97)]
        locs = []
        for ta, tb in pairs:
            loc = math.log(A[m][tb] / A[m][ta]) / math.log((1 + tb) / (1 + ta))
            locs.append(loc)
        print(f"   m={m}: " + ", ".join(f"{l:+.3f}" for l in locs))

    # ---- radius of convergence: coefficient ratios ------------------------
    # low-order Taylor coefficients from a symmetric least-squares fit on
    # |t|<=0.6 (well inside R); |c_n/c_{n-1}| -> 1/R = 1.
    print("\n=== Taylor coefficients about t=0 and ratio test (|c_n/c_{n-1}| -> 1/R) ===")
    win = [t for t in ts if abs(t) <= 0.6]
    dombsykes = {}
    for m in (1, 2):
        coef = np.polyfit(win, [A[m][t] for t in win], 6)[::-1]  # c0..c6
        ratios = [abs(coef[n] / coef[n - 1]) for n in range(1, 6)]
        dombsykes[m] = ratios
        print(f" m={m}: c0..c4 = " + ", ".join(f"{coef[i]:+.4f}" for i in range(5)))
        print(f"        |c_n/c_(n-1)|, n=1..5: " + ", ".join(f"{r:.3f}" for r in ratios))
    print(" ratios are O(1) and oscillate (two singularities of opposite sign at t=+-1);")
    print(" this is consistent with R~1 but the decisive evidence is the branch points.")

    # ---- figure ------------------------------------------------------------
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.5, 4.8))
    tarr = np.array(ts)
    colors = {1: "#0072B2", 2: "#D55E00"}
    for m in (1, 2):
        y = np.array([A[m][t] for t in ts])
        axL.semilogy(tarr, y, "o-", color=colors[m], ms=5, label=fr"$a_{m}(t)$")
        gm, Cm = gminuses[m]
        tt = np.linspace(-0.985, -0.55, 60)
        axL.semilogy(tt, Cm * (1 + tt) ** gm, "--", color=colors[m], lw=1.2, alpha=0.8)
        gp, Cp = gpluses[m]
        tu = np.linspace(0.55, 0.985, 60)
        axL.semilogy(tu, Cp * (1 - tu) ** gp, ":", color=colors[m], lw=1.4, alpha=0.9)
        axL.plot(0, 2 * math.sin(math.pi * m / 4), "*", color=colors[m], ms=14, zorder=5)
    for xv in (-1, 1):
        axL.axvline(xv, color="#aaaaaa", ls="-", lw=1.0)
    axL.axvline(0, color="#888888", ls="--", lw=0.8)
    axL.annotate(r"$t=-1$: $-(x^2-y^2)^2$" "\n" r"$a_m\sim(1{+}t)^{-1/2}$",
                 (-0.98, 9), fontsize=8, color="#333", ha="left")
    axL.annotate(r"$t=+1$: $8x^2y^2$" "\n" r"$a_m\sim(1{-}t)^{\gamma_+}$",
                 (0.98, 0.62), fontsize=8, color="#333", ha="right")
    axL.set_xlabel(r"anharmonicity  $t$"); axL.set_ylabel(r"admittance  $a_m$ (log)")
    axL.set_title("Full range: divergence at $t=-1$, slow decay at $t=+1$",
                  fontweight="bold", fontsize=10)
    axL.set_xlim(-1.05, 1.05); axL.grid(True, which="both", color="#eee")
    axL.spines[["top", "right"]].set_visible(False); axL.legend(frameon=False)

    # right panel: log-log endpoint fits (the branch points that fix R=1)
    for m in (1, 2):
        # t -> -1 branch: log a_m vs log(1+t)
        d = np.array([1 + t for t in tm]); y = np.array([A[m][t] for t in tm])
        gm, Cm = gminuses[m]
        axR.plot(np.log(d), np.log(y), "o", color=colors[m], ms=6,
                 label=fr"$m={m}$, $t\to-1$ (slope {gm:.2f})")
        xs = np.log(np.array([d.min(), d.max()]))
        axR.plot(xs, math.log(Cm) + gm * xs, "-", color=colors[m], lw=1.3)
        # t -> +1 branch: log a_m vs log(1-t)
        d2 = np.array([1 - t for t in tp]); y2 = np.array([A[m][t] for t in tp])
        gp, Cp = gpluses[m]
        axR.plot(np.log(d2), np.log(y2), "s", color=colors[m], ms=6, mfc="white",
                 label=fr"$m={m}$, $t\to+1$ (slope {gp:.2f})")
        xs2 = np.log(np.array([d2.min(), d2.max()]))
        axR.plot(xs2, math.log(Cp) + gp * xs2, ":", color=colors[m], lw=1.6)
    axR.set_xlabel(r"$\log(1\mp t)$  (distance to endpoint)")
    axR.set_ylabel(r"$\log a_m$")
    axR.set_title(r"Endpoint branch points at $t=\pm1$ $\Rightarrow$ $R=1$",
                  fontweight="bold", fontsize=10)
    axR.grid(True, color="#eee")
    axR.spines[["top", "right"]].set_visible(False)
    axR.legend(frameon=False, fontsize=7.5)

    fig.tight_layout()
    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_PATH, dpi=150)
    print(f"\nSaved {FIGURE_PATH}")


if __name__ == "__main__":
    main()
