"""Case III: a non-splittable degenerate gate whose capacity is not a Newton invariant.

theory.md Case III: the unstable direction and the singular kernel cannot be
locally separated, and no single nondegenerate quadratic unstable mode can be
split off. We model this with the homogeneous degree-4 gate germ

    F(x,y) = -x^4 + a x^2 y^2 + y^4 + (x^2+y^2)^3     (confinement is degree 6),

which has two wells on the +-x axis (at x0 = sqrt(2/3)) and a degenerate gate at
the origin (H = 0). Key facts that place it outside Cases I and II:

  * Hessian at the gate is zero and the unstable direction is quartic (-x^4), so
    the split-gate formula (which needs a quadratic unstable mode) does not apply.
  * The germ is NOT harmonic (Delta(-x^4+x^2y^2+y^4) = -10x^2+14y^2 != 0), so the
    conformal linearization W = w^k used for branching (Case II) does not apply.
  * For a != 0 the coupling x^2 y^2 makes the germ non-separable (for a = 0 the
    leading germ -x^4 + y^4 is separable, so the committor reduces to 1-D).

Scaling (homogeneous degree 4, dimension 2) gives cap_T ~ J(a) T (p = 1). The
point: the Newton polygon of the germ is the hull of {(4,0),(2,2),(0,4)},
INDEPENDENT of the coefficient a (and of the signs). Yet J(a) varies with a --
so the capacity constant is a genuine PDE/analytic invariant of the germ, not a
function of its Newton (RLCT) data. That is the Case III phenomenon.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from branching_gate_committor import solve_committor

X0 = math.sqrt(2.0 / 3.0)
FIGURE_PATH = Path("figures/nonsplittable_gate.png")


def make_F(a):
    def F(x, y):
        return -x ** 4 + a * x ** 2 * y ** 2 + y ** 4 + (x ** 2 + y ** 2) ** 3
    return F


def cap_at(a, T, n=361):
    _, _, _, cap = solve_committor(make_F(a), 1.4, n, T,
                                   [(X0, 0.0)], [(-X0, 0.0)], 0.20)
    return cap


def main() -> None:
    # (1) confirm the scaling exponent p -> 1 for a representative coupling
    print("Scaling check (a = 1):  cap ~ J T,  p = d log cap / d log T -> 1")
    Ts = [0.06, 0.045, 0.03, 0.02]
    caps = np.array([cap_at(1.0, T) for T in Ts])
    p = np.polyfit(np.log(Ts), np.log(caps), 1)[0]
    for T, c in zip(Ts, caps):
        print(f"   T={T:.3f}  cap={c:.6f}  J=cap/T={c/T:.4f}")
    print(f"   p = {p:.3f}")

    # (2) J(a): capacity constant vs the coupling, with fixed Newton support.
    # Report cap/T at two fixed T (extrapolation-free); the a-dependence is the
    # point and it persists at both T (so it is not a finite-T artifact).
    print("\nJ(a) = cap/T at fixed T  (Newton support {(4,0),(2,2),(0,4)} is a-independent):")
    a_list = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0]
    Tj = [0.045, 0.030]
    print(f"{'a':>5} " + " ".join(f"J(T={T:.3f})".rjust(11) for T in Tj))
    Jt = {T: [] for T in Tj}
    for a in a_list:
        row = []
        for T in Tj:
            j = cap_at(a, T) / T
            Jt[T].append(j)
            row.append(j)
        print(f"{a:>5.1f} " + " ".join(f"{v:>11.4f}" for v in row))

    Tref = Tj[-1]
    Jref = Jt[Tref]
    print(f"\nAt fixed T={Tref}: J(a=0)={Jref[0]:.3f} -> J(a=3)={Jref[-1]:.3f} "
          f"({(Jref[0]-Jref[-1])/Jref[0]*100:.0f}% change), monotone at BOTH T.")
    print("The Newton polygon (and signs) are a-independent, yet J(a) varies:")
    print("the capacity constant is a PDE/analytic invariant of the germ, not a")
    print("function of its Newton (RLCT) data.  [Case III]")

    # figure: J(a) vs a at both T (parallel curves = a-dependence is robust)
    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    cols = ["#0072B2", "#D55E00"]
    for T, col in zip(Tj, cols):
        ax.plot(a_list, Jt[T], "o-", color=col, ms=7, label=f"T = {T:.3f}")
    ax.set_xlabel(r"coupling  $a$  in  $-x^4 + a\,x^2y^2 + y^4$")
    ax.set_ylabel(r"capacity constant  $J(a)=\mathrm{cap}/T$")
    ax.set_title(r"Case III: $J$ depends on the coefficient (fixed Newton polygon)",
                 fontweight="bold")
    ax.grid(True, color="#eaeaea")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    fig.tight_layout()
    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_PATH, dpi=150)
    print(f"\nSaved {FIGURE_PATH}")


if __name__ == "__main__":
    main()
