"""Closed-form investigation of the Case III capacity constant J(a).

For the pure homogeneous gate germ F = -x^4 + a x^2 y^2 + y^4 the capacity is
scale invariant: cap_T = J(a) T for all T, so cap at T = 1 equals J(a). We solve
this directly (no confinement) with reservoirs deep in the two +-x valleys.

Analytic results (verified below):

  * J(0) = 1 exactly. At a = 0 the weight e^{x^4 - y^4} separates, the committor
    is one-dimensional h(x), and
        J(0) = [ ∫ e^{-x^4} dx ]^{-1} · ∫ e^{-y^4} dy = 1.

  * First-order slope (envelope theorem; h_0'(x) = e^{-x^4}/∫e^{-x^4}):
        dJ/da |_0 = - ( ∫ x^2 e^{-x^4} dx / ∫ e^{-x^4} dx )^2
                  = - ( Γ(3/4)/Γ(1/4) )^2  =  - 2π^2 / Γ(1/4)^4  ≈  -0.11424.

  * The second-order coefficient already needs the resolvent of the a=0 operator
    (a genuine 2-D PDE solve), so J(a) has no elementary closed form beyond first
    order — the low-order data are elementary, the full invariant is not. This is
    the precise sense in which J(a) is a PDE-algebraic invariant.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma

from branching_gate_committor import solve_committor

FIGURE_PATH = Path("figures/nonsplittable_Ja.png")
SLOPE0 = -(gamma(0.75) / gamma(0.25)) ** 2   # = -2 pi^2 / Gamma(1/4)^4


def J_num(a, L=2.0, n=451, rres=0.22):
    """Scale-invariant J(a) = cap at T=1 for the pure germ, reservoirs in valleys."""
    F = lambda x, y: -x ** 4 + a * x ** 2 * y ** 2 + y ** 4
    xr = L - rres - 0.1
    _, _, _, cap = solve_committor(F, L, n, 1.0, [(xr, 0.0)], [(-xr, 0.0)], rres)
    return cap


def main() -> None:
    print(f"Analytic:  J(0) = 1,   dJ/da|_0 = -(G(3/4)/G(1/4))^2 = -2pi^2/G(1/4)^4"
          f" = {SLOPE0:.5f}")

    print("\nJ(0) domain-independence (-> 1):")
    for L in [1.6, 1.8, 2.0, 2.2]:
        print(f"  L={L}: J(0) = {J_num(0.0, L):.5f}")

    j0, j05, j10 = J_num(0.0), J_num(0.05), J_num(0.10)
    tangent = (4 * j05 - 3 * j0 - j10) / 0.1   # O(a^2)-corrected slope at 0
    print(f"\nTangent slope at a=0 (numeric) = {tangent:.5f}   "
          f"(predicted {SLOPE0:.5f})")

    a_list = np.arange(0.0, 3.01, 0.25)
    Js = np.array([J_num(a) for a in a_list])
    print("\nFull curve J(a):")
    for a, J in zip(a_list, Js):
        print(f"  a={a:.2f}: J={J:.4f}")
    sm = a_list <= 1.0
    c2, c1, c0 = np.polyfit(a_list[sm], Js[sm], 2)
    print(f"\nsmall-a fit: J ≈ {c0:.4f} + ({c1:.4f}) a + ({c2:.4f}) a^2")
    print(f"  (c0=1, c1=-2pi^2/Γ(1/4)^4={SLOPE0:.4f} exact; c2 needs the 2-D resolvent)")

    # figure: J(a) with the exact tangent 1 + SLOPE0*a
    fig, ax = plt.subplots(figsize=(6.8, 4.9))
    aa = np.linspace(0, 3, 100)
    ax.plot(aa, 1 + SLOPE0 * aa, color="#888888", ls="--", lw=1.4,
            label=r"exact tangent  $1-\frac{2\pi^2}{\Gamma(1/4)^4}a$")
    ax.plot(a_list, Js, "o-", color="#0072B2", ms=6, label="J(a)  (pure germ)")
    ax.plot(0, 1, "s", color="#D55E00", ms=9, label="J(0)=1 (exact)")
    ax.set_xlabel(r"coupling  $a$  in  $-x^4+a\,x^2y^2+y^4$")
    ax.set_ylabel(r"capacity constant  $J(a)$")
    ax.set_title(r"Case III $J(a)$: $J(0)=1$, $J'(0)=-2\pi^2/\Gamma(1/4)^4$ exact",
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
