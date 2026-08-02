"""General classification of the non-splittable-gate capacity constant J.

For a homogeneous degree-d gate germ F = -x^d + y^d + (coupling), the
scale-invariant capacity J is governed by:

  * Separable backbone.  For the separable germ -U(x) + V(y) the committor is
    1-D and J0 = [∫ e^{-V}] / [∫ e^{-U}]; for U = V = (·)^d this is 1.

  * General first-order coupling law.  Adding W(x,y) = Σ c_pq x^p y^q,
        J'(0)/J0 = - Σ c_pq <x^p>_d <y^q>_d,
        <x^p>_d = ∫ x^p e^{-x^d} / ∫ e^{-x^d}
                = Γ((p+1)/d)/Γ(1/d)   (p even),   0   (p odd).
    Derivation: dJ/dε = -(1/Z_U^2) ∫∫ W e^{-U-V}, with h0'(x)=e^{-U}/Z_U.

  * Parity selection rule: any coupling with an odd exponent has zero first-order
    effect on J.

  * Beyond first order the coefficient requires the resolvent of the a=0 operator
    (a 2-D PDE solve), so J is non-elementary for every non-splittable germ; only
    the separable backbone and the first-order Gamma-moment response are elementary.

This script verifies the law across degrees 4, 6, 8 and several couplings.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma

from branching_gate_committor import solve_committor

FIGURE_PATH = Path("figures/nonsplittable_classification.png")


def moment(p, d):
    """<x^p>_d = Γ((p+1)/d)/Γ(1/d) for even p, 0 for odd p."""
    return 0.0 if p % 2 else gamma((p + 1) / d) / gamma(1.0 / d)


def J_num(d, coupling, L, n=451, rres=0.22):
    """Scale-invariant J at T=1 for F = -x^d + y^d + Σ c x^p y^q."""
    def F(x, y):
        val = -np.abs(x) ** d + np.abs(y) ** d
        for (p, q, c) in coupling:
            val = val + c * x ** p * y ** q
        return val
    xr = L - rres - 0.1
    _, _, _, cap = solve_committor(F, L, n, 1.0, [(xr, 0.0)], [(-xr, 0.0)], rres)
    return cap


def slope(d, p, q, L, eps=0.08):
    Jp = J_num(d, [(p, q, +eps)], L)
    Jm = J_num(d, [(p, q, -eps)], L)
    return (Jp - Jm) / (2 * eps)


def main() -> None:
    print("Separable backbone J0 = 1 for -x^d + y^d:")
    for d, L in [(4, 1.9), (6, 1.7), (8, 1.5)]:
        print(f"  d={d}: J0 = {J_num(d, [], L):.4f}")

    # (d, p, q, L) test cases; predicted slope = -<x^p>_d <y^q>_d
    cases = [
        (4, 2, 2, 1.9),
        (6, 4, 2, 1.7), (6, 2, 4, 1.7), (6, 3, 3, 1.7),
        (8, 6, 2, 1.5), (8, 4, 4, 1.5), (8, 2, 6, 1.5), (8, 5, 3, 1.5),
    ]
    print("\nFirst-order coupling law  J'(0) = -<x^p>_d <y^q>_d:")
    print(f"{'germ':>26} {'numeric':>10} {'predicted':>10}")
    num, pred = [], []
    for (d, p, q, L) in cases:
        s = slope(d, p, q, L)
        pr = -moment(p, d) * moment(q, d)
        num.append(s); pred.append(pr)
        tag = f"-x^{d}+y^{d}+eps x^{p}y^{q}"
        print(f"{tag:>26} {s:>10.5f} {pr:>10.5f}")

    # figure: numeric vs predicted (should lie on the diagonal)
    fig, ax = plt.subplots(figsize=(6.2, 6.0))
    lo = min(min(num), min(pred)) - 0.01
    ax.plot([lo, 0.005], [lo, 0.005], color="#888888", ls="--", lw=1.2,
            label="numeric = predicted")
    for (d, p, q, L), s, pr in zip(cases, num, pred):
        ax.plot(pr, s, "o", ms=8,
                color={4: "#0072B2", 6: "#D55E00", 8: "#009E73"}[d])
    # legend by degree
    for d, col in [(4, "#0072B2"), (6, "#D55E00"), (8, "#009E73")]:
        ax.plot([], [], "o", color=col, label=f"degree {d}")
    ax.set_xlabel(r"predicted  $-\langle x^p\rangle_d\langle y^q\rangle_d$")
    ax.set_ylabel(r"numeric  $J'(0)$")
    ax.set_title("Non-splittable J: first-order coupling law", fontweight="bold")
    ax.set_aspect("equal")
    ax.grid(True, color="#eaeaea")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    fig.tight_layout()
    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_PATH, dpi=150)
    print(f"\nParity selection rule: odd couplings (x^3y^3, x^5y^3) give J'(0)=0.")
    print(f"Saved {FIGURE_PATH}")


if __name__ == "__main__":
    main()
