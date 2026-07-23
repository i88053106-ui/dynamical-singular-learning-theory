"""Split-gate capacity classification via transverse Newton pole data.

For a split gate with local form

    F - H = -(mu/2) s^2 + K(v) + higher order        (s: unstable, v: transverse)

the derived capacity formula (docs/capacity-derivation.md)

    cap_T = e^{-H/T} I_K(T) sqrt(mu T / 2pi),  I_K(T) = int e^{-K/T} dv,

reduces the capacity prefactor to the transverse Laplace integral I_K. If

    I_K(T) ~ C T^{lambda_perp} (log 1/T)^{m_perp - 1},

then

    cap_T ~ e^{-H/T} C sqrt(mu/2pi) T^{lambda_perp + 1/2} (log 1/T)^{m_perp - 1},

i.e. the capacity asymptotic data are

    (kappa, r) = (lambda_perp + 1/2, m_perp).

For a Newton-nondegenerate germ K in two transverse variables the pole data
(lambda_perp, m_perp) are read off the Newton polyhedron: with d the Newton
distance (where the diagonal meets the boundary),

    lambda_perp = 1/d,
    m_perp = 2 if the diagonal meets a vertex (0-dim face), else 1 (edge).

This script computes the Newton prediction from each germ's monomial support
and independently verifies the Laplace scaling numerically. The decisive
contrast is y^4 + y^2 z^2 + z^4 (middle term collinear -> no log, m=1) versus
y^6 + y^2 z^2 + z^6 (middle term a vertex -> log, m=2).
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import integrate

OUTPUT_PATH = Path("data/split_gate_classification.csv")


# --------------------------------------------------------------------------
# Newton pole data from monomial support
# --------------------------------------------------------------------------
def _cross(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def newton_pole_data(support):
    """Return (lambda_perp, m_perp, d) for a two-variable germ.

    `support` is the list of exponent pairs of the monomials of K.
    """
    pts = sorted(set(support))
    # Non-dominated points toward the origin (candidate compact-face vertices).
    pareto = [p for p in pts
              if not any(q != p and q[0] <= p[0] and q[1] <= p[1] for q in pts)]
    pareto.sort()

    # Strict lower hull (drops collinear points).
    hull = []
    for p in pareto:
        while len(hull) >= 2 and _cross(hull[-2], hull[-1], p) <= 0:
            hull.pop()
        hull.append(p)

    if len(hull) == 1:
        d = hull[0][0]
        return 1.0 / d, 2, d

    for i in range(len(hull) - 1):
        p, q = hull[i], hull[i + 1]
        denom = (q[0] - p[0]) - (q[1] - p[1])
        if abs(denom) < 1e-12:
            continue
        sigma = (p[1] - p[0]) / denom
        if -1e-9 <= sigma <= 1 + 1e-9:
            d = p[0] + sigma * (q[0] - p[0])
            at_vertex = sigma < 1e-7 or sigma > 1 - 1e-7
            return 1.0 / d, (2 if at_vertex else 1), d

    raise ValueError(f"diagonal intersection not found for support {support}")


# --------------------------------------------------------------------------
# Independent numerical transverse Laplace integral
# --------------------------------------------------------------------------
def I_K(T, K, lo=-34.0, hi=3.0, epsrel=1e-9, epsabs=1e-14):
    """int_{R^2} exp(-K(y,z)/T) dy dz for K even in y and z (fourfold sym)."""
    def integrand(v, u):
        return 4.0 * math.exp(u + v) * math.exp(-K(math.exp(u), math.exp(v)) / T)

    breaks = [lo, -12.0, -8.0, -5.0, -3.0, -1.5, 0.0, hi]
    total = 0.0
    for i in range(len(breaks) - 1):
        for j in range(len(breaks) - 1):
            val, _ = integrate.dblquad(
                integrand, breaks[i], breaks[i + 1], breaks[j], breaks[j + 1],
                epsrel=epsrel, epsabs=epsabs,
            )
            total += val
    return total


# --------------------------------------------------------------------------
# Germ library
# --------------------------------------------------------------------------
GERMS = [
    ("y^2 + z^2",            lambda y, z: y**2 + z**2,                 [(2, 0), (0, 2)]),
    ("y^4 + z^2",            lambda y, z: y**4 + z**2,                 [(4, 0), (0, 2)]),
    ("y^4 + z^4",            lambda y, z: y**4 + z**4,                 [(4, 0), (0, 4)]),
    ("y^6 + z^6",            lambda y, z: y**6 + z**6,                 [(6, 0), (0, 6)]),
    ("y^4 + y^2 z^2 + z^4",  lambda y, z: y**4 + y**2 * z**2 + z**4,   [(4, 0), (2, 2), (0, 4)]),
    ("y^6 + y^2 z^2 + z^6",  lambda y, z: y**6 + y**2 * z**2 + z**6,   [(6, 0), (2, 2), (0, 6)]),
]

TEMPS = [1e-2, 1e-3, 1e-4]


def main() -> None:
    rows = []
    print(f"{'germ':<22}{'Newton (lam,m)':<16}{'(kappa,r)':<12}"
          f"{'r(T) = I_K / [T^lam (log1/T)^(m-1)]'}")
    print("-" * 92)

    for name, K, support in GERMS:
        lam, m, d = newton_pole_data(support)
        kappa, r_cap = lam + 0.5, m

        ratios = []
        for T in TEMPS:
            val = I_K(T, K)
            denom = T ** lam * (math.log(1.0 / T)) ** (m - 1)
            ratios.append(val / denom)
            rows.append(dict(germ=name, T=T, I_K=val,
                             lambda_perp=lam, m_perp=m,
                             kappa=kappa, r=r_cap, ratio=val / denom))

        ratio_str = "  ".join(f"{x:9.5f}" for x in ratios)
        print(f"{name:<22}({lam:.3f}, {m})      "
              f"({kappa:.3f}, {r_cap})   {ratio_str}")

    df = pd.DataFrame(rows)
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_PATH, index=False)

    print("-" * 92)
    print("r(T) -> const confirms the predicted (lambda_perp, m_perp).")
    print("Decisive contrast: quartic (m=1) has r(T) flat; sextic (m=2) needs the")
    print("log factor -- without it I_K/T^(1/2) would grow like log(1/T).")
    print(f"Saved {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
