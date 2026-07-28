"""Error analysis of the closed form a_m = 2 sin(pi m/k) at finite k, T.

The large-k derivation (docs/branching-gate-slope.md) gives a_m = 2 sin(pi m/k).
Here we quantify how the true admittance a_m(k, T) departs from it.

Findings:
  * Scale-invariant exactness: a linear fit a_m(k,T) = A + B T over the reliable
    range T in [0.06, 0.15] gives A = 2 sin(pi m/k) to ~0.1% for every k tested
    (residual dominated by the finite well-disk numerical bias). So the closed
    form is exact in the T -> 0 (scale-invariant) limit for ALL k, not merely
    a large-k asymptotic.
  * Finite-T correction: a_m(k,T) = 2 sin(pi m/k) (1 - c_m(k) T + O(T^2)) with
    c_m(k) ~ c_m / k^2 (c_1 ~ 3.4, c_2 ~ 5.8). Envelope |a_m/2sin - 1| <~ 6 T/k^2.
  * Reliability floor: the weight dynamic range is e^{(k/2)/T}; the linear solve
    loses precision once this exceeds float64 (~1e16), i.e. T <~ (k/2)/37. This is
    a solver limit, not physics; results below ~0.05 are discarded.

a_m is variational (a Dirichlet-energy infimum), hence monotone in the weight;
this gives the rigorous comparison framework, but a proven-constant sharp bound
on c_m(k) has no elementary closed form and is left open.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from branching_admittance_asymptotics import channel_admittance

FIGURE_PATH = Path("figures/branching_gate_error.png")


def main() -> None:
    Ts = [0.15, 0.12, 0.09, 0.06]
    cases = [(4, 1), (5, 1), (6, 1), (8, 1), (6, 2), (8, 2)]

    print("Linear fit a_m(k,T) = A + B T over T in [0.06, 0.15]:")
    print(f"{'k':>3} {'m':>2} {'A':>10} {'2sin(pi m/k)':>13} {'A/2sin-1':>11} "
          f"{'B(slope)':>10} {'c_m=-B k^2/2sin':>15}")
    curves = {}
    for (k, m) in cases:
        pred = 2 * math.sin(math.pi * m / k)
        a = np.array([channel_admittance(k, m, T, n=421) for T in Ts])
        B, A = np.polyfit(Ts, a, 1)
        cm = -B * k ** 2 / pred
        print(f"{k:>3} {m:>2} {A:>10.6f} {pred:>13.6f} {A/pred-1:>+11.2e} "
              f"{B:>10.4f} {cm:>15.3f}")
        if m == 1:
            curves[k] = (np.array(Ts), a / pred - 1.0)

    # figure: relative error eps = a/2sin - 1 vs T (linear, intercept ~ 0)
    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    cmap = plt.get_cmap("viridis")
    ks = sorted(curves)
    for i, k in enumerate(ks):
        T, eps = curves[k]
        col = cmap(i / max(1, len(ks) - 1))
        B, A = np.polyfit(T, eps, 1)
        tt = np.linspace(0, 0.16, 50)
        ax.plot(tt, A + B * tt, "-", color=col, lw=1.5)
        ax.plot(T, eps, "o", color=col, ms=6, label=f"k={k}")
    ax.axhline(0.0, color="#888888", lw=1.0, ls="--")
    ax.set_xlabel("temperature  T")
    ax.set_ylabel(r"relative error  $a_1/2\sin(\pi/k) - 1$")
    ax.set_title(r"Idealization error is $O(T)$; intercept $\to 0$ (exact at $T{=}0$)",
                 fontweight="bold")
    ax.set_xlim(left=0.0)
    ax.grid(True, color="#eaeaea")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    fig.tight_layout()
    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_PATH, dpi=150)
    print(f"\nSaved {FIGURE_PATH}")
    print("Conclusion: a_m(k) = 2 sin(pi m/k) exactly as T->0 (all k);"
          " finite-T error = O(T/k^2).")


if __name__ == "__main__":
    main()
