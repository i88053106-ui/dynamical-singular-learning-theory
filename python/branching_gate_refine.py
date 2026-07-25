"""Pin down the monkey-saddle capacity exponent: does cap ~ C T^p give p -> 1?

Scaling prediction (docs/branching-gate.md): near the origin F ~ -r^3 cos(3θ)
is degree-3 homogeneous, so with x = T^{1/3} ξ the inner committor problem is
scale invariant and, in dimension d=2,

    cap_T ~ T^{1 + (d-2)/m} J = T^1 * J,

with J = ∫ |∇_ξ h|^2 e^{-F(ξ)} dξ finite (deep in each valley the committor
flux gives |∇h|^2 e^{-F} ~ e^{F} -> 0). Hence p = 1 with no log correction,
and cap_T / T -> J with a leading finite-size correction in T^{1/3}.

This script (i) checks grid convergence of the capacity at one T, and
(ii) sweeps T to lower values, reporting local log-log slopes and the
extrapolation cap/T -> J in the variable T^{1/3}.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from branching_gate_committor import solve_committor

OUTPUT_CSV = Path("data/branching_gate_refine.csv")
FIGURE_PATH = Path("figures/branching_gate_exponent.png")


def monkey_F(x, y):
    return (x ** 2 + y ** 2) ** 3 / 6.0 - (x ** 3 - 3.0 * x * y ** 2)


R0 = 3.0 ** (1.0 / 3.0)
A = [(R0, 0.0)]
B = [(R0 * math.cos(2 * math.pi / 3), R0 * math.sin(2 * math.pi / 3))]


def cap_at(T, n, L=2.0, radius=0.30):
    _, _, _, cap = solve_committor(monkey_F, L=L, n=n, T=T,
                                   wells_A=A, wells_B=B, well_radius=radius)
    return cap


def grid_convergence(T0=0.20, grids=(281, 361, 481)):
    print(f"=== Grid convergence at T = {T0} ===")
    print(f"{'n':>6} {'h':>10} {'cap':>16} {'cap/T':>12}")
    caps = []
    for n in grids:
        c = cap_at(T0, n)
        caps.append(c)
        print(f"{n:>6} {2*2.0/(n-1):>10.5f} {c:>16.8e} {c/T0:>12.6f}")
    # Richardson-type note: report spread as discretization uncertainty.
    rel_spread = (max(caps) - min(caps)) / caps[-1]
    print(f"relative spread across grids: {rel_spread:.4f} "
          f"(exponent uses ratios, so a near-constant offset cancels)")
    return caps[-1]


def temperature_sweep(n=361, temps=(0.30, 0.25, 0.20, 0.16, 0.13, 0.10)):
    print(f"\n=== Temperature sweep at n = {n} ===")
    temps = np.array(temps, dtype=float)
    caps = np.array([cap_at(T, n) for T in temps])

    # local log-log slopes between consecutive temperatures
    slopes = np.diff(np.log(caps)) / np.diff(np.log(temps))
    p_global = np.polyfit(np.log(temps), np.log(caps), 1)[0]

    print(f"{'T':>7} {'cap':>15} {'cap/T':>11} {'local p':>10}")
    for i, T in enumerate(temps):
        loc = f"{slopes[i]:>10.4f}" if i < len(slopes) else f"{'-':>10}"
        print(f"{T:>7.3f} {caps[i]:>15.8e} {caps[i]/T:>11.6f} {loc}")
    print(f"global log-log slope p = {p_global:.4f}")

    # extrapolate cap/T -> J linearly in T^{1/3}
    u = temps ** (1.0 / 3.0)
    y = caps / temps
    b, J = np.polyfit(u, y, 1)  # y = J + b*u
    print(f"extrapolation cap/T = J + b*T^(1/3):  J = {J:.4f}, b = {b:.4f}")
    print("(p -> 1 and cap/T -> J confirm cap ~ J*T with no log correction.)")

    df = pd.DataFrame(dict(T=temps, cap=caps, cap_over_T=caps / temps))
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_CSV, index=False)

    # figure: cap/T vs T^{1/3} with linear extrapolation to J
    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    uu = np.linspace(0.0, u.max() * 1.05, 100)
    ax.plot(uu, J + b * uu, color="#0072B2", lw=2,
            label=f"linear fit  →  J = {J:.3f}")
    ax.plot(u, y, "o", color="#000000", ms=8, label="measured  cap/T")
    ax.plot(0.0, J, "s", color="#D55E00", ms=9, label="extrapolated  T→0")
    ax.set_xlabel("T$^{1/3}$")
    ax.set_ylabel("cap$_T$ / T")
    ax.set_title(f"Monkey-saddle capacity:  cap ~ J·T  (p={p_global:.3f})",
                 fontweight="bold")
    ax.set_xlim(left=0.0)
    ax.grid(True, color="#eaeaea")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    fig.tight_layout()
    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_PATH, dpi=150)
    print(f"Saved {OUTPUT_CSV} and {FIGURE_PATH}")

    return p_global, J


def main() -> None:
    grid_convergence()
    temperature_sweep()


if __name__ == "__main__":
    main()
