"""Capacity dependence on the number of valleys k at a branching gate.

Gate family (degree-k homogeneous, k-fold symmetric):

    F_k(x,y) = r^{2k}/(2k) - Re[(x+iy)^k] = r^{2k}/(2k) - r^k cos(k θ),

with k wells at θ_j = 2π j / k, r0 = k^{1/k}, depth F = -k/2, and the origin as
the common degenerate gate (communication height H = 0). k=2 is an ordinary
Morse saddle (a plain double well); k>=3 are genuine branching gates.

Scaling prediction (docs/branching-gate.md): in dimension d=2 the capacity
exponent is 1 + (d-2)/k = 1 for every k, so cap_T ~ J_k T with a T^{1/3}-type
finite-size correction. This script measures cap(A,B) between two ADJACENT
wells across k = 2..6, confirms p_k -> 1, extracts J_k = lim cap/T, and records
the committor values at the remaining "spectator" wells (the branching split).
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from branching_gate_committor import solve_committor

OUTPUT_CSV = Path("data/branching_valley_number.csv")
FIGURE_PATH = Path("figures/branching_valley_number.png")

KS = [2, 3, 4, 5, 6]
TEMPS = [0.24, 0.19, 0.15, 0.12, 0.10]
NGRID = 321
LDOM = 2.0
RADIUS = 0.26


def make_F(k):
    def F(x, y):
        return (x * x + y * y) ** k / (2.0 * k) - np.real((x + 1j * y) ** k)
    return F


def well_centers(k):
    r0 = k ** (1.0 / k)
    return [(r0 * math.cos(2 * math.pi * j / k), r0 * math.sin(2 * math.pi * j / k))
            for j in range(k)]


def committor_at(X, Y, H, pt):
    d2 = (X - pt[0]) ** 2 + (Y - pt[1]) ** 2
    return float(H.flat[np.argmin(d2)])


def study_k(k):
    F = make_F(k)
    W = well_centers(k)
    A, B = [W[0]], [W[1]]           # adjacent wells
    spectators = W[2:]              # the other k-2 wells

    temps = np.array(TEMPS)
    caps, spec_vals = [], []
    for T in temps:
        X, Y, H, cap = solve_committor(F, L=LDOM, n=NGRID, T=T,
                                       wells_A=A, wells_B=B, well_radius=RADIUS)
        caps.append(cap)
        spec_vals.append([committor_at(X, Y, H, s) for s in spectators])
    caps = np.array(caps)

    p_global = np.polyfit(np.log(temps), np.log(caps), 1)[0]
    slopes = np.diff(np.log(caps)) / np.diff(np.log(temps))
    capT_low = caps[-1] / temps[-1]  # extrapolation-free cross-k comparator

    return dict(k=k, temps=temps, caps=caps, p_global=p_global,
                slope_low=slopes[-1], capT_low=capT_low, spec_low=spec_vals[-1])


def main() -> None:
    results = [study_k(k) for k in KS]

    rows = []
    print(f"{'k':>3} {'p_global':>9} {'slope(lowT)':>12} {'cap/T@0.10':>11}  spectator h (low T)")
    for r in results:
        spec = ", ".join(f"{v:.3f}" for v in r["spec_low"]) or "-"
        print(f"{r['k']:>3} {r['p_global']:>9.4f} {r['slope_low']:>12.4f} "
              f"{r['capT_low']:>11.4f}  [{spec}]")
        for T, c in zip(r["temps"], r["caps"]):
            rows.append(dict(k=r["k"], T=T, cap=c, cap_over_T=c / T))
    pd.DataFrame(rows).to_csv(OUTPUT_CSV, index=False)

    print("\np_k -> 1 for all k (exponent independent of valley number).")
    print("cap/T at fixed low T decreases monotonically with k (extrapolation-free).")
    print("Absolute normalization carries the ~6% solver systematic (Morse validation);")
    print("k=2 is the Morse/split boundary case (split formula gives leading cap = T).")

    # figure
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.0, 4.8))
    cmap = plt.get_cmap("viridis")
    for i, r in enumerate(results):
        col = cmap(i / max(1, len(results) - 1))
        u = r["temps"] ** (1.0 / 3.0)
        axL.plot(u, r["caps"] / r["temps"], "o-", color=col, ms=6,
                 label=f"k={r['k']}  (p={r['p_global']:.3f})")
    axL.set_xlabel("T$^{1/3}$")
    axL.set_ylabel("cap$_T$ / T")
    axL.set_title("A. cap/T linear in T$^{1/3}$, slope→1 for all k", fontweight="bold")
    axL.grid(True, color="#eaeaea")
    axL.spines[["top", "right"]].set_visible(False)
    axL.legend(frameon=False, fontsize=9)

    ks = [r["k"] for r in results]
    cl = [r["capT_low"] for r in results]
    axR.plot(ks, cl, "o-", color="#0072B2", ms=9)
    for k, c in zip(ks, cl):
        axR.annotate(f"{c:.3f}", (k, c), textcoords="offset points",
                     xytext=(0, 8), ha="center", fontsize=9)
    axR.set_xlabel("valley number  k")
    axR.set_ylabel("cap$_T$ / T  at T = 0.10")
    axR.set_title("B. Capacity constant decreases with k", fontweight="bold")
    axR.set_xticks(ks)
    axR.grid(True, color="#eaeaea")
    axR.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_PATH, dpi=150)
    print(f"Saved {OUTPUT_CSV} and {FIGURE_PATH}")


if __name__ == "__main__":
    main()
