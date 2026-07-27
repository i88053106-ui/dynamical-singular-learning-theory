"""Closed form a_m = 2 sin(pi m/k) and the dispersion slope Lambda'(0+) = 1.

Large-k mechanism (docs/branching-gate-slope.md): as k -> infinity the interior
of the unit disk becomes a UNIFORM conductor (weight e^{rho^k cos k phi / T} -> 1
for rho < 1), with the deep wells acting as perfect reservoirs (contacts) on the
boundary rho = 1. The Dirichlet-to-Neumann (Steklov) eigenvalue of the unit disk
on mode e^{i m phi} is exactly |m|, so the boundary current density is |m| e^{i m phi};
integrating over one valley's angular territory [-pi/k, pi/k] gives

    a_m = \int_{-pi/k}^{pi/k} |m| cos(m phi) d phi = 2 sin(pi m / k),

hence the dispersion Lambda(q) = 2 sin(q/2) at q = 2 pi m/k, and Lambda'(0+) = 1.

This script verifies a_m -> 2 sin(pi m/k) and produces the collapse plot of the
measured admittances onto the curve 2 sin(q/2).
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from branching_admittance_asymptotics import channel_admittance

FIGURE_PATH = Path("figures/branching_steklov_dispersion.png")


def main() -> None:
    T = 0.10
    ks = [6, 8, 10, 12, 16]

    print("Verify a_m = 2 sin(pi m/k)   (ratio -> 1 as k grows):")
    print(f"{'k':>3} {'m':>3} {'q=2pi m/k':>10} {'a_m':>9} {'2sin(pi m/k)':>13} {'ratio':>8}")
    qs, ams = [], []
    adm = {k: {} for k in ks}
    for k in ks:
        for m in range(1, k // 2 + 1):
            a = channel_admittance(k, m, T)
            adm[k][m] = a
            adm[k][k - m] = a  # a_m = a_{k-m}
            q = 2 * math.pi * m / k
            pred = 2 * math.sin(math.pi * m / k)
            qs.append(q); ams.append(a)
            print(f"{k:>3} {m:>3} {q:>10.4f} {a:>9.5f} {pred:>13.5f} {a/pred:>8.4f}")

    # adjacent capacity constant J_k = 1 / [ (1/k) sum_m 4 sin^2(pi m/k) / a_m ]
    # with a_m = 2 sin(pi m/k) this is exactly (k/2) tan(pi/2k) -> pi/4.
    print("\nAdjacent constant J_k -> J_inf = pi/4:")
    print(f"{'k':>3} {'J_k(meas)':>10} {'(k/2)tan(pi/2k)':>16} {'pi/4':>8}")
    Jmeas = {}
    for k in ks:
        R01 = sum(4 * math.sin(math.pi * m / k) ** 2 / adm[k][m]
                  for m in range(1, k)) / k
        Jmeas[k] = 1.0 / R01
        ideal = (k / 2) * math.tan(math.pi / (2 * k))
        print(f"{k:>3} {Jmeas[k]:>10.5f} {ideal:>16.5f} {math.pi/4:>8.5f}")

    # figure: collapse onto Lambda(q) = 2 sin(q/2), with tangent slope 1 at 0
    fig, ax = plt.subplots(figsize=(6.8, 5.0))
    qq = np.linspace(0, math.pi, 200)
    ax.plot(qq, 2 * np.sin(qq / 2), color="#0072B2", lw=2,
            label=r"$\Lambda(q)=2\sin(q/2)$")
    ax.plot(qq, qq, color="#888888", lw=1.2, ls=":",
            label=r"tangent $\Lambda\simeq q$  (slope $\Lambda'(0^+)=1$)")
    ax.plot(qs, ams, "o", color="#D55E00", ms=7, label="measured $a_m$")
    ax.set_xlabel(r"wavenumber  $q = 2\pi m/k$")
    ax.set_ylabel(r"channel admittance  $a_m$")
    ax.set_title(r"Dispersion collapse:  $a_m = 2\sin(\pi m/k)$", fontweight="bold")
    ax.grid(True, color="#eaeaea")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    fig.tight_layout()
    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_PATH, dpi=150)
    print(f"\nSaved {FIGURE_PATH}")
    print("Lambda(q) = 2 sin(q/2)  =>  Lambda'(0+) = 1  (QED, via disk Steklov |m|).")

    # figure: J_k -> pi/4
    fig2, ax2 = plt.subplots(figsize=(6.6, 4.8))
    kk = np.linspace(3, max(ks) + 1, 200)
    ax2.axhline(math.pi / 4, color="#888888", ls="--", lw=1.2, label=r"$J_\infty=\pi/4$")
    ax2.plot(kk, (kk / 2) * np.tan(math.pi / (2 * kk)), color="#0072B2", lw=2,
             label=r"$\frac{k}{2}\tan\frac{\pi}{2k}$  (from $a_m=2\sin$)")
    ax2.plot(list(ks), [Jmeas[k] for k in ks], "o", color="#D55E00", ms=7,
             label="measured $J_k$")
    ax2.set_xlabel("valley number  k")
    ax2.set_ylabel(r"adjacent constant  $J_k=\mathrm{cap}/T$")
    ax2.set_title(r"$J_k \to J_\infty = \pi/4$", fontweight="bold")
    ax2.grid(True, color="#eaeaea")
    ax2.spines[["top", "right"]].set_visible(False)
    ax2.legend(frameon=False)
    fig2.tight_layout()
    out2 = Path("figures/branching_Jinfinity.png")
    fig2.savefig(out2, dpi=150)
    print(f"Saved {out2}")
    print(f"J_inf = pi/4 = {math.pi/4:.6f}")


if __name__ == "__main__":
    main()
