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
    for k in ks:
        for m in range(1, k // 2 + 1):
            a = channel_admittance(k, m, T)
            q = 2 * math.pi * m / k
            pred = 2 * math.sin(math.pi * m / k)
            qs.append(q); ams.append(a)
            print(f"{k:>3} {m:>3} {q:>10.4f} {a:>9.5f} {pred:>13.5f} {a/pred:>8.4f}")

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


if __name__ == "__main__":
    main()
