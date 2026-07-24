"""Visualize the Langevin exit-time verification of the capacity clock law.

Produces a three-panel figure summarizing the numerical evidence recorded in
research-log/2026-07.md:

  A. Measured mean first-passage time E[tau] vs T, bracketed by the leading
     and refined analytic predictions (shared exp(1/4T) growth).
  B. Crossover: measured/leading (< 1) and measured/refined (> 1) bracket 1,
     the expected behaviour since log(1/T) ~ 3 is comparable to gamma+6log2.
  C. The Newton-degenerate signature: the reduced quantity that ADDS
     -(-log log(1/T)) is flat (range 0.08), while the plain power-prefactor
     reduction drifts (range 0.29) -- the data prefer the log-corrected form.

Consolidated measurements are the authoritative values from the runs logged
on 2026-07-23 (vectorized Euler-Maruyama, dt=1e-3, N=600-800).
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

FIGURE_PATH = Path("figures/exit_time_verification.png")

# CVD-safe categorical hues (Okabe-Ito) + text tokens.
C_MEASURED = "#000000"
C_LEADING = "#0072B2"   # blue
C_REFINED = "#E69F00"   # orange
C_SIG = "#009E73"       # bluish green
C_POWER = "#D55E00"     # vermillion
INK = "#222222"
MUTED = "#888888"

# Consolidated measured separatrix first-passage times (dt=1e-3).
T = np.array([0.080, 0.070, 0.060, 0.055, 0.050, 0.045, 0.040])
TAU = np.array([20.476, 27.653, 43.567, 61.086, 87.216, 131.089, 250.018])
SE = np.array([0.800, 1.061, 1.429, 2.456, 2.900, 4.712, 8.486])

EK_CONSTANT = 3.0 * math.sqrt(2.0) * math.pi ** 1.5
G6 = np.euler_gamma + 6.0 * math.log(2.0)


def separatrix_leading(t):
    return 0.5 * EK_CONSTANT * np.sqrt(t) / np.log(1.0 / t) * np.exp(1.0 / (4.0 * t))


def separatrix_refined(t):
    return separatrix_leading(t) * np.log(1.0 / t) / (np.log(1.0 / t) + G6)


def main() -> None:
    lead = separatrix_leading(T)
    refined = separatrix_refined(T)

    g = np.log(TAU) - 1.0 / (4.0 * T) - 0.5 * np.log(T)
    with_log = g + np.log(np.log(1.0 / T))     # -> const if 1/log(1/T) present
    without = g                                 # -> const if plain power prefactor

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(15.0, 4.6))
    for ax in (axA, axB, axC):
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(colors=INK)
        ax.grid(True, color="#eaeaea", linewidth=0.8)
        ax.set_axisbelow(True)

    # --- Panel A: E[tau] vs T with prediction curves ---
    tt = np.linspace(T.min(), T.max(), 200)
    axA.semilogy(tt, separatrix_leading(tt), color=C_LEADING, lw=2,
                 label="leading prediction")
    axA.semilogy(tt, separatrix_refined(tt), color=C_REFINED, lw=2,
                 label="refined prediction")
    axA.errorbar(T, TAU, yerr=SE, fmt="o", color=C_MEASURED, ms=6,
                 capsize=3, lw=1.5, label="measured (Langevin)", zorder=5)
    axA.set_xlabel("temperature  T", color=INK)
    axA.set_ylabel("mean first-passage time  E[τ]", color=INK)
    axA.set_title("A. Exit time vs prediction", color=INK, fontweight="bold")
    axA.invert_xaxis()
    axA.legend(frameon=False, fontsize=9, loc="upper left")

    # --- Panel B: crossover ratios bracket 1 ---
    axB.axhline(1.0, color=MUTED, lw=1.5, ls="--")
    axB.text(T.min(), 1.02, "measured", color=MUTED, fontsize=9, va="bottom")
    axB.plot(T, TAU / lead, "o-", color=C_LEADING, lw=2, ms=6,
             label="measured / leading")
    axB.plot(T, TAU / refined, "s-", color=C_REFINED, lw=2, ms=6,
             label="measured / refined")
    axB.set_xlabel("temperature  T", color=INK)
    axB.set_ylabel("ratio to prediction", color=INK)
    axB.set_title("B. Crossover brackets the constant", color=INK, fontweight="bold")
    axB.invert_xaxis()
    axB.legend(frameon=False, fontsize=9, loc="center left")

    # --- Panel C: signature (flat) vs power-only (drifts), each minus its mean ---
    axC.axhline(0.0, color=MUTED, lw=1.0)
    axC.plot(T, with_log - with_log.mean(), "o-", color=C_SIG, lw=2, ms=6,
             label=f"with  −log log(1/T)  (range {np.ptp(with_log):.2f})")
    axC.plot(T, without - without.mean(), "s--", color=C_POWER, lw=2, ms=6,
             label=f"power prefactor only  (range {np.ptp(without):.2f})")
    axC.set_xlabel("temperature  T", color=INK)
    axC.set_ylabel("reduced quantity  − its mean", color=INK)
    axC.set_title("C. Newton-degenerate signature", color=INK, fontweight="bold")
    axC.invert_xaxis()
    axC.legend(frameon=False, fontsize=9, loc="upper center")

    fig.suptitle(
        "Langevin verification of the split-gate clock law   "
        "E[τ] ~ C · T^{1/2} / log(1/T) · e^{1/(4T)}",
        fontsize=12, fontweight="bold", color=INK,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_PATH, dpi=150)
    print(f"Saved {FIGURE_PATH}")


if __name__ == "__main__":
    main()
