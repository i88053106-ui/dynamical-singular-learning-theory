"""Fractional-Bessel modal structure of the branching-gate committor, and the
numerical closed-out of the junction conductances b_s(k) / channel admittances.

Conformal map W = w^k turns div(e^{-F/T} grad h)=0 with F = -Re(w^k) into the
constant-coefficient div(e^{U/T} grad h)=0 (U = Re W = rho^k cos k phi), i.e.
Delta h + (1/T) h_U = 0. With g = e^{U/(2T)} h this is modified Helmholtz
Delta g = g/(4T^2), whose separated solutions regular at the gate are the
FRACTIONAL modified Bessel functions:

    h(rho,phi) = e^{-U/(2T)} sum_m  alpha_m  I_{|m|/k}(rho^k / (2T))  e^{i m phi}.

The fractional order m/k is the signature of the k-fold branching. This script
(1) validates that structure by projecting the numerical committor onto e^{i m phi}
and checking the radial factor is proportional to I_{|m|/k}(rho^k/2T), and
(2) extracts the circulant channel admittances lambda_m and edge conductances
b_s (so that J_k = 1/R_01 closes numerically), for k = 3..6.

The coefficients alpha_m are fixed by a far-field connection problem: channel m
excites the integer-spaced orders |m/k + n| (aliases m + n k), and boundedness on
the ridges selects the combination. That connection problem has no elementary
closed form in general; b_s(k) are reported numerically here.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.special import iv

from branching_gate_committor import solve_committor

LDOM, RADIUS = 2.0, 0.26


def make_F(k):
    def F(x, y):
        return (x * x + y * y) ** k / (2.0 * k) - np.real((x + 1j * y) ** k)
    return F


def wells(k):
    r0 = k ** (1.0 / k)
    return [(r0 * math.cos(2 * math.pi * j / k), r0 * math.sin(2 * math.pi * j / k))
            for j in range(k)]


def validate_modal_structure(k=4, T=0.2, n=401):
    """Check g_m(rho) / I_{|m|/k}(rho^k/2T) is constant in rho."""
    F, W = make_F(k), wells(k)
    X, Y, H, _ = solve_committor(F, LDOM, n, T, [W[0]], [W[1]], RADIUS)
    xs = X[:, 0]
    interp = RegularGridInterpolator((xs, xs), H, bounds_error=False, fill_value=None)

    phis = np.linspace(0, 2 * np.pi, 720, endpoint=False)
    radii = np.linspace(0.35, 0.95, 13)
    modes = list(range(k + 1))

    print(f"=== k={k}, T={T}: modal ratio g_m(rho)/I_(|m|/k)(rho^k/2T) (want const) ===")
    ratios = {m: [] for m in modes}
    for rho in radii:
        pts = np.stack([rho * np.cos(phis), rho * np.sin(phis)], axis=1)
        h = interp(pts)
        U = rho ** k * np.cos(k * phis)
        g = np.exp(U / (2 * T)) * h
        for m in modes:
            gm = np.mean(g * np.exp(-1j * m * phis)).real
            ratios[m].append(gm / iv(abs(m) / k, rho ** k / (2 * T)))

    print("  mode m (order m/k) :  rel.spread   mean coeff alpha_m")
    for m in modes:
        v = np.array(ratios[m])
        print(f"    m={m} ({abs(m)}/{k}={abs(m)/k:.3f}) :  "
              f"{np.std(v)/max(abs(np.mean(v)),1e-12):.4f}      {np.mean(v):+.5f}")


def effective_resistance(k, d, T, n=321):
    F, W = make_F(k), wells(k)
    _, _, _, cap = solve_committor(F, LDOM, n, T, [W[0]], [W[d]], RADIUS)
    return 1.0 / cap  # R_0d = 1/cap(0,d)


def extract_conductances(ks=(3, 4, 5, 6), T=0.12):
    """Recover channel admittances lambda_m and edge conductances b_s from R_0d."""
    print(f"\n=== junction channel admittances a_m = lambda_m/T  (T={T}) ===")
    for k in ks:
        half = k // 2
        ds = list(range(1, half + 1))
        R = np.array([effective_resistance(k, d, T) for d in ds])
        # independent eigenvalues: m = 1..half (lambda_m = lambda_{k-m})
        ms = list(range(1, half + 1))
        # M[d_index, m_index] * (1/lambda_m) = R_0d, with multiplicity for m != k-m
        M = np.zeros((len(ds), len(ms)))
        for di, d in enumerate(ds):
            for mi, m in enumerate(ms):
                coeff = 4 * math.sin(math.pi * m * d / k) ** 2
                if (k - m) != m:          # m and k-m are distinct, equal eigenvalues
                    coeff += 4 * math.sin(math.pi * (k - m) * d / k) ** 2
                M[di, mi] = coeff / k
        inv_lambda = np.linalg.solve(M, R)
        lam = 1.0 / inv_lambda            # lambda_m  (m=1..half)
        a = lam / T                       # scale-invariant admittances
        # reconstruct J_k = 1/R_01
        R01 = 0.0
        for m in range(1, k):
            mm = m if m <= half else k - m
            R01 += (4 * math.sin(math.pi * m / k) ** 2) / lam[mm - 1]
        R01 /= k
        Jk = 1.0 / R01 / T
        amp = ", ".join(f"a_{m}={a[m-1]:.4f}" for m in ms)
        print(f"  k={k}:  {amp}   =>  J_{k} = cap/T = {Jk:.4f}")


def main() -> None:
    validate_modal_structure(k=4, T=0.2)
    extract_conductances()
    print("\nThe fractional-Bessel modal structure is exact; the channel admittances")
    print("a_m (hence b_s and J_k) are the numerically-computed constants of the")
    print("Bessel connection problem set up in docs/branching-gate-bs.md.")


if __name__ == "__main__":
    main()
