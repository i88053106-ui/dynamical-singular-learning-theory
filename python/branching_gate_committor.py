"""Branching-gate capacity via a numerical committor solve (third target).

The split-gate classification (docs/split-gate-classification.md) assumes the
unstable direction separates from the transverse singular geometry. That
reduction breaks at a *branching gate*, where several downhill valleys meet
and the committor is genuinely two-dimensional. The canonical example is the
monkey saddle.

We solve the equilibrium potential (committor) h of the reversible generator

    T Laplacian(h) - grad F . grad h = 0,   h = 1 on A, h = 0 on B,

by a finite-volume discretization of the self-adjoint form
div(e^{-F/T} grad h) = 0, which is symmetric positive definite (no upwinding
needed). The capacity is the Dirichlet energy of the discrete minimizer,

    cap_T(A,B) = T * sum_edges a_edge (h_i - h_j)^2,  a = e^{-F/T}.

Two experiments:
  1. VALIDATION on a Morse double well, where the derived split-gate formula
     gives cap_T = (T/omega) e^{-1/(4T)} exactly. Confirms the solver.
  2. MONKEY SADDLE three-well model F = (1/6) r^6 - (x^3 - 3 x y^2): the gate
     at the origin is degenerate (Hessian = 0), so the Morse/split formula is
     inapplicable. We measure cap_T(A,B) vs T, its effective exponent, and the
     committor value at the third well C (branching signature h_C = 1/2 by
     symmetry).
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

FIGURE_PATH = Path("figures/branching_gate_committor.png")


def solve_committor(F, L, n, T, wells_A, wells_B, well_radius):
    """Finite-volume committor and capacity on [-L, L]^2.

    wells_A, wells_B: lists of (cx, cy) centres held at h = 1 and h = 0.
    Returns (X, Y, H, capacity).
    """
    xs = np.linspace(-L, L, n)
    hh = xs[1] - xs[0]
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    Fv = F(X, Y)
    a = np.exp(-Fv / T)

    fixed = np.full((n, n), np.nan)
    r2 = well_radius ** 2
    for (cx, cy) in wells_A:
        fixed[(X - cx) ** 2 + (Y - cy) ** 2 <= r2] = 1.0
    for (cx, cy) in wells_B:
        fixed[(X - cx) ** 2 + (Y - cy) ** 2 <= r2] = 0.0

    is_fixed = ~np.isnan(fixed)
    idx = -np.ones((n, n), dtype=int)
    unknown = np.argwhere(~is_fixed)
    for k, (i, j) in enumerate(unknown):
        idx[i, j] = k
    n_unk = len(unknown)

    rows, cols, vals = [], [], []
    rhs = np.zeros(n_unk)

    def edge_a(i1, j1, i2, j2):
        return 0.5 * (a[i1, j1] + a[i2, j2])

    for k, (i, j) in enumerate(unknown):
        diag = 0.0
        for (di, dj) in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            ii, jj = i + di, j + dj
            if not (0 <= ii < n and 0 <= jj < n):
                continue  # no-flux Neumann outer boundary
            ae = edge_a(i, j, ii, jj)
            diag -= ae
            if is_fixed[ii, jj]:
                rhs[k] -= ae * fixed[ii, jj]
            else:
                rows.append(k)
                cols.append(idx[ii, jj])
                vals.append(ae)
        rows.append(k)
        cols.append(k)
        vals.append(diag)

    A = csr_matrix((vals, (rows, cols)), shape=(n_unk, n_unk))
    sol = spsolve(A, rhs)

    H = fixed.copy()
    for k, (i, j) in enumerate(unknown):
        H[i, j] = sol[k]

    # Capacity = T * sum over edges of a_edge (h_i - h_j)^2.
    cap = 0.0
    dh_x = H[1:, :] - H[:-1, :]
    a_x = 0.5 * (a[1:, :] + a[:-1, :])
    cap += np.sum(a_x * dh_x ** 2)
    dh_y = H[:, 1:] - H[:, :-1]
    a_y = 0.5 * (a[:, 1:] + a[:, :-1])
    cap += np.sum(a_y * dh_y ** 2)
    cap *= T

    return X, Y, H, cap


def run_validation():
    print("=== VALIDATION: Morse double well, cap should be (T/omega) e^{-1/4T} ===")
    F = lambda x, y: (x ** 2 - 1) ** 2 / 4.0 + 0.5 * y ** 2  # omega = 1
    print(f"{'T':>6} {'cap_numeric':>14} {'cap_formula':>14} {'rel':>10}")
    for T in [0.10, 0.08, 0.06]:
        _, _, _, cap = solve_committor(
            F, L=1.8, n=361, T=T,
            wells_A=[(1.0, 0.0)], wells_B=[(-1.0, 0.0)], well_radius=0.22,
        )
        formula = T * math.exp(-1.0 / (4.0 * T))
        print(f"{T:>6.3f} {cap:>14.6e} {formula:>14.6e} {abs(cap-formula)/formula:>10.3f}")


def run_monkey():
    print("\n=== MONKEY SADDLE: F = (1/6) r^6 - (x^3 - 3 x y^2), three wells ===")
    F = lambda x, y: (x ** 2 + y ** 2) ** 3 / 6.0 - (x ** 3 - 3.0 * x * y ** 2)
    r0 = 3.0 ** (1.0 / 3.0)
    A = [(r0 * math.cos(0.0), r0 * math.sin(0.0))]
    B = [(r0 * math.cos(2 * math.pi / 3), r0 * math.sin(2 * math.pi / 3))]
    Cwell = (r0 * math.cos(4 * math.pi / 3), r0 * math.sin(4 * math.pi / 3))

    temps = [0.30, 0.25, 0.20, 0.15]
    caps, hC = [], []
    fields = {}
    for T in temps:
        X, Y, H, cap = solve_committor(
            F, L=2.0, n=321, T=T, wells_A=A, wells_B=B, well_radius=0.30,
        )
        caps.append(cap)
        # committor at the third well C
        d2 = (X - Cwell[0]) ** 2 + (Y - Cwell[1]) ** 2
        hC.append(float(H.flat[np.argmin(d2)]))
        fields[T] = (X, Y, H)

    caps = np.array(caps)
    temps = np.array(temps)
    # effective exponent p in cap ~ C T^p (no exponential, since H = 0 here)
    p = np.polyfit(np.log(temps), np.log(caps), 1)[0]

    print(f"{'T':>6} {'cap(A,B)':>14} {'h_C (should be 0.5)':>20}")
    for T, c, h in zip(temps, caps, hC):
        print(f"{T:>6.3f} {c:>14.6e} {h:>20.4f}")
    print(f"effective exponent p in cap ~ C T^p : p = {p:.3f}")
    print("(Morse/split formula is inapplicable: Hessian at the gate is zero.)")

    # figure: committor field at the smallest T
    Tfig = temps[-1]
    X, Y, H = fields[Tfig]
    fig, ax = plt.subplots(figsize=(6.0, 5.6))
    cf = ax.contourf(X, Y, H, levels=np.linspace(0, 1, 21), cmap="RdBu_r")
    ax.contour(X, Y, H, levels=[0.5], colors="k", linewidths=1.5)
    for (cx, cy), lab, col in [
        (A[0], "A (h=1)", "#b2182b"),
        (B[0], "B (h=0)", "#2166ac"),
        (Cwell, f"C (h≈{hC[-1]:.2f})", "#111111"),
    ]:
        ax.plot(cx, cy, "o", color=col, ms=9, mec="white", mew=1.5)
        ax.annotate(lab, (cx, cy), textcoords="offset points", xytext=(8, 8),
                    fontsize=10, color=col, fontweight="bold")
    ax.plot(0, 0, "x", color="black", ms=10, mew=2)
    ax.annotate("monkey saddle", (0, 0), textcoords="offset points",
                xytext=(6, -16), fontsize=9)
    ax.set_aspect("equal")
    ax.set_title(f"Branching-gate committor  h(x,y)   (T = {Tfig})", fontweight="bold")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.colorbar(cf, ax=ax, label="committor h")
    fig.tight_layout()
    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_PATH, dpi=150)
    print(f"Saved {FIGURE_PATH}")


def main() -> None:
    run_validation()
    run_monkey()


if __name__ == "__main__":
    main()
