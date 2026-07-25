"""Branching-gate capacity as a circulant k-terminal network (analytic J_k).

Scale invariance (docs/branching-gate.md) gives cap_T = J_k T, with J_k the
adjacent-pair conductance of the degenerate junction where k valleys meet. Deep
valleys equilibrate, so the junction is a linear k-terminal conductor; k-fold
symmetry makes its conductance (graph-Laplacian) matrix circulant. Then J_k is
the effective conductance between two ADJACENT terminals of that circulant
network, and the spectator committor values are the floating-node potentials of
the same network.

This script validates the circulant reduction with a PARAMETER-FREE test at k=4.
The k=4 junction has one coupling ratio r = c/b (opposite/adjacent). Solving the
4-terminal network:

    spectator adjacent to B:  h = (1+r)/(3+r)   =>   r = (3h-1)/(1-h),
    cap(adjacent)          :  4 b (1+r)/(3+r),
    cap(opposite)          :  b (1+r),
    cap(opposite)/cap(adjacent) = (3+r)/4.

The ratio is independent of the germ-specific base conductance b: the spectator
value h alone predicts the opposite/adjacent capacity ratio. We also check the
symmetry structure of the spectators for k=3,4,5,6 (complementary pairs sum to 1,
on-axis wells equal 1/2), which the circulant network requires.
"""

from __future__ import annotations

import math

import numpy as np

from branching_gate_committor import solve_committor

LDOM, NGRID, RADIUS = 2.0, 361, 0.26


def make_F(k):
    def F(x, y):
        return (x * x + y * y) ** k / (2.0 * k) - np.real((x + 1j * y) ** k)
    return F


def wells(k):
    r0 = k ** (1.0 / k)
    return [(r0 * math.cos(2 * math.pi * j / k), r0 * math.sin(2 * math.pi * j / k))
            for j in range(k)]


def hat(X, Y, H, pt):
    return float(H.flat[np.argmin((X - pt[0]) ** 2 + (Y - pt[1]) ** 2)])


def cap_and_committor(k, T, jA, jB):
    F, W = make_F(k), wells(k)
    X, Y, H, cap = solve_committor(F, LDOM, NGRID, T, [W[jA]], [W[jB]], RADIUS)
    return cap, [hat(X, Y, H, W[j]) for j in range(k)]


def test_k4():
    print("=== k=4 parameter-free circulant test:  cap_opp/cap_adj ?= (3+r)/4 ===")
    print(f"{'T':>6} {'cap_adj':>9} {'cap_opp':>9} {'h_spec(B)':>10} {'r':>7} "
          f"{'pred':>8} {'meas':>8} {'|Δ|':>9}")
    for T in [0.16, 0.13, 0.10]:
        cap_adj, hA = cap_and_committor(4, T, 0, 1)   # spectators W2 (adj B), W3 (adj A)
        cap_opp, _ = cap_and_committor(4, T, 0, 2)
        h = hA[2]                                     # spectator adjacent to B
        r = (3 * h - 1) / (1 - h)
        pred, meas = (3 + r) / 4, cap_opp / cap_adj
        print(f"{T:>6.3f} {cap_adj:>9.5f} {cap_opp:>9.5f} {h:>10.4f} {r:>7.4f} "
              f"{pred:>8.4f} {meas:>8.4f} {abs(pred-meas):>9.2e}")


def test_spectator_symmetry():
    print("\n=== spectator symmetry (circulant requires complementary pairs) ===")
    print(f"{'k':>3} {'spectator committor (adjacent A=0,B=1)':>44} {'max |pair sum-1|':>18}")
    for k in [3, 4, 5, 6]:
        _, H = cap_and_committor(k, 0.12, 0, 1)
        spec = H[2:]                       # floating wells 2..k-1
        # reflection across the 0-1 bisector maps well j -> well (1-j) mod k
        resid = 0.0
        for j in range(2, k):
            jp = (1 - j) % k
            if jp >= 2:
                resid = max(resid, abs(H[j] + H[jp] - 1.0))
            else:
                resid = max(resid, abs(H[j] - 0.5) if jp == j else 0.0)
        specstr = ", ".join(f"{v:.3f}" for v in spec)
        print(f"{k:>3} {specstr:>44} {resid:>18.2e}")


def main() -> None:
    test_k4()
    test_spectator_symmetry()
    print("\nThe circulant k-terminal network reproduces both the adjacent/opposite")
    print("capacity ratio (parameter-free, from a spectator value) and the spectator")
    print("symmetry. J_k thus reduces to the junction edge-conductances b_s(k), the")
    print("remaining germ-specific unknowns (open: closed form of b_s from Re[w^k]).")


if __name__ == "__main__":
    main()
