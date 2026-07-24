"""Langevin mean first-passage time for the first-target potential.

Overdamped Langevin dynamics

    dW = -grad F(W) dt + sqrt(2 T) dB,

    F(x,y,z) = (x^2 - 1)^2/4 + x^2 (y^2 + z^2) + y^6 + y^2 z^2 + z^6,

with wells at (+-1, 0, 0) and a Newton-degenerate gate at the origin
(communication height H = 1/4). Trajectories start at the left well
(-1, 0, 0); tau is the first-passage time to the separatrix {x >= 0}.

The split-gate capacity derivation (docs/capacity-derivation.md) predicts
the well-to-well Eyring-Kramers time

    E[tau_EK] ~ 3 sqrt(2) pi^{3/2} sqrt(T) / log(1/T) * exp(1/(4T)).

By the classical 1D relation, first passage to the separatrix is half of
that:

    E[tau_{-> 0}] ~ (1/2) * 3 sqrt(2) pi^{3/2} sqrt(T)/log(1/T) * exp(1/(4T)).

This is an independent Monte-Carlo check of the derived clock law. At the
temperatures reachable by direct simulation the barrier is shallow, so
o(1) corrections are large; the test is the *trend* (does the raw ratio
approach the predicted constant, and is the (log(1/T))^{-1} signature
visible), not a machine-precision match.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

OUTPUT_PATH = Path("data/langevin_exit_time.csv")

# Predicted well-to-well prefactor constant: 3 * sqrt(2) * pi^{3/2}.
EK_CONSTANT = 3.0 * math.sqrt(2.0) * math.pi ** 1.5


def grad_F(W: np.ndarray) -> np.ndarray:
    """Gradient of F, W has shape (n, 3)."""
    x = W[:, 0]
    y = W[:, 1]
    z = W[:, 2]

    gx = x ** 3 - x + 2.0 * x * (y ** 2 + z ** 2)
    gy = 2.0 * x ** 2 * y + 6.0 * y ** 5 + 2.0 * y * z ** 2
    gz = 2.0 * x ** 2 * z + 6.0 * z ** 5 + 2.0 * z * y ** 2

    return np.stack((gx, gy, gz), axis=1)


def mean_first_passage(
    T: float,
    n_traj: int,
    dt: float,
    max_steps: int,
    rng: np.random.Generator,
) -> tuple[float, float, int, int]:
    """Return (mean tau, standard error, n_absorbed, n_traj)."""
    W = np.zeros((n_traj, 3))
    W[:, 0] = -1.0  # start at the left well

    tcross = np.full(n_traj, -1, dtype=np.int64)
    active = np.ones(n_traj, dtype=bool)

    noise_scale = math.sqrt(2.0 * T * dt)

    for step in range(1, max_steps + 1):
        idx = np.nonzero(active)[0]
        if idx.size == 0:
            break

        Wa = W[idx]
        Wa = Wa - grad_F(Wa) * dt + noise_scale * rng.standard_normal(Wa.shape)
        W[idx] = Wa

        crossed_local = Wa[:, 0] >= 0.0
        if crossed_local.any():
            crossed_idx = idx[crossed_local]
            tcross[crossed_idx] = step
            active[crossed_idx] = False

    absorbed = tcross >= 0
    n_absorbed = int(absorbed.sum())

    times = tcross[absorbed] * dt
    if n_absorbed == 0:
        return math.nan, math.nan, 0, n_traj

    mean_tau = float(times.mean())
    std_err = float(times.std(ddof=1) / math.sqrt(n_absorbed)) if n_absorbed > 1 else math.nan

    return mean_tau, std_err, n_absorbed, n_traj


def predicted_separatrix_time(T: float) -> float:
    """(1/2) * EK well-to-well time, leading order."""
    return 0.5 * EK_CONSTANT * math.sqrt(T) / math.log(1.0 / T) * math.exp(1.0 / (4.0 * T))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--temps", type=float, nargs="+",
                        default=[0.06, 0.055, 0.05, 0.045])
    parser.add_argument("--n-traj", type=int, default=500)
    parser.add_argument("--dt", type=float, default=1e-3)
    parser.add_argument("--max-steps", type=int, default=4_000_000)
    parser.add_argument("--seed", type=int, default=20260723)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    rows = []
    for T in args.temps:
        mean_tau, std_err, n_abs, n_tot = mean_first_passage(
            T, args.n_traj, args.dt, args.max_steps, rng
        )
        predicted = predicted_separatrix_time(T)
        # g(T) = log E[tau] - 1/(4T) - (1/2) log T ; predicted -> -log log(1/T) + const.
        g = math.log(mean_tau) - 1.0 / (4.0 * T) - 0.5 * math.log(T)
        g_plus_loglog = g + math.log(math.log(1.0 / T))
        ratio = mean_tau / predicted

        rows.append(dict(
            T=T, mean_tau=mean_tau, std_err=std_err,
            n_absorbed=n_abs, n_traj=n_tot,
            predicted=predicted, ratio=ratio, g_plus_loglog=g_plus_loglog,
        ))

        print(
            f"T={T:6.4f}  E[tau]={mean_tau:12.3f} +/- {std_err:8.3f}  "
            f"pred={predicted:12.3f}  ratio={ratio:7.4f}  "
            f"g+loglog={g_plus_loglog:8.4f}  ({n_abs}/{n_tot} absorbed)"
        )

    df = pd.DataFrame(rows)
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_PATH, index=False)

    print()
    print("ratio -> 1 tests the derived prefactor constant.")
    print("g+loglog -> const (flat in T) tests the (log(1/T))^{-1} signature.")
    print(f"Saved {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
