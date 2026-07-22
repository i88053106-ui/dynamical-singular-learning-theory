from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


INPUT_PATH = Path("data/transverse_integral_julia.csv")
FIGURE_PATH = Path("figures/transverse_integral_ratios.png")


def main() -> None:
    if not INPUT_PATH.exists():
        raise FileNotFoundError(
            f"{INPUT_PATH} does not exist. Run the Julia experiment first."
        )

    data = pd.read_csv(INPUT_PATH).sort_values("T", ascending=False)

    print(data.to_string(index=False))

    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(8, 5))
    plt.semilogx(
        data["T"],
        data["raw_ratio"],
        marker="o",
        label="Raw ratio",
    )
    plt.semilogx(
        data["T"],
        data["corrected_ratio"],
        marker="s",
        label="Corrected ratio",
    )
    plt.axhline(
        data["target"].iloc[0],
        linestyle="--",
        label="sqrt(pi) / 3",
    )

    plt.gca().invert_xaxis()
    plt.xlabel("T")
    plt.ylabel("Scaled integral")
    plt.title("Transverse singular integral scaling")
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIGURE_PATH, dpi=180)

    print()
    print(f"Saved {FIGURE_PATH}")


if __name__ == "__main__":
    main()
