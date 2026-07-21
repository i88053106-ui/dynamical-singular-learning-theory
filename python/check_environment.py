from __future__ import annotations

import math
import sys

import matplotlib
import mpmath
import numpy as np
import pandas as pd
import scipy


def main() -> None:
    print("Python:", sys.version.split()[0])
    print("NumPy:", np.__version__)
    print("SciPy:", scipy.__version__)
    print("pandas:", pd.__version__)
    print("Matplotlib:", matplotlib.__version__)
    print("mpmath:", mpmath.__version__)
    print("sqrt(pi) / 3 =", math.sqrt(math.pi) / 3.0)
    print("Python environment OK")


if __name__ == "__main__":
    main()
