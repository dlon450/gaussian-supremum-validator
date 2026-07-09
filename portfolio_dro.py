"""Demo: mean-risk (CVaR) portfolio under a Wasserstein ambiguity set.

Generates synthetic yearly-return data and solves the mean-risk portfolio DRO
formulation (see ``formulations.mean_risk_portfolio_dro``) over a grid of
ambiguity radii, printing the optimal weights for each radius.

Note
----
The original version of this file also carried its own broken copies of the
validation routines (``algo1/2/3``, ``GS``, ``perform_stats``): they referenced
undefined names (``pd``, ``np.matlib``, ``N``, ``K``, ``B``, ``c``, ``b``), used
an unindexable ``set`` for the parameter grid, and wrote to a mislabelled CSV, so
that experiment never ran. Wiring the mean-risk portfolio into the chance-
constraint validators (:mod:`gsv.validators`) would require defining cost/threshold
semantics that were never specified, so this file now just exercises the solver
directly. The solver itself has been moved to :mod:`gsv.formulations` for reuse.

Run with:  python portfolio_dro.py
"""

import numpy as np

from gsv.formulations import mean_risk_portfolio_dro

n, d = 100, 10
alpha = 0.2
beta = 0.1

# Synthetic historical data of yearly returns.
i = np.arange(1, d + 1)
np.random.seed(1)
phi = 0.02 * np.random.normal(size=(n, d))
zeta = 0.03 * i + 0.025 * i * np.random.normal(size=(n, d))
zhat = np.maximum(phi + zeta, -1)

# Grid of Wasserstein ambiguity radii.
delta = np.array(sorted({base * 10 ** exp for base in range(10) for exp in [-3, -2, -1]}))

solution = mean_risk_portfolio_dro(delta, zhat, n=n, d=d, alpha=alpha)
for k, radius in enumerate(delta):
    print(f'radius = {radius:g}: weights = {np.round(solution[:, k], 4)}')
