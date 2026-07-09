"""Gaussian-supremum validation for data-driven chance-constrained programs.

Shared library for the experiments accompanying the paper (arXiv:1909.06477):

- :mod:`gsv.formulations` — base optimization formulations (the solvers).
- :mod:`gsv.validators`   — data-driven validation / parameter-selection methods.

Submodules are imported explicitly (e.g. ``from gsv.formulations import CCP_SO``)
so that importing :mod:`gsv` does not pull in the optional heavy solver
dependencies (gurobipy, cvxpy, rsome).
"""
