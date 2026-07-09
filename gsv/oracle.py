"""Oracle feasibility and path-oracle benchmarks.

The evaluation criterion needs a ground-truth notion of feasibility for a
candidate solution ``x`` under the *true* distribution:

* :func:`gaussian_feasibility` — exact closed form ``P(xi'x <= b) =
  Phi((b - mu'x)/sqrt(x'Sigma x))`` when ``xi ~ N(mu, Sigma)``.
* :func:`large_sample_feasibility` — an essentially-exact estimate from a large
  independent sample; valid for any DGP (half-normal, multivariate-t, ...).
* :func:`true_feasibility` — dispatches to the closed form for Gaussian, else the
  large sample.

On a fixed solution path ``{x*(s_j)}`` these give the **path oracle**
(:func:`path_oracle`): the lowest-objective candidate that is *truly* feasible,
and the least-conservative feasible ``s`` — used for objective-gap and
excess-conservativeness metrics.

:func:`exact_gaussian_ccp_solution` (the *global* oracle SOCP) needs a conic
solver (cvxpy, lazily imported) and therefore only runs in a solver environment.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import norm

from .config import DGP
from . import dgp as _dgp

__all__ = [
    "gaussian_feasibility", "large_sample_feasibility", "true_feasibility",
    "path_oracle", "exact_gaussian_ccp_solution",
]

_EPS = 1e-12


def _as_2d(x: np.ndarray) -> tuple[np.ndarray, bool]:
    x = np.asarray(x, dtype=float)
    if x.ndim == 1:
        return x[:, None], True
    return x, False


def gaussian_feasibility(x: np.ndarray, mu: np.ndarray, Sigma: np.ndarray, b: float) -> np.ndarray | float:
    """Exact ``P(xi'x <= b)`` for ``xi ~ N(mu, Sigma)``. ``x`` is (d,) or (d, p)."""
    X, scalar = _as_2d(x)
    mx = mu @ X                                   # (p,)
    quad = np.einsum("ij,ij->j", Sigma @ X, X)    # x'Sigma x per column
    denom = np.sqrt(np.maximum(quad, 0.0))
    out = np.empty(X.shape[1])
    nz = denom > _EPS
    out[nz] = norm.cdf((b - mx[nz]) / denom[nz])
    # degenerate variance: deterministic constraint
    out[~nz] = (mx[~nz] <= b).astype(float)
    return float(out[0]) if scalar else out


def large_sample_feasibility(x: np.ndarray, eval_sample: np.ndarray, b: float) -> np.ndarray | float:
    """Empirical ``P(xi'x <= b)`` over a large independent ``eval_sample`` (m, d)."""
    X, scalar = _as_2d(x)
    vals = eval_sample @ X                         # (m, p)
    out = np.mean(vals <= b, axis=0)
    return float(out[0]) if scalar else out


def true_feasibility(x: np.ndarray, dgp: DGP, d: int, b: float,
                     eval_sample: np.ndarray | None = None) -> np.ndarray | float:
    """Ground-truth feasibility: closed form for Gaussian, large sample otherwise."""
    if dgp.kind == "gaussian":
        mu, Sigma = _dgp.moments(dgp, d)
        return gaussian_feasibility(x, mu, Sigma, b)
    if eval_sample is None:
        raise ValueError(f"large-sample oracle needs eval_sample for DGP {dgp.kind!r}")
    return large_sample_feasibility(x, eval_sample, b)


def path_oracle(s_values: np.ndarray, objectives: np.ndarray, feasibilities: np.ndarray,
                target: float) -> dict:
    """Best achievable outcome on a fixed solution path under the true feasibility.

    Returns the lowest-objective *truly feasible* candidate (the path oracle
    solution), and the least-conservative feasible ``s`` for excess-conservativeness.
    """
    s_values = np.asarray(s_values, dtype=float)
    objectives = np.asarray(objectives, dtype=float)
    feasibilities = np.asarray(feasibilities, dtype=float)
    feasible = feasibilities >= target
    if not feasible.any():
        return {"any_feasible": False, "best_obj_idx": None, "best_obj": np.nan,
                "min_feasible_s": np.nan, "min_feasible_s_idx": None}
    idx_feasible = np.flatnonzero(feasible)
    best = idx_feasible[np.argmin(objectives[idx_feasible])]
    min_s = idx_feasible[np.argmin(s_values[idx_feasible])]
    return {"any_feasible": True, "best_obj_idx": int(best), "best_obj": float(objectives[best]),
            "min_feasible_s": float(s_values[min_s]), "min_feasible_s_idx": int(min_s)}


def exact_gaussian_ccp_solution(c, mu, Sigma, b, alpha, ub=1.0, lb=0.0):
    """Global oracle: exact CCP solution under Gaussian xi.

    Solves ``min c'x  s.t.  mu'x + z_{1-alpha} ||Sigma^{1/2} x||_2 <= b, lb <= x <= ub``
    (an SOCP). Requires cvxpy (solver environment).
    """
    import cvxpy as cp  # lazy: solver env only

    c = np.asarray(c, float); mu = np.asarray(mu, float); Sigma = np.asarray(Sigma, float)
    d = len(c)
    z = float(norm.ppf(1.0 - alpha))
    root = np.real(np.linalg.cholesky(Sigma + 1e-12 * np.eye(d)))
    x = cp.Variable(d)
    constraints = [mu @ x + z * cp.norm(root.T @ x, 2) <= b, x >= lb, x <= ub]
    prob = cp.Problem(cp.Minimize(c @ x), constraints)
    prob.solve()
    return np.asarray(x.value, float)
