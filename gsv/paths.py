"""Solution-path construction: per-formulation meshes + solves.

``build_path`` turns Phase-1 data into a discretized solution path
``{x*(s_j)}`` using the paper's mesh for each formulation and the existing
(paper-faithful, Gurobi) solvers in :mod:`gsv.formulations`. Returns
``(xx, s_values)`` with ``xx`` of shape (d, p). ``benchmark_solution`` returns the
formulation-specific literature calibration.

Meshes follow the manuscript (author-confirmed, moderate p):
  RO   : s_hat = (1-alpha)-quantile of Mahalanobis dist; s_j = (s_hat + off) * j/p
  SO   : s in {1..n1} (integer scenario counts; optionally subsampled to p)
  SAA  : s_j = step*j; num_constr = ceil(n1*(1-alpha+s_j))
  W-DRO: s in {3/2^j : j=0..9} u {(3/512)/1.2^j : j=1..15}
  FAST : x*(s) = (1-s) x_o + s x_fast, s = (j-1)/10, j=1..11
"""

from __future__ import annotations

import numpy as np
import gurobipy as gp

from . import formulations as F
from . import solvers_oss as OSS

__all__ = ["build_path", "benchmark_solution", "WASSERSTEIN_GRID"]

WASSERSTEIN_GRID = np.array([3 / 2 ** j for j in range(10)] + [(3 / 512) / 1.2 ** j for j in range(1, 16)])

# Gurobi (paper-faithful) primary; open-source fallback when the size-limited
# license rejects the model (d>=100, n>=2000). Convex optima are solver-identical.
_GUROBI = {"ro": F.CCP_RO_ellipsoid, "so": F.CCP_SO, "saa": F.CCP_SAA, "wass": F.CCP_DRO_wasserstein}
_OSS = {"ro": OSS.ro_ellipsoid, "so": OSS.so, "saa": OSS.saa, "wass": OSS.dro_wasserstein}


def _solve(name, *args):
    try:
        return _GUROBI[name](*args)
    except gp.GurobiError:
        return _OSS[name](*args)


def _mahalanobis_quantile(data, alpha):
    mu = data.mean(axis=0)
    Sigma = np.cov(data, rowvar=False)
    d = data.shape[1]
    inv = np.linalg.pinv(Sigma + 1e-10 * np.eye(d))
    diff = data - mu
    md = np.einsum("ij,jk,ik->i", diff, inv, diff)
    return float(np.quantile(md, 1 - alpha))


def build_path(formulation, data, c, b, alpha, d, n, mesh, rng=None):
    """Return (xx (d, p), s_values (p,)) for one Phase-1 dataset ``data`` (n, d)."""
    if formulation == "ro_ellipsoid":
        p = mesh.p or 25
        offset = float(mesh.params.get("offset", 20))
        s_hat = _mahalanobis_quantile(data, alpha)
        s_values = (s_hat + offset) * np.arange(1, p + 1) / p
        xx = _solve("ro", s_values, c, b, data, alpha, d, n)
        return xx, s_values

    if formulation == "so":
        s_values = np.arange(1, n + 1)
        if mesh.p and mesh.p < len(s_values):        # optional subsample to cap compute
            idx = np.unique(np.linspace(0, len(s_values) - 1, mesh.p).astype(int))
            s_values = s_values[idx]
        xx = _solve("so", s_values, c, b, data, alpha, d, n)
        return xx, s_values.astype(float)

    if formulation == "saa":
        # SAA tolerance-offset mesh s in (0, alpha]: the constraint is
        # (1/n1) sum 1(xi'x <= b) >= 1 - alpha + s, so s = alpha means "satisfy all
        # constraints" (num_constr = n1) at the conservative end. The paper used a
        # fixed step 0.002 that reaches s=alpha only when alpha=0.05; here we scale
        # the step to alpha/p so the conservative end always reaches all-constraints.
        p = mesh.p or 25
        step = float(mesh.params.get("step") or (alpha / p))
        s_values = step * np.arange(1, p + 1)
        num_constr = np.ceil(n * (1 - alpha + s_values)).astype(int)
        num_constr = np.clip(num_constr, 1, n)
        xx = _solve("saa", num_constr, c, b, data, alpha, d, n)
        return xx, s_values

    if formulation == "dro_wasserstein":
        s_values = WASSERSTEIN_GRID
        xx = _solve("wass", s_values, c, b, data, alpha, d, n)
        return xx, s_values

    if formulation == "fast":
        p = mesh.p or 11
        x_fast = _solve("so", np.array([n]), c, b, data, alpha, d, n)   # impose all n1 training constraints
        x_o = np.zeros(d)
        s_values = np.arange(p) / (p - 1)
        xx = np.stack([(1 - s) * x_o + s * x_fast for s in s_values], axis=1)
        return xx, s_values

    raise ValueError(f"no path builder for formulation {formulation!r}")


def benchmark_solution(benchmark, *, c, b, alpha, d, mu_true, sigma_true, full_data, n1):
    """Formulation-specific literature calibration (the 'benchmark' column)."""
    if benchmark is None:
        return None
    if benchmark == "SCA":
        try:
            from sca import SCA                            # true-moment SCA (Gurobi)
            return SCA(c, b, alpha, mu_true, sigma_true, lb=0.0, ub=1.0, output_flag=0)
        except gp.GurobiError:
            return OSS.sca(c, b, alpha, mu_true, sigma_true)
    if benchmark == "SO_all":
        n = full_data.shape[0]
        return _solve("so", np.array([n]), c, b, full_data, alpha, d, n)
    if benchmark == "FAST":
        from fast import SG_fast
        return SG_fast(c, b, full_data, N1=n1, output_flag=0)
    if benchmark == "chi2":
        raise NotImplementedError("moment-DRO chi2 benchmark needs cvxpy (solver env)")
    raise ValueError(f"unknown benchmark {benchmark!r}")
