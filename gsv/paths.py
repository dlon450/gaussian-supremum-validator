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
try:
    import gurobipy as gp
    _HAS_GUROBI = True
    _GErr = gp.GurobiError
except Exception:                        # Gurobi optional: fall back to OSS backend
    _HAS_GUROBI = False
    class _GErr(Exception):
        pass

from . import formulations as F
from . import solvers_oss as OSS

# The open-source conic solver (CLARABEL) spawns one thread per core, ignoring the
# OMP/BLAS env caps -> catastrophic oversubscription when many workers each fall back
# to it (e.g. RO d>=100, which exceeds the size-limited Gurobi license). Pin it to a
# single thread. Guarded import so the module still loads without threadpoolctl.
try:
    from threadpoolctl import threadpool_limits
except Exception:  # pragma: no cover
    from contextlib import contextmanager
    @contextmanager
    def threadpool_limits(limits=None):
        yield

__all__ = ["build_path", "benchmark_solution", "WASSERSTEIN_GRID"]

WASSERSTEIN_GRID = np.array([3 / 2 ** j for j in range(10)] + [(3 / 512) / 1.2 ** j for j in range(1, 16)])

# Gurobi (paper-faithful) primary; open-source fallback when the size-limited
# license rejects the model (d>=100, n>=2000). Convex optima are solver-identical.
_GUROBI = {"ro": F.CCP_RO_ellipsoid, "so": F.CCP_SO, "saa": F.CCP_SAA, "wass": F.CCP_DRO_wasserstein}
_OSS = {"ro": OSS.ro_ellipsoid, "so": OSS.so, "saa": OSS.saa, "wass": OSS.dro_wasserstein}


def _solve(name, *args):
    if _HAS_GUROBI:
        try:
            return _GUROBI[name](*args)
        except _GErr:                          # size-limited license -> OSS backend
            pass
    with threadpool_limits(limits=1):          # cap CLARABEL/HiGHS oversubscription
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

    if formulation == "dro_moment":
        # Moment-based DRO (delta-method confidence region). The conservativeness
        # parameter s is the uncertainty-set radius rho; the path spans 0..1.5*hat_s
        # where hat_s = sqrt(chi2_{0.95, q}) is the benchmark's 95% radius over the
        # q = d + d(d+1)/2 moment parameters (paper Section 6.2). SDP via cvxpy;
        # thread-capped to avoid the conic-solver oversubscription.
        p = mesh.p or 25
        scale = float(mesh.params.get("scale", 1.5))
        # radius uses the moment-estimation sample size n (= n1 here). Include s=0
        # (moments as equalities) at the aggressive end of the path.
        rho95 = _moment_chi2_radius(d, n)
        s_values = np.concatenate([[0.0], scale * rho95 * np.arange(1, p + 1) / p])
        with threadpool_limits(limits=1):
            xx = F.CCP_DRO_moment(s_values, c, b, data, alpha, d, n)
        return xx, s_values

    raise ValueError(f"no path builder for formulation {formulation!r}")


def _moment_chi2_radius(d: int, n: int) -> float:
    """95% radius hat_s = sqrt(chi2_{0.95, q} / n) of the delta-method moment
    confidence region over q = d + d(d+1)/2 parameters (mean + lower-triangular
    second moments). The 1/sqrt(n) factor scales the single-observation moment
    covariance V_est to the covariance of the sample-mean moment estimate; it
    matches the legacy implementation (DRO2.m: rho = sqrt(chi2inv(1-beta,q)/N))."""
    from scipy.stats import chi2
    q = d + d * (d + 1) // 2
    return float(np.sqrt(chi2.ppf(0.95, q) / n))


def benchmark_solution(benchmark, *, c, b, alpha, d, mu_true, sigma_true, full_data, n1):
    """Formulation-specific literature calibration (the 'benchmark' column)."""
    if benchmark is None:
        return None
    if benchmark == "SCA":
        if _HAS_GUROBI:
            try:
                from sca import SCA                        # true-moment SCA (Gurobi)
                return SCA(c, b, alpha, mu_true, sigma_true, lb=0.0, ub=1.0, output_flag=0)
            except _GErr:
                pass
        with threadpool_limits(limits=1):                 # cap CLARABEL oversubscription
            return OSS.sca(c, b, alpha, mu_true, sigma_true)
    if benchmark == "SO_all":
        n = full_data.shape[0]
        return _solve("so", np.array([n]), c, b, full_data, alpha, d, n)
    if benchmark == "FAST":
        from fast import SG_fast
        return SG_fast(c, b, full_data, N1=n1, output_flag=0)
    if benchmark == "chi2":
        # literature calibration: rho = 95% radius of the moment confidence region,
        # estimated from all n observations (delta-method, sqrt(chi2/n)).
        n = full_data.shape[0]
        with threadpool_limits(limits=1):
            x = F.CCP_DRO_moment(np.array([_moment_chi2_radius(d, n)]), c, b, full_data, alpha, d, n)
        return np.asarray(x, float)
    raise ValueError(f"unknown benchmark {benchmark!r}")
