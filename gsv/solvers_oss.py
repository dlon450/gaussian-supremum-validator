"""Open-source (license-free) solver backend, mirroring gsv.formulations.

Used to run sizes beyond the free Gurobi license cap (d>=100, n>=2000). For the
convex formulations (RO SOCP, SO/FAST LP) the optimal objective is
solver-independent, so CLARABEL/HiGHS reproduce Gurobi's results to solver
tolerance. For the MILPs (SAA, Wasserstein) HiGHS also returns the global optimum
when it converges, but may be slower and can pick a different optimal vertex.

Each function matches the corresponding CCP_* signature
``solve(para, c, b, data, alpha, d, n) -> (d,) or (d, p)`` so it is a drop-in
backend for gsv.paths.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import sqrtm

__all__ = ["ro_ellipsoid", "so", "saa", "dro_wasserstein", "sca"]


def _solve_socp(prob):
    """Solve a convex SOCP, preferring the single-threaded ECOS to avoid the
    massive thread oversubscription CLARABEL causes under many parallel workers
    (RO d>=100 falls back here). Falls back to CLARABEL if ECOS is unavailable/fails."""
    import cvxpy as cp
    try:
        prob.solve(solver=cp.ECOS)
        if prob.status in ("optimal", "optimal_inaccurate"):
            return
    except Exception:
        pass
    prob.solve(solver=cp.CLARABEL)


def sca(c, b, alpha, mu, sigma, lb=0.0, ub=1.0, output_flag=0):
    """SCA benchmark (true-moment SOCP), open-source drop-in for sca.SCA."""
    import cvxpy as cp
    c = np.asarray(c, float); mu = np.asarray(mu, float); sigma = np.asarray(sigma, float)
    d = len(c)
    root = np.real(sqrtm(sigma))
    phi = np.sqrt(2.0 * np.log(1.0 / alpha))
    x = cp.Variable(d)
    prob = cp.Problem(cp.Minimize(c @ x),
                      [mu @ x + phi * cp.norm(root @ x, 2) <= b, x >= lb, x <= ub])
    _solve_socp(prob)
    return x.value


def _finish(sols, scalar):
    out = np.array(sols, dtype=float).T
    # Match the Gurobi formulations' convention: a single-parameter mesh returns a
    # 1-D (d,) vector, not (d, 1). (The Gurobi solvers use `... else solution[:, 0]`.)
    # Without this, a 1-element `para` (e.g. the SO_all benchmark passing [n]) that
    # falls back to OSS at n>=2000 returns (d, 1) and breaks downstream `float(c @ x)`.
    return out[:, 0] if (scalar or out.shape[1] == 1) else out


def ro_ellipsoid(para, c, b, data, alpha, d, n):
    import cvxpy as cp
    mu = np.mean(data, axis=0)
    root = np.real(sqrtm(np.cov(data, rowvar=False)))
    scalar = isinstance(para, (float, np.floating))
    para = np.atleast_1d(para)
    sols = []
    for s in para:
        x = cp.Variable(d)
        prob = cp.Problem(cp.Minimize(c @ x),
                          [mu @ x + np.sqrt(s) * cp.norm(root @ x, 2) <= b, x >= 0, x <= 1])
        _solve_socp(prob)
        sols.append(x.value)
    return _finish(sols, scalar)


def so(para, c, b, data, alpha, d, n):
    import cvxpy as cp
    scalar = np.ndim(para) == 0
    para = np.atleast_1d(para).astype(int)
    sols = []
    for s in para:
        x = cp.Variable(d)
        cons = [x >= 0, x <= 1]
        if s > 0:
            cons.append(data[:s] @ x <= b)
        prob = cp.Problem(cp.Minimize(c @ x), cons)
        prob.solve(solver=cp.HIGHS)
        sols.append(x.value)
    return _finish(sols, scalar)


def saa(para, c, b, data, alpha, d, n):
    import cvxpy as cp
    M = float(np.ceil(np.abs(b) + np.max(np.sum(np.abs(data), axis=1))))
    scalar = np.ndim(para) == 0
    num_constr = np.atleast_1d(para).astype(int)
    sols = []
    prev = None
    for k, nc in enumerate(num_constr):
        if k > 0 and nc == num_constr[k - 1]:
            sols.append(prev); continue
        x = cp.Variable(d); z = cp.Variable(n, boolean=True)
        cons = [x >= 0, x <= 1, data @ x <= b + M * (1 - z), cp.sum(z) >= nc]
        prob = cp.Problem(cp.Minimize(c @ x), cons)
        prob.solve(solver=cp.HIGHS)
        prev = x.value if x.value is not None else np.full(d, np.nan)
        sols.append(prev)
    return _finish(sols, scalar)


def dro_wasserstein(para, c, b, data, alpha, d, n):
    import cvxpy as cp
    ksi = data
    M = float(abs(b) + max(np.sum(ksi, 1)))
    scalar = isinstance(para, (float, np.floating))
    para = np.atleast_1d(para)
    sols = []
    for s in para:
        x = cp.Variable(d); y = cp.Variable(n, boolean=True)
        z = cp.Variable(n, nonneg=True); sl = cp.Variable(n, nonneg=True)
        v = cp.Variable(nonneg=True); r = cp.Variable(nonneg=True)
        cons = [x >= 0, x <= 1,
                z + sl >= r,
                s * v + cp.sum(z) / n <= alpha * r,
                b - ksi @ x + M * (1 - y) >= sl,
                sl <= M * y,
                v >= x]                       # ||x||_* (l_inf here) <= v via v >= x_i (x>=0)
        prob = cp.Problem(cp.Minimize(c @ x), cons)
        prob.solve(solver=cp.HIGHS)
        sols.append(x.value if x.value is not None else np.full(d, np.nan))
    return _finish(sols, scalar)
