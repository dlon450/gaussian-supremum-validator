"""CVaR-constrained portfolio: the paper's *second* problem class (a general
stochastic constraint rather than a chance constraint).

Problem (paper eq. (cvar)):   min  c'x   s.t.  CVaR_alpha(-xi'x) <= gamma,
                                    sum(x)=1, x>=0,   with c = -E[xi].

The data-driven DRO reformulation over a type-1 Wasserstein ball of radius s
(paper eq. (cvar:reformulation), Esfahani-Kuhn Cor. 5.1) is *convex* (no big-M):

    min c'x
    s.t. s*v + (1/n) sum_i z_i <= alpha * r
         r <= z_i + gamma + xi_i' x      for i=1..n
         v >= x_j                        for j=1..d     (||x||_inf <= v, x>=0)
         sum(x)=1, x>=0, v,r,z >= 0

The conservativeness parameter is the radius s (s=0 => sample CVaR / aggressive;
larger s => more robust => lower true CVaR but lower return). The two-phase
validators pick the least-conservative s whose CVaR is certified <= gamma on
Phase-2 data, exactly as in the chance-constraint case but with the CVaR estimator
and its influence-function standard error in place of the feasibility indicator.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import norm

try:
    from threadpoolctl import threadpool_limits
except Exception:  # pragma: no cover
    from contextlib import contextmanager
    @contextmanager
    def threadpool_limits(limits=None):
        yield

__all__ = ["dro_cvar_path", "true_cvar_gaussian", "empirical_cvar", "cvar_phase2_stats",
           "select_cvar", "WASSERSTEIN_RADII"]

# Radius grid on the scale of the returns (~0.1): s=0 is the aggressive sample-CVaR
# end; larger s is more robust (lower true CVaR, lower return).
WASSERSTEIN_RADII = np.concatenate([[0.0], np.geomspace(1e-4, 0.15, 24)])


def _cvar_reformulation(s, c, gamma, data, alpha, d, n):
    """Solve the Wasserstein-DRO CVaR-constrained portfolio at radius s (convex)."""
    import cvxpy as cp
    x = cp.Variable(d); v = cp.Variable(nonneg=True); r = cp.Variable(nonneg=True)
    z = cp.Variable(n, nonneg=True)
    cons = [s * v + cp.sum(z) / n <= alpha * r,
            r <= z + gamma + data @ x,           # r <= z_i + gamma + xi_i'x  (vectorized)
            v >= x, cp.sum(x) == 1, x >= 0]
    prob = cp.Problem(cp.Minimize(c @ x), cons)
    try:
        prob.solve(solver=cp.ECOS)
        if prob.status not in ("optimal", "optimal_inaccurate") or x.value is None:
            raise cp.error.SolverError
    except Exception:
        prob.solve(solver=cp.CLARABEL)
    return x.value if x.value is not None else np.full(d, 1.0 / d)


def dro_cvar_path(s_values, c, gamma, data, alpha, d, n):
    """Solution path x*(s) over the Wasserstein radius grid; returns (d, p)."""
    s_values = np.atleast_1d(s_values)
    with threadpool_limits(limits=1):
        cols = [_cvar_reformulation(float(s), c, gamma, data, alpha, d, n) for s in s_values]
    return np.array(cols, float).T


def true_cvar_gaussian(x, mu, Sigma, alpha):
    """Closed-form CVaR_alpha of the loss L = -xi'x for xi ~ N(mu, Sigma).

    L ~ N(-mu'x, x'Sigma x); CVaR_alpha(L) = mean_L + std_L * phi(z_{1-alpha})/alpha.
    ``x`` is (d,) or (d,p) -> scalar or (p,)."""
    X = np.atleast_2d(x.T).T if x.ndim == 1 else x
    mL = -(mu @ X)                                   # (p,)
    sL = np.sqrt(np.maximum(np.einsum("ij,ij->j", Sigma @ X, X), 0.0))
    zq = norm.ppf(1 - alpha)
    cvar = mL + sL * norm.pdf(zq) / alpha
    return float(cvar[0]) if x.ndim == 1 else cvar


def empirical_cvar(losses, alpha):
    """Rockafellar-Uryasev empirical CVaR and its influence-function values.

    CVaR_hat = VaR_hat + mean((L - VaR_hat)_+)/alpha,  VaR_hat = (1-alpha)-quantile.
    Influence fn IF_i = VaR_hat + (L_i - VaR_hat)_+/alpha - CVaR_hat (mean 0), whose
    sample std / sqrt(n) is the CVaR standard error."""
    L = np.asarray(losses, float)
    var = np.quantile(L, 1 - alpha)
    exc = np.maximum(L - var, 0.0)
    cvar = var + exc.mean() / alpha
    IF = var + exc / alpha - cvar
    return cvar, IF


def cvar_phase2_stats(path, ksi_val, alpha):
    """Per-candidate CVaR estimate, standard error, and IF matrix on Phase-2 data.

    ``path`` (d,p), ``ksi_val`` (n2,d). Loss for candidate j is -ksi_val @ x_j.
    Returns (cvar_hat (p,), se (p,), IFmat (n2,p))."""
    losses = -(ksi_val @ path)                        # (n2, p)
    n2, p = losses.shape
    cvar = np.empty(p); IFmat = np.empty((n2, p))
    for j in range(p):
        c_j, IF_j = empirical_cvar(losses[:, j], alpha)
        cvar[j] = c_j; IFmat[:, j] = IF_j
    se = IFmat.std(axis=0, ddof=1) / np.sqrt(n2)
    return cvar, se, IFmat


def select_cvar(cvar_hat, se, IFmat, objectives, *, gamma, beta, n2, rng, sim_num=2000):
    """Pick the lowest-objective candidate whose CVaR is certified <= gamma.

    NV : cvar_hat <= gamma                  (no margin)
    UG : cvar_hat + z_{1-beta}*se <= gamma  (univariate Gaussian)
    NGS: cvar_hat + qhat*se <= gamma        (normalized Gaussian supremum)
    UNGS: cvar_hat + qhat_u/sqrt(n2) <= gamma (unnormalized)
    Feasibility here means "safely below the risk threshold", so the margin is
    ADDED to the CVaR estimate (mirrors subtracting it from a feasibility level)."""
    obj = np.asarray(objectives, float)
    Sigma = np.cov(IFmat.T)                            # (p,p) IF covariance
    Sigma = np.atleast_2d(Sigma)
    sd = np.sqrt(np.clip(np.diag(Sigma), 0, None))     # = se * sqrt(n2)
    nz = sd > 0
    z = float(norm.ppf(1 - beta))
    if np.any(nz):
        S = Sigma[np.ix_(nz, nz)]
        Z = rng.multivariate_normal(np.zeros(S.shape[0]), S, sim_num)
        tmat = np.tile(sd[nz], (sim_num, 1))
        qhat = float(np.percentile(np.amax(Z / tmat, axis=1), (1 - beta) * 100))
        qhat_u = float(np.percentile(np.amax(Z, axis=1), (1 - beta) * 100))
    else:
        qhat = qhat_u = 0.0

    def pick(margin):
        feas = cvar_hat + margin <= gamma
        if np.any(feas):
            idx = np.flatnonzero(feas); return int(idx[np.argmin(obj[idx])]), True
        return int(np.argmin(cvar_hat + margin)), False   # fallback: "safest"

    out = {}
    out["NV"] = pick(np.zeros_like(cvar_hat))
    out["UG"] = pick(z * se)
    out["NGS"] = pick(qhat * se)
    out["UNGS"] = pick(qhat_u / np.sqrt(n2) * np.ones_like(cvar_hat))
    return out
