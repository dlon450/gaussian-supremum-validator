"""Pure validator selection rules on a shared solution path.

Separates *selection* from *path construction* (unlike the monolithic
``gsv.validators.gaussian_supremum``): given the Phase-2 sample statistics of a
fixed candidate path, each rule returns the chosen candidate index. This lets the
orchestrator build the path once and evaluate every validator plus the oracle on
the same candidates, and it makes the RNG explicit (injected ``Generator``)
instead of NumPy's global state.

Math mirrors the validated ``gaussian_supremum`` exactly:
  * UG   (Algorithm 4, focal): margin z_{1-beta} * sigma_j / sqrt(n2)
  * NGS  (Algorithm 3): margin qhat * sigma_j / sqrt(n2)
  * UNGS (Algorithm 2): margin qhat_u / sqrt(n2)
  * NV   (naive): no margin
Each picks the lowest-objective candidate clearing its feasibility margin; if none
clear it, it falls back to the "most feasible" candidate (argmax of the margin-
adjusted feasibility), matching the authors' documented fallback.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import norm

__all__ = ["phase2_stats", "gs_quantiles", "select_all"]


def phase2_stats(path: np.ndarray, ksi_val: np.ndarray, b: float):
    """Per-candidate feasibility mean, std, and full covariance on Phase-2 data.

    ``path`` is (d, p); ``ksi_val`` is (n2, d). Returns (P_hat (p,), sigma_hat (p,),
    Sigma_hat (p, p)) where the indicator is 1(xi'x <= b). Matches the sample
    covariance convention used by the legacy validator (``np.cov`` default ddof=1).
    """
    ind = (ksi_val @ path <= b).astype(float)      # (n2, p)
    P_hat = ind.mean(axis=0)
    Sigma_hat = np.cov(ind.T)                        # (p, p), ddof=1
    Sigma_hat = np.atleast_2d(Sigma_hat)
    sigma_hat = np.sqrt(np.clip(np.diag(Sigma_hat), 0.0, None))
    return P_hat, sigma_hat, Sigma_hat


def gs_quantiles(sigma_hat, Sigma_hat, beta, rng, sim_num=2000):
    """Gaussian-supremum critical values (normalized qhat, unnormalized qhat_u).

    Uses only the non-degenerate (positive-variance) candidates, drawing the
    multivariate-normal simulation from the injected ``rng`` for reproducibility.
    """
    nz = sigma_hat > 0
    if not np.any(nz):
        return 0.0, 0.0
    S = Sigma_hat[np.ix_(nz, nz)]
    l = S.shape[0]
    Z = rng.multivariate_normal(np.zeros(l), S, sim_num)
    tmat = np.tile(sigma_hat[nz], (sim_num, 1))
    qhat = float(np.percentile(np.amax(Z / tmat, axis=1), (1 - beta) * 100))
    qhat_u = float(np.percentile(np.amax(Z, axis=1), (1 - beta) * 100))
    return qhat, qhat_u


def _pick(objectives, feasible_mask, fallback_score):
    """Lowest-objective feasible candidate, else argmax of the fallback score."""
    if np.any(feasible_mask):
        idx = np.flatnonzero(feasible_mask)
        return int(idx[np.argmin(objectives[idx])]), True
    return int(np.argmax(fallback_score)), False


def select_all(P_hat, sigma_hat, Sigma_hat, objectives, *, alpha, beta, n2, rng, sim_num=2000):
    """Return {method: (index, met_margin)} for NV, UG, NGS, UNGS."""
    P_hat = np.asarray(P_hat, float); sigma_hat = np.asarray(sigma_hat, float)
    objectives = np.asarray(objectives, float)
    z = float(norm.ppf(1 - beta))
    root_n2 = np.sqrt(n2)
    qhat, qhat_u = gs_quantiles(sigma_hat, Sigma_hat, beta, rng, sim_num)

    out = {}
    out["NV"] = _pick(objectives, P_hat >= 1 - alpha, P_hat)
    m_ug = z * sigma_hat / root_n2
    out["UG"] = _pick(objectives, P_hat >= 1 - alpha + m_ug, P_hat - m_ug)
    m_ngs = qhat * sigma_hat / root_n2
    out["NGS"] = _pick(objectives, P_hat >= 1 - alpha + m_ngs, P_hat - m_ngs)
    m_ungs = qhat_u / root_n2
    out["UNGS"] = _pick(objectives, P_hat >= 1 - alpha + m_ungs, P_hat - m_ungs)
    return out
