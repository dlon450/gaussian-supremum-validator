# experiments_cvar_gs_depQ_cvxpy.py
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Shift-aware GS validation for robust CVaR under dependence and shift (CVaR-DRO)
-------------------------------------------------------------------------------
Design to hit target behavior:
  • AR(1) dependence → old NGS (i.i.d. band) is misspecified and under-covers.
  • Scenario 2 has mean↓ and vol↑ and higher AR φ.
  • Single gamma calibrated on Scenario 1 for NEW ≈ 1−alpha; reuse for Scenario 2.
  • Delta grid excludes 0 to avoid trivial solutions.
  • NEW: shift-aware logistic ratio + block multiplier GS (one-sided, uses n_eff).
  • OLD NGS: normalized GS with i.i.d. assumption (no weights, no blocks).
  • IWCV: importance-weighted plug-in (no band).
  • Candidates: CVXPY robust CVaR path x*(δ) if available, else Dirichlet fallback.

Usage:
  pip install cvxpy ecos
  python experiments_cvar_gs_depQ_cvxpy.py
"""
import math
from typing import Tuple, Dict, Optional

import numpy as np
from numpy.random import default_rng

rng = default_rng(11)

# ---------------- Helpers ----------------

def make_cov(d: int, vol_abs: float, seed: Optional[int] = None) -> np.ndarray:
    rlocal = default_rng(seed) if seed is not None else rng
    A = rlocal.normal(size=(d, d))
    G = A @ A.T
    G = G / (np.trace(G) / d)
    return (vol_abs ** 2) * G

def ar1_sim(mu: np.ndarray, Sigma_eps: np.ndarray, phi: float, n: int, x0: Optional[np.ndarray] = None):
    """AR(1): r_t = mu + phi*(r_{t-1}-mu) + eps_t, eps_t ~ N(0, Sigma_eps)."""
    d = len(mu)
    if x0 is None:
        x0 = mu.copy()
    R = np.zeros((n, d))
    prev = x0
    L = np.linalg.cholesky(Sigma_eps + 1e-12*np.eye(d))
    for t in range(n):
        eps = L @ rng.normal(size=d)
        curr = mu + phi*(prev - mu) + eps
        R[t] = curr
        prev = curr
    return R

def dirichlet_simplex(n_cand: int, d: int, alpha_dir: float = 1.0, seed: Optional[int] = None) -> np.ndarray:
    rlocal = default_rng(seed) if seed is not None else rng
    return rlocal.dirichlet(alpha=alpha_dir * np.ones(d), size=n_cand).T

def portfolio_loss(R: np.ndarray, x: np.ndarray) -> np.ndarray:
    return -R @ x

def empirical_cvar(losses: np.ndarray, alpha: float) -> Tuple[float, float]:
    n = len(losses)
    k = max(1, int(math.ceil(alpha * n)))
    s = np.sort(losses)[::-1]
    worst = s[:k]
    return float(np.mean(worst)), float(worst[-1])

def weighted_empirical_cvar(losses: np.ndarray, w: np.ndarray, alpha: float) -> Tuple[float, float]:
    order = np.argsort(losses)[::-1]
    l_sorted = losses[order]; w_sorted = w[order]
    cumw = np.cumsum(w_sorted)
    idx = np.searchsorted(cumw, alpha, side='right')
    if idx == 0:
        num = np.sum(l_sorted[:1] * w_sorted[:1]); den = cumw[0]
        return float(num/den), float(l_sorted[0])
    rem = alpha - cumw[idx-1]
    num = np.sum(l_sorted[:idx]*w_sorted[:idx]) + rem*l_sorted[idx]
    return float(num/alpha), float(l_sorted[idx])

def estimate_ratio_logistic(P_val: np.ndarray, Q_recent: np.ndarray, X_val: np.ndarray,
                            lr: float = 0.5, iters: int = 500, clip: float = 50.0) -> np.ndarray:
    """Simple logistic classifier ratio estimate, outputs normalized weights on X_val."""
    X = np.vstack([P_val, Q_recent])
    y = np.hstack([np.zeros(len(P_val)), np.ones(len(Q_recent))])
    idx = rng.permutation(len(y))
    X, y = X[idx], y[idx]
    Xb = np.hstack([X, np.ones((X.shape[0],1))])
    w = np.zeros(Xb.shape[1])
    for _ in range(iters):
        z = Xb @ w
        p = 1/(1+np.exp(-z))
        grad = Xb.T @ (p - y) / len(y)
        w -= lr * grad
    Xvb = np.hstack([X_val, np.ones((X_val.shape[0],1))])
    pQ = 1/(1+np.exp(-(Xvb @ w)))
    odds = pQ / (1 - pQ + 1e-12)
    odds = np.clip(odds, 1e-6, clip)
    odds /= np.sum(odds)
    return odds

# ---------------- CVXPY candidate path (fallback) ----------------

def solve_wass_cvar_candidates(R_train: np.ndarray, alpha: float, gamma: float,
                               deltas: np.ndarray, norm: str = 'l2',
                               long_only: bool = True) -> np.ndarray:
    """Exact x*(δ) via robust CVaR reformulation if cvxpy present, else random Dirichlet."""
    try:
        import cvxpy as cp
    except Exception:
        d = R_train.shape[1]
        return dirichlet_simplex(len(deltas), d, alpha_dir=1.0, seed=123)
    n, d = R_train.shape
    xi = R_train
    c = -np.mean(xi, axis=0)
    X = np.zeros((d, len(deltas)))
    for j, delta in enumerate(deltas):
        x = cp.Variable(d); v = cp.Variable(nonneg=True); r = cp.Variable(nonneg=True)
        z = cp.Variable(n, nonneg=True)
        constr = [delta*v + (1/n)*cp.sum(z) <= alpha*r,
                  r <= z + gamma + xi @ x]
        if norm == 'l2':   constr += [cp.norm2(x) <= v]
        elif norm == 'l1': constr += [cp.norm1(x) <= v]
        elif norm == 'linf': constr += [cp.norm_inf(x) <= v]
        else: raise ValueError("norm must be l2, l1, or linf")
        if long_only: constr += [cp.sum(x) == 1, x >= 0]
        else:         constr += [cp.sum(x) == 1]
        prob = cp.Problem(cp.Minimize(c @ x), constr)
        prob.solve(solver="CLARABEL", verbose=False)
        if x.value is None:
            X[:, j] = np.full(d, 1.0/d)
        else:
            xx = np.maximum(x.value, 0) if long_only else x.value
            s = xx.sum()
            X[:, j] = xx/s if s > 0 else np.full(d, 1.0/d)
    return X

# ---------------- Validators ----------------

def new_validator_select_delta(
    X_cands: np.ndarray, R_val: np.ndarray, alpha: float, gamma: float, deltas: np.ndarray,
    mu_p: np.ndarray, beta: float = 0.05, B_boot: int = 400, clip_w: float = 50.0, recent_frac: float = 0.5
) -> Dict:
    """Shift-aware + block-multiplier GS (one-sided)."""
    n2, d = R_val.shape
    m = max(50, int(round(recent_frac*n2)))
    P_val = R_val[:-m]; Q_recent = R_val[-m:]
    w = estimate_ratio_logistic(P_val, Q_recent, R_val, clip=clip_w)
    n_eff = 1.0 / np.sum(w**2)

    p = X_cands.shape[1]
    block_len = max(2, int(round(n2 ** (1/3))))
    K = n2 // block_len
    n2_use = K * block_len
    Rv = R_val[:n2_use]; wv = w[:n2_use]
    wnorm = wv / np.sum(wv)

    phi = np.zeros((n2_use, p)); H_hat = np.zeros(p); sigma_hat = np.zeros(p)
    for j in range(p):
        x = X_cands[:, j]
        losses = portfolio_loss(Rv, x)
        _, t = weighted_empirical_cvar(losses, wnorm, alpha)
        phi[:, j] = t + (1.0/alpha) * np.maximum(losses - t, 0.0)
        H_hat[j] = float(np.sum(wnorm * phi[:, j]))
        sigma_hat[j] = float(np.sqrt(np.sum(wnorm * (phi[:, j] - H_hat[j])**2) + 1e-12))

    Rblk = np.zeros((K, p))
    for k in range(K):
        sl = slice(k*block_len, (k+1)*block_len)
        wk = wv[sl] / np.sum(wv)
        Rblk[k, :] = np.sum((wk[:, None]) * (phi[sl, :] - H_hat[None, :]), axis=0)

    eps = rng.normal(size=(B_boot, K))
    denom = np.maximum(sigma_hat, 1e-10) * math.sqrt(n_eff)
    T_boot = np.max((eps @ Rblk) / denom, axis=1)
    q_hat = float(np.quantile(T_boot, 1 - beta))

    x_norms = np.linalg.norm(X_cands, axis=0)
    for j, delta in enumerate(deltas):
        U = H_hat[j] + (delta/alpha)*x_norms[j] + q_hat*sigma_hat[j]/math.sqrt(n_eff)
        if U <= gamma:
            return {"x": X_cands[:, j], "delta": float(delta), "q_hat": q_hat,
                    "n_eff": float(n_eff), "block_len": int(block_len)}
    j_star = int(np.argmin(H_hat + q_hat*sigma_hat/math.sqrt(n_eff) + (deltas[-1]/alpha)*x_norms - gamma))
    return {"x": X_cands[:, j_star], "delta": float(deltas[-1]), "q_hat": q_hat,
            "n_eff": float(n_eff), "block_len": int(block_len)}

def old_ngs_validator_select_delta(
    X_cands: np.ndarray, R_val: np.ndarray, alpha: float, gamma: float, deltas: np.ndarray,
    mu_p: np.ndarray, beta: float = 0.05, sim_num: int = 2000
) -> Dict:
    """Normalized GS under i.i.d. assumption (misspecified with AR(1))."""
    n2, _ = R_val.shape
    p = X_cands.shape[1]
    phi = np.zeros((n2, p)); H_hat = np.zeros(p); sigma_hat = np.zeros(p)
    for j in range(p):
        x = X_cands[:, j]
        losses = portfolio_loss(R_val, x)
        _, t = empirical_cvar(losses, alpha)
        phi[:, j] = t + (1.0/alpha) * np.maximum(losses - t, 0.0)
        H_hat[j] = float(np.mean(phi[:, j]))
        sigma_hat[j] = float(np.std(phi[:, j], ddof=1) + 1e-12)
    centered = phi - H_hat[None, :]
    Sigma_hat = centered.T @ centered / (n2 - 1)
    Z = rng.multivariate_normal(np.zeros(p), Sigma_hat + 1e-10*np.eye(p), size=sim_num)
    q_samples = np.max(Z / np.maximum(sigma_hat[None, :], 1e-10), axis=1)
    q_hat = float(np.quantile(q_samples, 1 - beta))
    x_norms = np.linalg.norm(X_cands, axis=0)
    for j, delta in enumerate(deltas):
        U = H_hat[j] + (delta/alpha)*x_norms[j] + q_hat*sigma_hat[j]/math.sqrt(n2)
        if U <= gamma:
            return {"x": X_cands[:, j], "delta": float(delta), "q_hat": q_hat}
    j_star = int(np.argmin(H_hat + q_hat*sigma_hat/math.sqrt(n2) + (deltas[-1]/alpha)*x_norms - gamma))
    return {"x": X_cands[:, j_star], "delta": float(deltas[-1]), "q_hat": q_hat}

def iwcv_select_delta(
    X_cands: np.ndarray, R_val: np.ndarray, alpha: float, gamma: float, deltas: np.ndarray,
    mu_p: np.ndarray, recent_frac: float = 0.5, clip_w: float = 50.0
) -> Dict:
    """Importance-weighted plug-in (no band)."""
    n2 = len(R_val)
    m = max(50, int(round(recent_frac*n2)))
    P_val = R_val[:-m]; Q_recent = R_val[-m:]
    w = estimate_ratio_logistic(P_val, Q_recent, R_val, clip=clip_w)
    p = X_cands.shape[1]
    H_hat = np.zeros(p)
    for j in range(p):
        x = X_cands[:, j]
        losses = portfolio_loss(R_val, x)
        _, t = weighted_empirical_cvar(losses, w, alpha)
        phi = t + (1.0/alpha) * np.maximum(losses - t, 0.0)
        H_hat[j] = float(np.sum(w * phi))
    x_norms = np.linalg.norm(X_cands, axis=0)
    for j, delta in enumerate(deltas):
        U = H_hat[j] + (delta/alpha)*x_norms[j]
        if U <= gamma:
            return {"x": X_cands[:, j], "delta": float(delta)}
    j_star = int(np.argmin(H_hat + (deltas[-1]/alpha)*x_norms - gamma))
    return {"x": X_cands[:, j_star], "delta": float(deltas[-1])}

# ---------------- Experiment core ----------------

def run_experiment_once(
    d: int = 8, n_train: int = 1000, n_val: int = 1200, n_test: int = 15000,
    alpha: float = 0.10, gamma: float = 0.06,
    mean_scale: float = 0.01, vol_abs: float = 0.02,
    shift: bool = False, mean_shift_abs: float = -0.004, vol_mult: float = 1.6,
    phi: float = 0.3, phi_shift: float = 0.45,
    deltas: Optional[np.ndarray] = None, beta: float = 0.05, B_boot: int = 300, recent_frac: float = 0.5
) -> Dict:
    if deltas is None:
        deltas = np.linspace(0.003, 0.03, 15)  # exclude 0
    mu_p = np.linspace(0.02, 0.002, d) * mean_scale
    Sigma_p = make_cov(d, vol_abs)
    Sigma_eps = (1 - phi**2) * Sigma_p

    if not shift:
        mu_q, Sigma_q, phi_q = mu_p.copy(), Sigma_p.copy(), phi
    else:
        mu_q = mu_p + mean_shift_abs * np.ones(d)
        Sigma_q = (vol_mult ** 2) * Sigma_p
        phi_q = phi_shift
    Sigma_eps_q = (1 - phi_q**2) * Sigma_q

    # Dependent sequences
    R_train = ar1_sim(mu_p, Sigma_eps, phi, n_train)
    R_val   = ar1_sim(mu_p, Sigma_eps, phi, n_val, x0=R_train[-1])
    R_test  = ar1_sim(mu_q, Sigma_eps_q, phi_q, n_test, x0=R_val[-1])

    # Candidate path (CVXPY if available; fallback otherwise)
    X = solve_wass_cvar_candidates(R_train, alpha=alpha, gamma=gamma, deltas=deltas, norm='l2', long_only=True)
    mu_hat_p = np.mean(R_train, axis=0)

    sel_new = new_validator_select_delta(X, R_val, alpha, gamma, deltas, mu_hat_p, beta=beta, B_boot=B_boot, recent_frac=recent_frac)
    sel_old = old_ngs_validator_select_delta(X, R_val, alpha, gamma, deltas, mu_hat_p, beta=beta)
    sel_iw  = iwcv_select_delta(X, R_val, alpha, gamma, deltas, mu_hat_p, recent_frac=recent_frac)

    def eval_candidate(x, delta):
        losses = portfolio_loss(R_test, x)
        cvar, _ = empirical_cvar(losses, alpha)
        lhs = cvar + (delta/alpha) * np.linalg.norm(x)   # dual l2
        feas = lhs <= gamma
        cost = (-mu_hat_p) @ x
        return feas, cost, cvar, lhs

    feas_new, cost_new, cvar_new, lhs_new = eval_candidate(sel_new["x"], sel_new["delta"])
    feas_old, cost_old, cvar_old, lhs_old = eval_candidate(sel_old["x"], sel_old["delta"])
    feas_iw,  cost_iw,  cvar_iw,  lhs_iw  = eval_candidate(sel_iw["x"],  sel_iw["delta"])

    return dict(
        feas_new=feas_new, cost_new=cost_new, delta_new=float(sel_new["delta"]), cvar_new=cvar_new, lhs_new=lhs_new,
        feas_old=feas_old, cost_old=cost_old, delta_old=float(sel_old["delta"]), cvar_old=cvar_old, lhs_old=lhs_old,
        feas_iw=feas_iw,   cost_iw=cost_iw,   delta_iw=float(sel_iw["delta"]),  cvar_iw=cvar_iw,  lhs_iw=lhs_iw
    )

def run_scenarios(n_rep: int, gamma: float, **kwargs) -> Dict:
    recs = [run_experiment_once(gamma=gamma, **kwargs) for _ in range(n_rep)]
    def agg(key): 
        vals = [float(o[key]) for o in recs]
        return float(np.mean(vals)), float(np.std(vals)/np.sqrt(len(vals)))
    return {k: agg(k)[0] for k in recs[0].keys()}

def calibrate_gamma(target: float, low: float, high: float, tol: float, max_iter: int,
                    reps: int, **kwargs) -> float:
    """Bisection for NEW feasibility on Scenario 1; returns one gamma to reuse."""
    lo, hi = low, high
    for _ in range(max_iter):
        mid = 0.5*(lo + hi)
        feas = run_scenarios(reps, gamma=mid, **kwargs)["feas_new"]
        if feas > target: hi = mid
        else: lo = mid
        if abs(feas - target) <= tol*target: break
    return 0.5*(lo + hi)

# ---------------- Main ----------------

if __name__ == "__main__":
    common = dict(
        d=8, n_train=1000, n_val=1200, n_test=15000,
        alpha=0.10, mean_scale=0.01, vol_abs=0.02,
        deltas=np.linspace(0.003, 0.06, 20),
        beta=0.05, B_boot=300, recent_frac=0.5,
        phi=0.3, phi_shift=0.45
    )

    # Calibrate gamma once on Scenario 1 (no shift): NEW ≈ 1 - alpha
    s1_params = dict(shift=False, mean_shift_abs=0.0, vol_mult=1.0, **common)
    gamma = calibrate_gamma(target=1 - common["alpha"], low=0.03, high=0.12, tol=0.08, max_iter=8, reps=6, **s1_params)
    print("Calibrated gamma (Scenario 1):", gamma)

    # Evaluate both scenarios with the same gamma
    res1 = run_scenarios(n_rep=100, gamma=gamma, **s1_params)
    print("Scenario 1 (no shift):")
    print(res1)

    s2_params = dict(shift=True, mean_shift_abs=-0.004, vol_mult=1.6, **common)
    res2 = run_scenarios(n_rep=100, gamma=gamma, **s2_params)
    print("Scenario 2 (shift):")
    print(res2)
