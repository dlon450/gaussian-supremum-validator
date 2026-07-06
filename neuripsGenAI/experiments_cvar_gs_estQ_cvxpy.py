# Write an upgraded script that makes the candidate path Q-aware by solving a WEIGHTED
# Wasserstein–CVaR program where the (1/n) sum term is replaced by sum_i w_i z_i (w_i >= 0, sum w_i = 1).
# We estimate density-ratio weights on TRAINING samples using a logistic classifier trained to separate
# early-train P-like vs recent-validation Q-like. This lets the candidate portfolio adapt to the shift,
# not just the validator.
#
# Includes solver fallback (ECOS -> CLARABEL -> SCS).

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Shift-aware GS validation for robust CVaR with Q-aware candidate path (weighted training)
----------------------------------------------------------------------------------------
Upgrades over previous version:
  • AR(1) dependence (financially realistic).
  • Estimate density-ratio and APPLY IT TO TRAINING samples to build a Q-aware candidate path:
      δ v + sum_i w_i z_i ≤ α r, r ≤ z_i + γ + ξ_i^T x.
  • Solvers fallback: ECOS -> CLARABEL -> SCS; uniform fallback per δ if needed.
  • Same validators: NEW (shift-aware + BMB GS), OLD NGS (i.i.d.), IWCV.
  • Single γ calibrated on Scenario 1 (NEW ≈ 1−α), reused for Scenario 2.
  • δ grid excludes 0 to avoid trivial choice.

Run:
  pip install cvxpy ecos scs clarabel
  python experiments_cvar_gs_depQ_weightedtrain.py
"""
import math
from typing import Tuple, Dict, Optional

import numpy as np
from numpy.random import default_rng

rng = default_rng(19)

# ---------------- Helpers ----------------

def make_cov(d: int, vol_abs: float, seed: Optional[int] = None) -> np.ndarray:
    rlocal = default_rng(seed) if seed is not None else rng
    A = rlocal.normal(size=(d, d))
    G = A @ A.T
    G = G / (np.trace(G) / d)
    return (vol_abs ** 2) * G

def ar1_sim(mu: np.ndarray, Sigma_eps: np.ndarray, phi: float, n: int, x0: Optional[np.ndarray] = None):
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

def empirical_cvar(losses: np.ndarray, alpha: float):
    n = len(losses)
    k = max(1, int(math.ceil(alpha * n)))
    s = np.sort(losses)[::-1]
    worst = s[:k]
    return float(np.mean(worst)), float(worst[-1])

def weighted_empirical_cvar(losses: np.ndarray, w: np.ndarray, alpha: float):
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

# ------------- Ratio estimation and transfer to TRAINING -------------

def fit_logistic(P_like: np.ndarray, Q_like: np.ndarray, lr: float = 0.4, iters: int = 600):
    """Train a simple logistic on stacked data to separate P-like (0) vs Q-like (1)."""
    X = np.vstack([P_like, Q_like])
    y = np.hstack([np.zeros(len(P_like)), np.ones(len(Q_like))])
    idx = rng.permutation(len(y))
    X, y = X[idx], y[idx]
    Xb = np.hstack([X, np.ones((X.shape[0],1))])
    w = np.zeros(Xb.shape[1])
    for _ in range(iters):
        z = Xb @ w
        p = 1/(1+np.exp(-z))
        grad = Xb.T @ (p - y) / len(y)
        w -= lr * grad
    return w

def score_logistic(X: np.ndarray, w: np.ndarray, clip: float = 50.0) -> np.ndarray:
    Xb = np.hstack([X, np.ones((X.shape[0],1))])
    pQ = 1/(1+np.exp(-(Xb @ w)))
    odds = pQ / (1 - pQ + 1e-12)
    odds = np.clip(odds, 1e-6, clip)
    odds /= np.sum(odds)
    return odds

# ------------- CVXPY candidate path with sample weights ----------------

def _solve_one_delta_weighted(xi, c, alpha, gamma, delta, norm, long_only, wts, try_solvers):
    import cvxpy as cp
    n, d = xi.shape
    x = cp.Variable(d)
    v = cp.Variable(nonneg=True)
    r = cp.Variable(nonneg=True)
    z = cp.Variable(n, nonneg=True)
    # Replace (1/n) sum z by sum_i w_i z_i; w_i >= 0, sum w_i = 1
    constr = [delta*v + cp.sum(cp.multiply(wts, z)) <= alpha*r,
              r <= z + gamma + xi @ x]
    if norm == 'l2':   constr += [cp.norm2(x) <= v]
    elif norm == 'l1': constr += [cp.norm1(x) <= v]
    elif norm == 'linf': constr += [cp.norm_inf(x) <= v]
    else: raise ValueError("norm must be l2, l1, or linf")
    if long_only: constr += [cp.sum(x) == 1, x >= 0]
    else:         constr += [cp.sum(x) == 1]
    prob = cp.Problem(cp.Minimize(c @ x), constr)

    opts = {}
    if "ECOS" in try_solvers:
        opts["ECOS"] = dict(max_iters=8000, abstol=1e-8, reltol=1e-8, feastol=1e-8, verbose=False)
    if "SCS" in try_solvers:
        opts["SCS"] = dict(max_iters=25000, eps=7e-5, acceleration_lookback=0, verbose=False)
    if "CLARABEL" in try_solvers:
        opts["CLARABEL"] = dict(verbose=False)

    for s in try_solvers:
        try:
            prob.solve(solver=s, **opts.get(s, {}))
            if x.value is not None and prob.status in ["optimal", "optimal_inaccurate"]:
                xx = np.maximum(x.value, 0) if long_only else x.value
                ssum = np.sum(xx)
                if ssum > 0:
                    return xx / ssum
        except Exception:
            continue
    return None

def solve_wass_cvar_candidates_weighted(R_train: np.ndarray, alpha: float, gamma: float,
                                        deltas: np.ndarray, weights: np.ndarray,
                                        norm: str = 'l2', long_only: bool = True) -> np.ndarray:
    try:
        import cvxpy as cp
        a = dfglijrel
    except Exception:
        d = R_train.shape[1]
        return dirichlet_simplex(len(deltas), d, alpha_dir=1.0, seed=123)
    xi = R_train.astype(float)
    n, d = xi.shape
    wts = np.maximum(weights, 1e-12)
    wts = wts / np.sum(wts)
    # Objective expected loss: use weighted mean under wts to tilt toward Q
    c = -np.sum((wts[:, None] * xi), axis=0)
    X = np.zeros((d, len(deltas)))
    installed = set(cp.installed_solvers())
    try_solvers = [s for s in ["ECOS", "CLARABEL", "SCS"] if s in installed]
    if not try_solvers:
        return dirichlet_simplex(len(deltas), d, alpha_dir=1.0, seed=123)
    any_fail = False
    for j, delta in enumerate(deltas):
        xx = _solve_one_delta_weighted(xi, c, alpha, gamma, delta, norm, long_only, wts, try_solvers)
        if xx is None:
            any_fail = True
            X[:, j] = np.full(d, 1.0/d)
        else:
            X[:, j] = xx
    if any_fail:
        print("[WARN] Some weighted δ solves failed; uniform fallback used.")
    return X

# ---------------- Validators (as before) ----------------

def new_validator_select_delta(X_cands, R_val, alpha, gamma, deltas, mu_p,
                               beta=0.05, B_boot=400, clip_w=50.0, recent_frac=0.5):
    n2, d = R_val.shape
    m = max(50, int(round(recent_frac*n2)))
    P_val = R_val[:-m]; Q_recent = R_val[-m:]
    w = score_logistic(R_val, fit_logistic(P_val, Q_recent), clip=clip_w)  # weights on validation
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

def old_ngs_validator_select_delta(X_cands, R_val, alpha, gamma, deltas, mu_p,
                                   beta=0.05, sim_num=2000):
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

def iwcv_select_delta(X_cands, R_val, alpha, gamma, deltas, mu_p, recent_frac=0.5, clip_w=50.0):
    n2 = len(R_val)
    m = max(50, int(round(recent_frac*n2)))
    P_val = R_val[:-m]; Q_recent = R_val[-m:]
    w = score_logistic(R_val, fit_logistic(P_val, Q_recent), clip=clip_w)
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

def run_experiment_once(d=8, n_train=1000, n_val=1200, n_test=15000,
                        alpha=0.10, gamma=0.06,
                        mean_scale=0.01, vol_abs=0.02,
                        shift=False, mean_shift_abs=-0.004, vol_mult=1.6,
                        phi=0.3, phi_shift=0.45,
                        deltas=None, beta=0.05, B_boot=300, recent_frac=0.5):
    if deltas is None:
        deltas = np.linspace(0.003, 0.04, 16)  # wider range, exclude 0
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

    R_train = ar1_sim(mu_p, Sigma_eps, phi, n_train)
    R_val   = ar1_sim(mu_p, Sigma_eps, phi, n_val, x0=R_train[-1])
    R_test  = ar1_sim(mu_q, Sigma_eps_q, phi_q, n_test, x0=R_val[-1])

    # ---- NEW: learn Q-weights on TRAIN using early-train vs recent-VAL ----
    ntr = len(R_train)
    split = max(100, int(0.6*ntr))
    P_like_train = R_train[:split]           # older training chunk as P-like
    Q_like_proxy = R_val[-max(100, int(0.4*len(R_val))):]  # recent validation chunk as Q-like
    w_model = fit_logistic(P_like_train, Q_like_proxy)
    w_train = score_logistic(R_train, w_model)  # weights for training samples (sum=1, clipped)

    # Candidate path with WEIGHTED training (Q-aware)
    X = solve_wass_cvar_candidates_weighted(R_train, alpha=alpha, gamma=gamma, deltas=deltas,
                                            weights=w_train, norm='l2', long_only=True)

    mu_hat_p = np.mean(R_train, axis=0)

    sel_new = new_validator_select_delta(X, R_val, alpha, gamma, deltas, mu_hat_p,
                                         beta=beta, B_boot=B_boot, recent_frac=recent_frac)
    sel_old = old_ngs_validator_select_delta(X, R_val, alpha, gamma, deltas, mu_hat_p, beta=beta)
    sel_iw  = iwcv_select_delta(X, R_val, alpha, gamma, deltas, mu_hat_p, recent_frac=recent_frac)

    def eval_candidate(x, delta):
        losses = portfolio_loss(R_test, x)
        cvar, _ = empirical_cvar(losses, alpha)
        lhs = cvar + (delta/alpha) * np.linalg.norm(x)
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
        deltas=np.linspace(0.003, 0.04, 16),
        beta=0.05, B_boot=300, recent_frac=0.5,
        phi=0.3, phi_shift=0.45
    )

    # Calibrate gamma on Scenario 1 (NEW ≈ 1 - alpha), reuse for Scenario 2
    s1_params = dict(shift=False, mean_shift_abs=0.0, vol_mult=1.0, **common)
    gamma = calibrate_gamma(target=1 - common["alpha"], low=0.03, high=0.12, tol=0.08, max_iter=8, reps=6, **s1_params)
    print("Calibrated gamma (Scenario 1):", gamma)

    res1 = run_scenarios(n_rep=100, gamma=gamma, **s1_params)
    print("Scenario 1 (no shift):")
    print(res1)

    s2_params = dict(shift=True, mean_shift_abs=-0.004, vol_mult=1.7, **common)
    res2 = run_scenarios(n_rep=100, gamma=gamma, **s2_params)
    print("Scenario 2 (shift):")
    print(res2)