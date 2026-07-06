#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fair-gamma CVaR–DRO with Clarabel; shift-aware candidates; FIXED selector fallback
-------------------------------------------------------------------------------
- Gamma fixed ex-ante from equal-weight CVaR on Scenario-1 train; SAME γ in both scenarios.
- Candidates x*(δ) solved via CVXPY (Clarabel -> ECOS -> SCS).
- Candidate path is shift-aware (train weights blended with uniform).
- Validators: NEW (shift-aware + BMB GS), OLD (i.i.d. NGS), IW (plug-in).
- IMPORTANT FIX: when no δ is feasible, return the (x, δ) that MINIMIZES U(δ) instead of forcing δ_max.
"""

from typing import Optional, Dict, Tuple
import numpy as np
import numpy.random as npr
import math

_rng = npr.default_rng(123)

# ---------------- Data gen ----------------
def make_cov(d: int, vol_abs: float, seed: Optional[int] = None) -> np.ndarray:
    rng = _rng if seed is None else npr.default_rng(seed)
    A = rng.normal(size=(d, d))
    G = A @ A.T
    G /= (np.trace(G) / d)
    return (vol_abs ** 2) * G

def ar1_sim(mu: np.ndarray, Sigma_eps: np.ndarray, phi: float, n: int, x0: Optional[np.ndarray] = None) -> np.ndarray:
    d = len(mu)
    if x0 is None: x0 = mu.copy()
    R = np.zeros((n, d)); prev = x0
    L = np.linalg.cholesky(Sigma_eps + 1e-12*np.eye(d))
    for t in range(n):
        eps = L @ _rng.normal(size=d)
        curr = mu + phi*(prev - mu) + eps
        R[t] = curr; prev = curr
    return R

# ---------------- Risk ----------------
def portfolio_loss(R: np.ndarray, x: np.ndarray) -> np.ndarray:
    return -R @ x

def empirical_cvar(losses: np.ndarray, alpha: float):
    n = len(losses); k = max(1, int(math.ceil(alpha*n)))
    s = np.sort(losses)[::-1]; tail = s[:k]
    return float(np.mean(tail)), float(tail[-1])

def weighted_empirical_cvar(losses: np.ndarray, w: np.ndarray, alpha: float):
    order = np.argsort(losses)[::-1]
    l_sorted = losses[order]; w_sorted = w[order]
    cumw = np.cumsum(w_sorted)
    idx = np.searchsorted(cumw, alpha, side='right')
    if idx == 0:
        num = np.sum(l_sorted[:1]*w_sorted[:1]); den = cumw[0]
        return float(num/den), float(l_sorted[0])
    rem = alpha - cumw[idx-1]
    num = np.sum(l_sorted[:idx]*w_sorted[:idx]) + rem*l_sorted[idx]
    return float(num/alpha), float(l_sorted[idx])

# ---------------- Density ratio ----------------
def fit_logistic(P_like: np.ndarray, Q_like: np.ndarray, lr: float = 0.4, iters: int = 600) -> np.ndarray:
    X = np.vstack([P_like, Q_like]); y = np.hstack([np.zeros(len(P_like)), np.ones(len(Q_like))])
    idx = _rng.permutation(len(y)); X, y = X[idx], y[idx]
    Xb = np.hstack([X, np.ones((X.shape[0],1))]); w = np.zeros(Xb.shape[1])
    for _ in range(iters):
        z = Xb @ w; p = 1/(1+np.exp(-z)); grad = Xb.T @ (p - y) / len(y); w -= lr * grad
    return w

def score_logistic(X: np.ndarray, w: np.ndarray, clip: float = 50.0) -> np.ndarray:
    Xb = np.hstack([X, np.ones((X.shape[0],1))]); pQ = 1/(1+np.exp(-(Xb @ w)))
    odds = pQ / (1 - pQ + 1e-12); odds = np.clip(odds, 1e-6, clip); odds /= np.sum(odds)
    return odds

# ---------------- CVXPY candidates (Clarabel first) ----------------
def _solve_one_delta_weighted(Xi: np.ndarray, c: np.ndarray, alpha: float, gamma: float,
                              delta: float, wts: np.ndarray, norm: str,
                              long_only: bool, try_solvers) -> Optional[np.ndarray]:
    import cvxpy as cp
    n, d = Xi.shape
    x = cp.Variable(d); v = cp.Variable(nonneg=True); r = cp.Variable(nonneg=True); z = cp.Variable(n, nonneg=True)
    constr = [delta*v + cp.sum(cp.multiply(wts, z)) <= alpha*r, r <= z + gamma + Xi @ x]
    if norm == 'l2':   constr += [cp.norm2(x) <= v]
    elif norm == 'l1': constr += [cp.norm1(x) <= v]
    elif norm == 'linf': constr += [cp.norm_inf(x) <= v]
    else: raise ValueError("norm must be l2/l1/linf")
    constr += [cp.sum(x) == 1, x >= 0] if long_only else [cp.sum(x) == 1]
    prob = cp.Problem(cp.Minimize(c @ x), constr)
    opts = {"CLARABEL": dict(verbose=False),
            "ECOS": dict(max_iters=8000, abstol=1e-8, reltol=1e-8, feastol=1e-8, verbose=False),
            "SCS": dict(max_iters=25000, eps=7e-5, acceleration_lookback=0, verbose=False)}
    for s in try_solvers:
        try:
            prob.solve(solver=s, **opts.get(s, {}))
            if x.value is not None and prob.status in ["optimal","optimal_inaccurate"]:
                xx = np.maximum(x.value,0) if long_only else x.value
                ssum = float(np.sum(xx))
                if ssum > 0: return xx/ssum
        except Exception:
            continue
    return None

def solve_wass_cvar_candidates_weighted(R_train: np.ndarray, alpha: float, gamma: float,
                                        deltas: np.ndarray, weights: np.ndarray,
                                        norm: str = 'l2', long_only: bool = True) -> np.ndarray:
    import cvxpy as cp
    Xi = R_train.astype(float); n, d = Xi.shape
    wts = np.maximum(weights, 1e-12); wts /= np.sum(wts)
    c = - (Xi.T @ wts)                              # Q-aware mean loss
    X = np.zeros((d, len(deltas)))
    try_solvers = [s for s in ["CLARABEL","ECOS","SCS"] if s in cp.installed_solvers()]
    if not try_solvers: return np.full((d, len(deltas)), 1.0/d)
    for j, delta in enumerate(deltas):
        xx = _solve_one_delta_weighted(Xi, c, alpha, gamma, float(delta), wts, norm, long_only, try_solvers)
        X[:, j] = xx if xx is not None else np.full(d, 1.0/d)
    return X

# ---------------- Validators (with FIXED fallback) ----------------
def _pick_minU(feasible_mask, U_vals, X_cands, deltas):
    if np.any(feasible_mask):
        j = int(np.argmin(U_vals[feasible_mask]))
        idx = np.where(feasible_mask)[0][j]
        return {"x": X_cands[:, idx], "delta": float(deltas[idx])}
    # no feasible: pick j with MIN U and return its δ (NOT δ_max!)
    j = int(np.argmin(U_vals))
    return {"x": X_cands[:, j], "delta": float(deltas[j])}

def new_validator_select_delta(X_cands: np.ndarray, R_val: np.ndarray, alpha: float, gamma: float,
                               deltas: np.ndarray, beta: float = 0.05, B_boot: int = 600,
                               recent_frac: float = 0.6, clip_w: float = 50.0) -> Dict:
    n2, d = R_val.shape
    m = max(50, int(round(recent_frac*n2)))
    P_val = R_val[:-m]; Q_recent = R_val[-m:]
    w = score_logistic(R_val, fit_logistic(P_val, Q_recent), clip=clip_w)
    n_eff = 1.0 / np.sum(w**2)

    p = X_cands.shape[1]
    phi = np.zeros((n2, p)); H_hat = np.zeros(p); sigma_hat = np.zeros(p)
    for j in range(p):
        x = X_cands[:, j]; losses = portfolio_loss(R_val, x)
        _, t = weighted_empirical_cvar(losses, w, alpha)
        phi[:, j] = t + (1.0/alpha)*np.maximum(losses - t, 0.0)
        H_hat[j] = float(np.sum(w * phi[:, j]))
        sigma_hat[j] = float(np.sqrt(np.sum(w * (phi[:, j] - H_hat[j])**2) + 1e-12))

    # BMB
    block_len = max(2, int(round(n2 ** (1/3)))); K = n2 // block_len; n2_use = K*block_len
    Rblk = np.zeros((K, p))
    for k in range(K):
        sl = slice(k*block_len, (k+1)*block_len)
        wk = w[sl] / np.sum(w[sl])
        Rblk[k, :] = np.sum((wk[:,None]) * (phi[sl,:] - H_hat[None,:]), axis=0)
    eps = _rng.normal(size=(B_boot, K))
    denom = np.maximum(sigma_hat, 1e-10) * math.sqrt(n_eff)
    T_boot = np.max((eps @ Rblk) / denom, axis=1)
    q_hat = float(np.quantile(T_boot, 1 - beta))

    x_norms = np.linalg.norm(X_cands, axis=0)     # dual l2
    U = H_hat + q_hat*sigma_hat/np.sqrt(n_eff) + (deltas/alpha)*x_norms
    feas = U <= gamma
    return _pick_minU(feas, U, X_cands, deltas)

def old_ngs_validator_select_delta(X_cands: np.ndarray, R_val: np.ndarray, alpha: float, gamma: float,
                                   deltas: np.ndarray, beta: float = 0.05) -> Dict:
    n2, _ = R_val.shape; p = X_cands.shape[1]
    phi = np.zeros((n2, p)); H_hat = np.zeros(p); sigma_hat = np.zeros(p)
    for j in range(p):
        x = X_cands[:, j]; losses = portfolio_loss(R_val, x)
        _, t = empirical_cvar(losses, alpha)
        phi[:, j] = t + (1.0/alpha)*np.maximum(losses - t, 0.0)
        H_hat[j] = float(np.mean(phi[:, j])); sigma_hat[j] = float(np.std(phi[:, j], ddof=1) + 1e-12)
    centered = phi - H_hat[None,:]
    Sigma_hat = centered.T @ centered / (n2 - 1)
    Z = _rng.multivariate_normal(np.zeros(p), Sigma_hat + 1e-10*np.eye(p), size=2000)
    q_samples = np.max(Z / np.maximum(sigma_hat[None,:], 1e-10), axis=1)
    q_hat = float(np.quantile(q_samples, 1 - beta))
    x_norms = np.linalg.norm(X_cands, axis=0)
    U = H_hat + q_hat*sigma_hat/np.sqrt(n2) + (deltas/alpha)*x_norms
    feas = U <= gamma
    return _pick_minU(feas, U, X_cands, deltas)

def iwcv_select_delta(X_cands: np.ndarray, R_val: np.ndarray, alpha: float, gamma: float,
                      deltas: np.ndarray, recent_frac: float = 0.6, clip_w: float = 50.0) -> Dict:
    n2 = len(R_val); m = max(50, int(round(recent_frac*n2)))
    P_val = R_val[:-m]; Q_recent = R_val[-m:]
    w = score_logistic(R_val, fit_logistic(P_val, Q_recent), clip=clip_w)
    p = X_cands.shape[1]
    H_hat = np.zeros(p)
    for j in range(p):
        x = X_cands[:, j]; losses = portfolio_loss(R_val, x)
        _, t = weighted_empirical_cvar(losses, w, alpha)
        phi = t + (1.0/alpha)*np.maximum(losses - t, 0.0)
        H_hat[j] = float(np.sum(w * phi))
    x_norms = np.linalg.norm(X_cands, axis=0)
    U = H_hat + (deltas/alpha)*x_norms
    feas = U <= gamma
    return _pick_minU(feas, U, X_cands, deltas)

# ---------------- Train weights for candidates ----------------
def estimate_train_weights_for_candidates(R_train: np.ndarray, R_val: np.ndarray,
                                          recent_frac: float = 0.6, blend_lambda: float = 0.3,
                                          clip: float = 50.0) -> np.ndarray:
    ntr = len(R_train); split = max(100, int(0.6*ntr))
    P_like_train = R_train[:split]
    m = max(50, int(round(recent_frac*len(R_val))))
    Q_like_proxy = R_val[-m:]
    w_model = fit_logistic(P_like_train, Q_like_proxy)
    w_train = score_logistic(R_train, w_model, clip=clip)
    w_blend = (1.0 - blend_lambda) * (np.ones(ntr)/ntr) + blend_lambda * w_train
    w_blend /= np.sum(w_blend)
    return w_blend

# ---------------- Core experiment ----------------
def solve_candidates(R_train, R_val, alpha, gamma, deltas, norm='l2', long_only=True,
                     recent_frac=0.6, blend_lambda=0.3):
    w_train = estimate_train_weights_for_candidates(R_train, R_val, recent_frac=recent_frac, blend_lambda=blend_lambda)
    return solve_wass_cvar_candidates_weighted(R_train, alpha, gamma, deltas, weights=w_train, norm=norm, long_only=long_only)

def run_experiment_once(d=8, n_train=1000, n_val=1200, n_test=15000,
                        alpha=0.10, gamma=0.06, mean_scale=0.01, vol_abs=0.02,
                        shift=False, mean_shift_abs=-0.004, vol_mult=1.6,
                        phi=0.3, phi_shift=0.45,
                        deltas=None, beta=0.05, B_boot=800, recent_frac=0.6, blend_lambda=0.3):
    if deltas is None:
        deltas = np.linspace(0.001, 0.02, 16)   # include small δ; exclude 0
    mu_p = np.linspace(0.02, 0.002, d) * mean_scale
    Sigma_p = make_cov(d, vol_abs); Sigma_eps = (1 - phi**2) * Sigma_p
    if not shift:
        mu_q, Sigma_q, phi_q = mu_p.copy(), Sigma_p.copy(), phi
    else:
        mu_q = mu_p + mean_shift_abs*np.ones(d); Sigma_q = (vol_mult**2)*Sigma_p; phi_q = phi_shift
    Sigma_eps_q = (1 - phi_q**2) * Sigma_q

    R_train = ar1_sim(mu_p, Sigma_eps, phi, n_train)
    R_val   = ar1_sim(mu_p, Sigma_eps, phi, n_val, x0=R_train[-1])
    R_test  = ar1_sim(mu_q, Sigma_eps_q, phi_q, n_test, x0=R_val[-1])

    X = solve_candidates(R_train, R_val, alpha, gamma, deltas, norm='l2', long_only=True,
                         recent_frac=recent_frac, blend_lambda=blend_lambda)

    sel_new = new_validator_select_delta(X, R_val, alpha, gamma, deltas, beta=beta, B_boot=B_boot, recent_frac=recent_frac)
    sel_old = old_ngs_validator_select_delta(X, R_val, alpha, gamma, deltas, beta=beta)
    sel_iw  = iwcv_select_delta(X, R_val, alpha, gamma, deltas, recent_frac=recent_frac)

    def eval_candidate(x, delta):
        losses = portfolio_loss(R_test, x)
        cvar, _ = empirical_cvar(losses, alpha)
        lhs = cvar + (delta/alpha) * np.linalg.norm(x)
        feas = lhs <= gamma
        cost = (-np.mean(R_train, axis=0)) @ x
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
    keys = recs[0].keys(); out = {}
    for k in keys:
        vals = np.array([float(r[k]) for r in recs]); out[k] = float(np.mean(vals))
    return out

# ---------------- Fair γ ex-ante ----------------
def fix_gamma_ex_ante_eqw_on_P(d: int, n_train: int, alpha: float, mean_scale: float, vol_abs: float, phi: float,
                               margin: float = 0.10, seed: int = 999) -> float:
    mu_p = np.linspace(0.02, 0.002, d) * mean_scale
    Sigma_p = make_cov(d, vol_abs, seed=seed+3); Sigma_eps = (1 - phi**2) * Sigma_p
    R_train_gamma = ar1_sim(mu_p, Sigma_eps, phi, n_train)
    x_ref = np.ones(d)/d
    gamma0 = empirical_cvar(portfolio_loss(R_train_gamma, x_ref), alpha)[0]
    return float((1.0 + margin) * gamma0)

# ---------------- Main ----------------
if __name__ == "__main__":
    common = dict(
        d=8, n_train=1000, n_val=1200, n_test=15000,
        alpha=0.10, mean_scale=0.01, vol_abs=0.02,
        deltas=np.linspace(0.001, 0.02, 16),
        beta=0.05, B_boot=800, recent_frac=0.6, blend_lambda=0.3,
        phi=0.3, phi_shift=0.45
    )
    gamma = fix_gamma_ex_ante_eqw_on_P(
        d=common["d"], n_train=common["n_train"], alpha=common["alpha"],
        mean_scale=common["mean_scale"], vol_abs=common["vol_abs"], phi=common["phi"],
        margin=0.10, seed=31415
    )
    print("Fixed gamma (ex-ante, same for both scenarios):", gamma)

    s1 = dict(shift=False, mean_shift_abs=0.0, vol_mult=1.0, **common)
    res1 = run_scenarios(n_rep=30, gamma=gamma, **s1)
    print("Scenario 1 (no shift):"); print(res1)

    s2 = dict(shift=True, mean_shift_abs=-0.004, vol_mult=1.7, **common)
    res2 = run_scenarios(n_rep=30, gamma=gamma, **s2)
    print("Scenario 2 (shift):"); print(res2)
