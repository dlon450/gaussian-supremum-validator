"""Safe convex approximation (SCA) of the chance constraint.

Provides the Bonferroni/Nemirovski-style deterministic surrogate

    mu^T x + sqrt(2 log(1/alpha)) * || Sigma^{1/2} x ||_2 <= b

as a standalone baseline solver. Running the module prints the SCA solution
under the true (half-normal) distribution.

Run with:  python sca.py
"""

import numpy as np
from gurobipy import Model, GRB, quicksum
from numpy.linalg import eigh

def _sqrt_psd(sigma, eps=1e-12):
    """Symmetric PSD square root via eigen-decomposition; numerically stable."""
    # force symmetry
    S = 0.5 * (sigma + sigma.T)
    w, V = eigh(S)
    w_clipped = np.clip(w, 0, None)
    return (V * np.sqrt(w_clipped)) @ V.T

def SCA(c, b, alpha, mu, sigma, lb=0.0, ub=1.0, output_flag=1):
    """
    Safe Convex Approximation:
        minimize   c^T x
        subject to mu^T x + sqrt(2 log(1/alpha)) * || Sigma^{1/2} x ||_2 <= b
                   lb <= x_j <= ub

    Returns:
        x_opt (np.ndarray) if optimal, else NaNs.
    """
    d = len(c)
    c = np.asarray(c, dtype=float).ravel()
    mu = np.asarray(mu, dtype=float).ravel()
    sigma = np.asarray(sigma, dtype=float)

    # factor in the SCA 'radius' term
    phi = np.sqrt(2.0 * np.log(1.0 / alpha))

    # Sigma^{1/2} (stable)
    sigma_rt = _sqrt_psd(sigma)

    m = Model()
    m.Params.OutputFlag = output_flag

    # decision variables (match your RO code: 0 <= x <= 1)
    x = m.addVars(d,
                  lb=(lb if lb is not None else -GRB.INFINITY),
                  ub=(ub if ub is not None else  GRB.INFINITY),
                  vtype=GRB.CONTINUOUS, name="x")

    # objective: minimize c^T x
    m.setObjective(quicksum(float(c[j]) * x[j] for j in range(d)), GRB.MINIMIZE)

    # build linear expressions y = Sigma^{1/2} x
    y = m.addVars(d, lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name="y")
    for i in range(d):
        m.addConstr(y[i] == quicksum(float(sigma_rt[i, j]) * x[j] for j in range(d)))

    # z = || y ||_2
    z = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="z")
    m.addGenConstrNorm(z, [y[i] for i in range(d)], 2.0, "norm_y")

    # SCA constraint: mu^T x + phi * z <= b
    m.addConstr(quicksum(float(mu[j]) * x[j] for j in range(d)) + float(phi) * z <= float(b), name="sca")

    m.optimize()
    if m.Status == GRB.OPTIMAL:
        return np.array([x[j].X for j in range(d)])
    else:
        return np.full(d, np.nan, dtype=float)

def SCA_from_data(c, b, alpha, data, lb=0.0, ub=1.0, output_flag=0):
    """Convenience wrapper using empirical mean/cov from a training set."""
    mu_hat = np.mean(data, axis=0)
    sigma_hat = np.cov(data, rowvar=False)
    return SCA(c, b, alpha, mu_hat, sigma_hat, lb=lb, ub=ub, output_flag=output_flag)

if __name__ == "__main__":
    np.random.seed(2)
    d = 10
    alpha = 0.05
    beta = 0.1
    c = -np.ones(d)
    b = .833 * d

    mu_true = np.sqrt(2/np.pi) * np.ones(d)
    sigma_true = (1 - 2/np.pi) * np.eye(d)
    ksi_v = abs(np.random.normal(size = (1000000,d)))
    x_SCA_true = SCA(c, b, alpha, mu_true, sigma_true, lb=0.0, ub=1.0)
    obj_SCA_true = float(np.dot(c, x_SCA_true))
    P_hat_SCA_true = float(np.mean((ksi_v @ x_SCA_true) <= b))
    feasible_SCA_true = (P_hat_SCA_true >= 1 - alpha)
    print(f"SCA (true dist): obj {obj_SCA_true:.4f}, P_hat {P_hat_SCA_true:.4f}, feasible {feasible_SCA_true}")