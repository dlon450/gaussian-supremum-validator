"""Base optimization formulations for the data-driven chance-constrained program.

Every solver approximately enforces the chance constraint

    P( xi^T x <= b ) >= 1 - alpha

for a cost vector ``c`` (the objective ``min c^T x``) using a different
uncertainty model, and each is parameterised by a scalar (or 1-D array)
``para`` controlling how conservative the approximation is.

To make the solvers interchangeable inside the validation routines
(:mod:`gsv.validators`), they all share the signature

    solve(para, c, b, data, alpha, d, n) -> np.ndarray

where

    para   scalar or 1-D array of uncertainty-set / conservativeness levels
    c      cost vector in the objective                      (length d)
    b      right-hand threshold in the constraint            (scalar)
    data   observation matrix, one sample per row            (n x d)
    alpha  tolerance level of the original chance constraint (scalar)
    d      dimension of the decision variable x
    n      number of observations (rows of ``data`` actually used)

The return value is a length-``d`` solution vector when ``para`` is scalar,
or a ``d x len(para)`` matrix (one column per parameter value) otherwise.

The heavyweight optional dependencies (``cvxpy`` for the moment-DRO solver
and ``rsome`` for the portfolio solver) are imported lazily so that the
Gurobi-only solvers can be used without them installed.
"""

import os
import numpy as np
from scipy.linalg import sqrtm
from scipy.sparse import eye

_MILP_TL = float(os.environ.get("GSV_MILP_TIMELIMIT", "5"))     # MILP wall-clock cap (s)
_MILP_GAP = float(os.environ.get("GSV_MILP_GAP", "1e-2"))       # MILP relative gap
try:                                            # Gurobi is optional: the OSS backend
    from gurobipy import Model, GRB, quicksum, LinExpr   # (gsv.solvers_oss) covers the
except Exception:                               # convex formulations without it.
    Model = GRB = quicksum = LinExpr = None


def CCP_DRO_moment(para, c, b, data, alpha, d, n):
    """Moment-based distributionally robust chance constraint.

    Uses the first two empirical moments of ``data`` together with a
    (semidefinite) uncertainty set of radius ``para`` estimated from the
    sampling covariance of those moments. Solved with cvxpy.
    """
    import cvxpy as cp

    # Construct the Jacobian of the moment map (mean and lower-triangular
    # second moments) so its sampling covariance sizes the uncertainty set.
    d_aug = d * (d + 1) // 2
    mu_hat = np.mean(data, axis=0)
    grad_A = np.sqrt(alpha / (1 - alpha)) * eye(d).toarray()
    grad_C = eye(d_aug).toarray()
    grad_B = np.zeros((0, d))

    # Lower-left block, grad_B, of the Jacobian matrix.
    for col in range(d):
        block_size = d - col
        new_block = np.zeros((block_size, d))
        new_block[:, col:] = -mu_hat[col]
        new_block[:, col] += -mu_hat[col:d]
        grad_B = np.vstack([grad_B, new_block])

    grad = np.vstack([np.hstack([grad_A, np.zeros((d, d_aug))]),
                      np.hstack([grad_B, grad_C])])

    # Sampling covariance of the moment estimates.
    tril_ind = np.tril_indices(d)
    row_ind, col_ind = tril_ind
    data_moment = np.hstack([data, data[:, row_ind] * data[:, col_ind]])
    moment_sigma = np.cov(data_moment, rowvar=False, bias=True)
    V_est = grad @ moment_sigma @ grad.T

    Sigma_hat = np.cov(data, rowvar=False, bias=True)
    tilde_c = -np.sqrt(alpha / (1 - alpha)) * b
    sqrt_cov = sqrtm(V_est)
    svec_operator = np.sqrt(2) * np.ones((d, d)) + (1 - np.sqrt(2)) * np.eye(d)
    svec_multiplier = svec_operator[tril_ind]
    A = np.diag(np.hstack([np.ones(d), svec_multiplier]))

    if isinstance(para, np.floating):
        para = np.array([para])
    rho = para
    mesh_size = len(rho)
    solution = np.zeros((d, mesh_size))

    for k in range(mesh_size):
        # Decision variables.
        x_dro = cp.Variable(d)
        W = cp.Variable((d, d), symmetric=True)
        vec_q = cp.Variable(d + d_aug)
        psd = cp.Variable((d + 1, d + 1), symmetric=True)
        eta = cp.Variable()

        objective = cp.Minimize(c @ x_dro)
        constraints = [
            np.sqrt(alpha / (1 - alpha)) * mu_hat @ x_dro
            + cp.trace(Sigma_hat @ W)
            + rho[k] * cp.norm((A @ sqrt_cov).T @ vec_q)
            + tilde_c + eta / 4 <= 0,
            vec_q == cp.hstack([x_dro, cp.multiply(W[tril_ind], svec_multiplier)]),
            psd[:d, :d] == W,
            psd[d, :d] == x_dro.T,
            psd[:d, d] == x_dro,
            psd[d, d] == eta,
            psd >> 0,
            0 <= x_dro,
            x_dro <= 1,
        ]

        prob = cp.Problem(objective, constraints)
        prob.solve()
        solution[:, k] = x_dro.value if x_dro.value is not None else np.nan

    return solution if mesh_size > 1 else solution[:, 0]


def CCP_RO_ellipsoid(para, c, b, data, alpha, d, n):
    """Robust ellipsoidal approximation of the chance constraint.

    Enforces ``mu_hat^T x + sqrt(para) * || Sigma_hat^{1/2} x ||_2 <= b`` over an
    ellipsoidal uncertainty set of squared radius ``para``.
    """
    mu_hat = np.mean(data, axis=0)
    sigma_hat = np.cov(data, rowvar=False)
    sigma_hat_rt = np.real(sqrtm(sigma_hat))

    if isinstance(para, np.floating):
        para = np.array([para])
    radius = para
    mesh_size = len(radius)
    solution = np.zeros((d, mesh_size))

    for k in range(mesh_size):
        model = Model()
        model.Params.OutputFlag = 0
        model.Params.Threads = 1        # one solver thread per worker (avoid oversubscription)
        x_ro = model.addVars(d, lb=0, ub=1, vtype=GRB.CONTINUOUS, name="x_ro")
        model.setObjective(quicksum(c[j] * x_ro[j] for j in range(d)), GRB.MINIMIZE)

        # Auxiliary variables holding Sigma_hat^{1/2} x, used by the norm constraint.
        aux = []
        for i in range(d):
            lexpr = LinExpr(0)
            for j in range(d):
                lexpr.add(sigma_hat_rt[i, j] * x_ro[j])
            v = model.addVar(lb=-GRB.INFINITY, name="aux%d" % i)
            model.addConstr(v == lexpr, name="aux_constr%d" % i)
            aux.append(v)

        normaux = model.addVar(name="normaux")
        model.addGenConstrNorm(normaux, aux, 2.0, "normconstr")

        model.addConstr(quicksum(mu_hat[j] * x_ro[j] for j in range(d))
                        + np.sqrt(radius[k]) * normaux - b <= 0)
        model.optimize()
        solution[:, k] = np.array([x_ro[j].X for j in range(d)])

    return solution if mesh_size > 1 else solution[:, 0]


def CCP_SO(para, c, b, data, alpha, d, n):
    """Scenario optimization: enforce the constraint on the first ``para`` samples.

    ``para`` is an integer (or 1-D integer array) giving the number of scenarios.
    """
    if type(para) in [np.int32, np.int64]:
        para = np.array([para])
    mesh_size = len(para)
    solution = np.zeros((d, mesh_size))

    for k in range(mesh_size):
        model = Model()
        model.Params.OutputFlag = 0
        model.Params.Threads = 1        # one solver thread per worker (avoid oversubscription)
        x_gen = model.addVars(d, lb=0, ub=1, vtype=GRB.CONTINUOUS, name="x_gen")
        model.setObjective(quicksum(c[j] * x_gen[j] for j in range(d)), GRB.MINIMIZE)
        for i in range(para[k]):
            model.addConstr(quicksum(data[i, j] * x_gen[j] for j in range(d)) <= b)
        model.optimize()
        solution[:, k] = np.array([x_gen[j].X for j in range(d)])

    return solution if mesh_size > 1 else solution[:, 0]


def CCP_SAA(para, c, b, data, alpha, d, n):
    """Sample-average approximation with a big-M chance-constraint relaxation.

    ``para`` is the required number of satisfied constraints (an integer, or a
    1-D integer array). At least ``para`` of the ``n`` sampled constraints must
    hold; the rest are relaxed by the big-M term.
    """
    M = np.ceil(np.abs(b) + np.max(np.sum(np.abs(data), axis=1)))
    num_constr = para
    if type(num_constr) in [np.int32, np.int64]:
        num_constr = np.array([num_constr])
    assert (n >= num_constr).all()
    mesh_size = len(num_constr)
    solution = np.zeros((d, mesh_size))

    for k in range(mesh_size):
        if k > 0 and num_constr[k] == num_constr[k - 1]:
            solution[:, k] = solution[:, k - 1]
        else:
            model = Model()
            model.Params.OutputFlag = 0
            model.Params.Threads = 1    # single-threaded: avoid oversubscription
            model.Params.Seed = 0       # deterministic MILP search (reproducible across runs)
            model.Params.MIPGap = _MILP_GAP
            model.Params.TimeLimit = _MILP_TL

            x_SAA = model.addVars(d, lb=0, ub=1, vtype=GRB.CONTINUOUS, name="x_SAA")
            z = model.addVars(n, vtype=GRB.BINARY, name="z")
            model.setObjective(quicksum(c[j] * x_SAA[j] for j in range(d)), GRB.MINIMIZE)
            for i in range(n):
                model.addConstr(quicksum(data[i, j] * x_SAA[j] for j in range(d))
                                <= b + M * (1 - z[i]))
            model.addConstr(quicksum(z[i] for i in range(n)) >= num_constr[k])
            model.optimize()

            try:
                solution[:, k] = np.array([x_SAA[j].X for j in range(d)])
            except AttributeError:
                # solver returned no incumbent: mark this candidate invalid (NaN) so it
                # is masked out downstream, rather than silently left as x=0 (which is
                # trivially chance-feasible and would inflate coverage).
                solution[:, k] = np.nan

    return solution if mesh_size > 1 else solution[:, 0]


def CCP_DRO_wasserstein(para, c, b, data, alpha, d, n):
    """Wasserstein-ball distributionally robust chance constraint (big-M MILP).

    ``para`` is the Wasserstein radius (a scalar, or 1-D array of radii). This
    consolidates the legacy ``solution`` (scalar radius) and ``solution_``
    (vector of radii) routines into a single interface; the models built and
    solved are identical to the originals.
    """
    ksi = data
    M = abs(b) + max(np.sum(ksi, 1))

    scalar = isinstance(para, (float, np.floating))
    para = np.atleast_1d(para)
    mesh_size = len(para)
    solution = np.zeros((d, mesh_size))

    for k in range(mesh_size):
        model = Model("NewModel")
        model.Params.OutputFlag = 0
        model.Params.Threads = 1        # single-threaded: avoid oversubscription
        model.Params.Seed = 0           # deterministic MILP search (reproducible across runs)
        model.Params.MIPGap = _MILP_GAP
        model.Params.TimeLimit = _MILP_TL

        y = {}
        z = {}
        s = {}
        x = {}
        for i in range(n):
            y[i] = model.addVar(vtype=GRB.BINARY)
            z[i] = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)
            s[i] = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)
        for i in range(d):
            x[i] = model.addVar(lb=0, ub=1, vtype=GRB.CONTINUOUS)
        v = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)
        r = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)

        model.setObjective(quicksum(c[i] * x[i] for i in range(d)), GRB.MINIMIZE)
        model.addConstrs(z[i] + s[i] >= r for i in range(n))
        model.addConstr(para[k] * v + quicksum(z[i] for i in range(n)) / n <= alpha * r)
        for i in range(n):
            model.addConstr(b - quicksum(ksi[i, j] * x[j] for j in range(d))
                            + M * (1 - y[i]) >= s[i])
        model.addConstrs(M * y[i] >= s[i] for i in range(n))
        model.addConstrs(v >= x[i] for i in range(d))

        model.optimize()
        # under a TimeLimit an incumbent is virtually always found quickly; if not,
        # fall back to the trivially feasible conservative solution x=0.
        if model.SolCount > 0:
            solution[:, k] = np.array([x[j].x for j in range(d)])
        else:
            solution[:, k] = np.nan

    return solution[:, 0] if scalar else solution


def mean_risk_portfolio_dro(para, zhat, n=30, d=10, alpha=0.2):
    """Mean-risk (CVaR) portfolio under a type-1 Wasserstein ambiguity set.

    Distinct from the chance-constrained solvers above: it returns portfolio
    weights that sum to one and optimises a mean-CVaR objective, so it keeps its
    native signature and is adapted to the ``solve(...)`` interface by the caller.
    Solved with rsome.
    """
    from rsome import dro, E
    import rsome as rso

    rho = 10                                  # risk-aversion coefficient
    a1, b1 = -1, rho                          # coefficients of the piecewise expression
    a2, b2 = -1 - rho / alpha, rho - rho / alpha

    if isinstance(para, np.floating):
        para = np.array([para])
    radius = para
    mesh_size = len(radius)
    solution = np.zeros((d, mesh_size))

    for k in range(mesh_size):
        model = dro.Model(n)
        x = model.dvar(d)
        tau = model.dvar()

        z = model.rvar(d)
        u = model.rvar()
        fset = model.ambiguity()
        for s in range(n):
            fset[s].suppset(rso.norm(z - zhat[s], 1) <= u, z >= -1)
        fset.exptset(E(u) <= radius[k])
        pr = model.p
        fset.probset(pr == 1 / n)

        r = z @ x
        model.minsup(E(rso.maxof(a1 * r + b1 * tau, a2 * r + b2 * tau)), fset)
        model.st(x.sum() == 1, x >= 0)
        model.solve()
        solution[:, k] = x.get()

    return solution if mesh_size > 1 else solution[:, 0]
