import numpy as np
import numpy.matlib
import pandas as pd
from gurobipy import Model, GRB, quicksum, LinExpr
from scipy.linalg import sqrtm
from scipy.stats import norm
from scipy.sparse import eye, tril, block_diag, vstack
from gurobipy import Model, GRB

def CCP_DRO_moment(para, c, b, data, alpha, d, n):
    # para: a sequence of parameter values specifying the size of the uncertainty set
    # c: cost vector in the objective
    # b: right-hand threshold in the constraint
    # alpha: tolerance level in the original chance constraint
    # data: data matrix, each row is an observation

    # Construct Jacobian matrix
    d_aug = d * (d + 1) // 2
    mu_hat = np.mean(data, axis=0)
    grad_A = np.sqrt(alpha / (1 - alpha)) * eye(d).toarray()
    grad_C = eye(d_aug).toarray()
    grad_B = np.zeros((0, d))

    # Construct the lower-left block, grad_B, of the Jacobian matrix
    for col in range(d):
        block_size = d - col
        new_block = np.zeros((block_size, d))
        new_block[:, col:] = -mu_hat[col]
        new_block[:, col] += -mu_hat[col:d]
        grad_B = np.vstack([grad_B, new_block])

    grad = np.vstack([np.hstack([grad_A, np.zeros((d, d_aug))]), np.hstack([grad_B, grad_C])])

    # Compute covariance matrix
    tril_ind = np.tril_indices(d)
    row_ind, col_ind = tril_ind
    data_moment = np.hstack([data, data[:, row_ind] * data[:, col_ind]])
    moment_sigma = np.cov(data_moment, rowvar=False, bias=True)
    V_est = grad @ moment_sigma @ grad.T

    # Solver
    Sigma_hat = np.cov(data, rowvar=False, bias=True)
    tilde_c = -np.sqrt(alpha / (1 - alpha)) * b
    sqrt_cov = sqrtm(V_est)
    svec_operator = np.sqrt(2) * np.ones((d, d)) + (1 - np.sqrt(2)) * np.eye(d)
    svec_multiplier = svec_operator[tril_ind]
    A = np.diag(np.hstack([np.ones(d), svec_multiplier]))

    if isinstance(para, np.floating): para = np.array([para])
    rho = para
    mesh_size = len(rho)
    solution = np.zeros((d, mesh_size))

    for k in range(mesh_size):
        # Create a new model
        model = Model()

        # Add variables
        x_dro = model.addMVar(d, lb=0, ub=1, name="x_dro")
        W = model.addMVar((d, d), name="W")
        vec_q = model.addMVar(d + d_aug, name="vec_q")
        eta = model.addVar(name="eta")

        # Set objective
        model.setObjective(c @ x_dro, GRB.MINIMIZE)
        _sqrt_cov = A @ sqrt_cov

        l = {}
        for i in range(d + d_aug):
            lexpr = LinExpr(0)
            for j in range(d + d_aug):
                lexpr.add(_sqrt_cov[j, i] * vec_q[j])        
            l[i] = lexpr 
        
        aux = []
        for i in range(d + d_aug):
            v = model.addVar(lb=-GRB.INFINITY, name="aux%d"%(i))
            model.addConstr(v == l[i], name="aux_constr%d"%(i))
            aux.append(v)
        
        normaux = model.addVar(name="normaux")
        model.addGenConstrNorm(normaux, aux, 2.0, "normconstr")

        # Add constraints
        model.addConstr(np.sqrt(alpha / (1 - alpha)) * quicksum([mu_hat[i] * x_dro[i] for i in range(d)])
                        + quicksum([Sigma_hat[i, i] * W[i, i] for i in range(d)]) 
                        + rho[k] * normaux + tilde_c + eta / 4 <= 0)
        model.addConstrs(vec_q[i] == x_dro[i] for i in range(d))
        model.addConstrs(vec_q[i + d] == W[tril_ind][i] * svec_multiplier[i] for i in range(d_aug))
        # model.addConstr(vstack([np.hstack([W, x_dro[:, None]]), np.hstack([x_dro, np.array([eta])])]).toarray() >> 0)

        # Add positive semidefinite constraint
        for i in range(d):
            for j in range(i, d):
                if i == j:
                    model.addQConstr(W[i, j] >= eta * (x_dro[i] ** 2))
                else:
                    model.addQConstr(W[i, j] + W[j, i] >= 2 * eta * (x_dro[i] * x_dro[j]))

        # Add constraints for symmetry and PSD
        for i in range(d):
            model.addConstr(W[i, d] == x_dro[i])
            model.addConstr(W[d, i] == x_dro[i])
        model.addConstr(W[d, d] == eta)

        # Optimize model
        model.optimize()

        # Retrieve the solution
        if model.status == GRB.OPTIMAL:
            solution[:, k] = x_dro.X

    return solution if mesh_size > 1 else solution[:, 0]

def CCP_RO_ellipsoid(para, c, b, data, alpha, d, n):
    """
    Parameters:
    para: a sequence of parameter values specifying the size of the uncertainty set
    c: cost vector in the objective
    b: right-hand threshold in the constraint
    data: data matrix, each row is an observation
    """
    mu_hat = np.mean(data, axis=0)
    sigma_hat = np.cov(data, rowvar=False)
    # order1 = binom.ppf(1-beta, n, 1-alpha)
    # if order1 > n - 1:
    #     raise ValueError("too small N1")
    
    # mu_hat_span = np.tile(mu_hat, (n, 1))
    # t_value = np.diag((data - mu_hat_span) @ np.linalg.inv(sigma_hat) @ (data - mu_hat_span).T)
    # t_sorted = np.sort(t_value)
    # t_select = t_sorted[int(order1) + 1]  # the selected parameter
    sigma_hat_rt = np.real(sqrtm(sigma_hat))
    
    if isinstance(para, np.floating): para = np.array([para])
    radius = para
    mesh_size = len(radius)
    solution = np.zeros((d, mesh_size))

    for k in range(mesh_size):
        model = Model()
        model.Params.OutputFlag = 0
        x_ro = model.addVars(d, lb=0, ub=1, vtype=GRB.CONTINUOUS, name="x_ro")
        model.setObjective(quicksum(c[j] * x_ro[j] for j in range(d)), GRB.MINIMIZE)

        l = {}
        for i in range(d):
            lexpr = LinExpr(0)
            for j in range(d):
                lexpr.add(sigma_hat_rt[i, j] * x_ro[j])        
            l[i] = lexpr 
        
        aux = []
        for i in range(d):
            v = model.addVar(lb=-GRB.INFINITY, name="aux%d"%(i))
            model.addConstr(v == l[i], name="aux_constr%d"%(i))
            aux.append(v)
        
        normaux = model.addVar(name="normaux")
        model.addGenConstrNorm(normaux, aux, 2.0, "normconstr")

        model.addConstr(quicksum(mu_hat[j] * x_ro[j] for j in range(d)) + np.sqrt(radius[k]) * normaux - b <= 0)
        model.optimize()
        solution[:, k] = np.array([x_ro[j].X for j in range(d)])
    
    return solution if mesh_size > 1 else solution[:, 0]

def CCP_SO(para, c, b, data, alpha, d, n):
    """
    Parameters
    ----------
    para: a sequence of integers representing the number of scenarios to use
    c: cost vector in the objective
    b: right-hand threshold in the constraint
    data: data matrix, each row is an observation
    d: dimension of x
    """
    if type(para) in [np.int32, np.int64]: para = np.array([para])
    mesh_size = len(para)
    solution = np.zeros((d, mesh_size))

    for k in range(mesh_size):
        model = Model()
        model.Params.OutputFlag = 0
        x_gen = model.addVars(d, lb=0, ub=1, vtype=GRB.CONTINUOUS, name="x_gen")
        model.setObjective(quicksum(c[j] * x_gen[j] for j in range(d)), GRB.MINIMIZE)        
        for i in range(para[k]):
            model.addConstr(quicksum(data[i, j] * x_gen[j] for j in range(d)) <= b)
        model.optimize()
        solution[:, k] = np.array([x_gen[j].X for j in range(d)])

    return solution if mesh_size > 1 else solution[:, 0]

def CCP_SAA(para, c, b, data, alpha, d, n):
    """
    Parameters
    ----------
    para: a sequence of parameter values specifying the tolerance level
    c: cost vector in the objective
    b: right-hand threshold in the constraint
    data: data matrix, each row is an observation
    alpha: tolerance level
    d: dimension of x
    n: rows in data
    """
    M = np.ceil(np.abs(b) + np.max(np.sum(np.abs(data), axis=1)))
    num_constr = para
    if type(num_constr) in [np.int32, np.int64]: num_constr = np.array([num_constr])
    assert((n >= num_constr).all())
    mesh_size = len(num_constr)
    solution = np.zeros((d, mesh_size))
    
    for k in range(mesh_size):
        if k > 0 and num_constr[k] == num_constr[k - 1]:
            solution[:, k] = solution[:, k - 1]
        else:
            model = Model()
            model.Params.OutputFlag = 0
            
            x_SAA = model.addVars(d, lb=0, ub=1, vtype=GRB.CONTINUOUS, name="x_SAA")
            z = model.addVars(n, vtype=GRB.BINARY, name="z")
            model.setObjective(quicksum(c[j] * x_SAA[j] for j in range(d)), GRB.MINIMIZE)
            for i in range(n):
                model.addConstr(quicksum(data[i, j] * x_SAA[j] for j in range(d)) <= b + M * (1 - z[i]))
            model.addConstr(quicksum(z[i] for i in range(n)) >= num_constr[k])
            model.optimize()

            try:
                solution[:, k] = np.array([x_SAA[j].X for j in range(d)])
            except AttributeError:
                continue
    
    return solution if mesh_size > 1 else solution[:, 0]