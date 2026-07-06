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
    num_constr = np.ceil(n * (1 - alpha + np.array(para))).astype(int)  # transform the tolerance level into the number of satisfied constraints
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

def algo1(delta, c, b, alpha,n,d,ksi,beta,K):
    m = len(delta)
    l = int(n / K)
    V = np.zeros(shape = (K,m))
    C = np.zeros(shape = (1,m))
    ksi_shuffle = np.random.permutation(ksi)
    for i in range(m):
        cost = 0
        for j in range(K):
            ksi_validation = ksi_shuffle[j * l:(j + 1) * l]
            ksi_train = np.vstack((ksi_shuffle[:j * l],ksi_shuffle[(j + 1) * l:]))
            # x_ = solution(delta[i], c, b, alpha, n - l, d,ksi_train)
            x_ = CCP_RO_ellipsoid(delta[i], c, b, ksi_train, alpha, d, n - l)
            P = np.dot(ksi_validation, x_ ) - b
            cost += np.dot(c,x_)
            p = len(P[P < 0]) / l
            V[j,i] = p
        C[0,i] = cost / K
    qualified = np.zeros(shape = (1,m))
    avg_ind = np.zeros(shape = (1,m))
    for i in range(m):
        avg_ind[0,i] = len(V[:,i][V[:,i] >= 1 - alpha]) / K
        if avg_ind[0,i] >= 1 - beta:
            qualified[0,i] = 1
    if sum(sum(qualified)) == 0:
        min_ind = np.where(avg_ind == np.max(avg_ind))[1][0]
    else:
        min_val = np.Infinity
        for i in range(m):
            if qualified[0,i] == 1 and C[0,i] < min_val:
                min_val = C[0,i]
                min_ind = i
    delta_star = delta[min_ind]
    # x = solution(delta_star, c, b, alpha, n, d, ksi)
    x = CCP_RO_ellipsoid(delta_star, c, b, ksi, alpha, d, n)
    return x 
            
def algo2(delta, c, b, alpha,n,d,ksi,beta,B):
    R = np.random.randint(0,n,size = (B,n))
    setup = np.zeros(shape = (B,n))
    setup_c = np.zeros(shape = (B,n))
    for i in range(B):
        for j in range(n):
            if j in R[i]:
                setup[i,j] = 1
            else:
                setup_c[i,j] = 1
    m = len(delta)
    denom = np.sum(setup_c,axis = 1)
    V = np.zeros(shape = (B,m))
    C = np.zeros(shape = (1,m))
    for i in range(m):
        cost = 0
        for j in range(B):
            ksi_validation = ksi[setup_c[j] == 1]
            ksi_train = ksi[R[j,:],:]
            # x_ = solution(delta[i], c, b, alpha, n, d,ksi_train)
            x_ = CCP_RO_ellipsoid(delta[i], c, b, ksi_train, alpha, d, n)
            P = np.dot(ksi_validation, x_ ) - b
            cost += np.dot(c,x_)
            p = len(P[P < 0]) / denom[j]
            V[j,i] = p
        C[0,i] = cost / B
    qualified = np.zeros(shape = (1,m))
    avg_ind = np.zeros(shape = (1,m)) 
    for i in range(m):
        avg_ind[0,i] = len(V[:,i][V[:,i] >= 1 - alpha]) / B
        if avg_ind[0,i] >= 1 - beta:
            qualified[0,i] = 1
    if sum(sum(qualified)) == 0:
        min_ind = np.where(avg_ind == np.max(avg_ind))[1][0]
    else:
        min_val = np.Infinity
        for i in range(m):
            if qualified[0,i] == 1 and C[0,i] < min_val:
                min_val = C[0,i]
                min_ind = i
    delta_star = delta[min_ind]
    # x = solution(delta_star, c, b, alpha, n, d,ksi)
    x = CCP_RO_ellipsoid(delta_star, c, b, ksi, alpha, d, n)
    return x

def algo3(delta, c, b, alpha,n,d,ksi,beta,n1,n2,K):
    ksi_shuffle = np.random.permutation(ksi)
    ksi_train = ksi_shuffle[:n1]
    ksi_v = ksi_shuffle[n1:]
    m = len(delta)
    l = int(n2 / K)
    V = np.zeros(shape = (K,m))
    C = np.zeros(shape = (1,m))
    xx = np.zeros(shape = (d,m))
    for i in range(m):
        x_ = CCP_RO_ellipsoid(delta[i], c, b, ksi_train, alpha, d, n1)
        xx[:,i] = x_
        for j in range(K):
            ksi_validation = ksi_v[j * l:(j + 1) * l]
            P = np.dot(ksi_validation, x_ ) - b
            p = len(P[P < 0]) / l
            V[j,i] = p
        C[0,i] = np.dot(c,x_)
    qualified = np.zeros(shape = (1,m))
    avg_ind = np.zeros(shape = (1,m)) 
    for i in range(m):
        avg_ind[0,i] = len(V[:,i][V[:,i] >= 1 - alpha]) / K
        if avg_ind[0,i] >= 1 - beta:
            qualified[0,i] = 1
    if sum(sum(qualified)) == 0:
        min_ind = np.where(avg_ind == np.max(avg_ind))[1][0]
    else:
        min_val = np.Infinity
        for i in range(m):
            if qualified[0,i] == 1 and C[0,i] < min_val:
                min_val = C[0,i]
                min_ind = i
    x = xx[:,min_ind]
    return x

def GS(delta, c, b, alpha, n, d, ksi, beta, n1, n2):
    ksi_shuffle = np.random.permutation(ksi)
    ksi_train = ksi_shuffle[:n1]
    ksi_v = ksi_shuffle[n1:]
    xx = CCP_RO_ellipsoid(delta, c, b, ksi_train, alpha, d, n1)
    C = np.dot(c,xx)
    res = np.zeros(shape=(d,4))

    indifunc = (np.dot(ksi_v,xx) <= b).astype(int)
    Sigma_hat = np.cov(indifunc.T)
    sigma_hat = np.sqrt(np.diag(Sigma_hat))
    ind_nonzero = sigma_hat > 0
    P_hat = np.mean(indifunc,axis = 0)

    # Univariate Gaussian
    feasible_UG = P_hat >= 1 - alpha + norm.ppf(1 - beta) * sigma_hat / np.sqrt(n2)

    if np.sum(feasible_UG) > 0:
        xx_feasible = xx[:, feasible_UG]
        index1 = np.argmin(C[feasible_UG])
        X_UG = xx_feasible[:, index1]
    else:
        ind = np.argmax(P_hat - norm.ppf(1 - beta) * sigma_hat / np.sqrt(n2))
        X_UG = xx[:, ind]
    
    # Normalized and Unnormalized Gaussian Supremum 
    if np.any(ind_nonzero):
        Sigma_hat = Sigma_hat[ind_nonzero][:,ind_nonzero]
        sim_num = 2000
        l = len(Sigma_hat)
        Z = np.random.multivariate_normal(np.zeros(l),Sigma_hat,sim_num)
        tmat = np.matlib.repmat(sigma_hat[ind_nonzero], sim_num,1)
        q = np.amax(Z / tmat, axis = 1)
        qu = np.amax(Z, axis = 1)
        q_hat = np.percentile(q,(1 - beta) * 100)
        qu_hat = np.percentile(qu,(1 - beta) * 100)
    else:
        q_hat, qu_hat = 0, 0
    
    feasible_NGS = (P_hat >= 1 - alpha + q_hat * sigma_hat / np.sqrt(n2))
    if np.any(feasible_NGS):
        xx_feasible = xx[:,feasible_NGS]
        C_feasible = C[feasible_NGS]
        ind = np.argmin(C_feasible)
        X_NGS = xx_feasible[:,ind]
    else:
        ind = np.argmax(P_hat - q_hat * sigma_hat/ np.sqrt(n2))
        X_NGS = xx[:,ind]
        
    feasible_UNGS = (P_hat >= 1 - alpha + qu_hat / np.sqrt(n2))
    if np.any(feasible_UNGS):
        xx_feasible = xx[:,feasible_UNGS]
        C_feasible = C[feasible_UNGS]
        ind = np.argmin(C_feasible)
        X_UNGS = xx_feasible[:,ind]
    else:
        ind = np.argmax(P_hat - qu_hat / np.sqrt(n2))
        X_UNGS = xx[:,ind]

    # Baseline 
    feasible_NV = (P_hat >= 1 - alpha)
      
    if np.any(feasible_NV):
        xx_feasible = xx[:,feasible_NV]
        C_feasible = C[feasible_NV]
        ind = np.argmin(C_feasible)
        X_NV = xx_feasible[:,ind]
    else:
        ind = np.argmax(P_hat)
        X_NV = xx[:,ind]
        
    res[:,0],res[:,1],res[:,2],res[:,3] = X_NGS, X_UNGS, X_NV, X_UG
    return res

def perform_stats(c,N,d,n,alpha,beta,K,B,n1,n2,delta,b,dim=2):
    x = np.zeros(shape=(d,dim))
    collect1 = np.zeros(shape=(N,dim))
    collect2 = np.zeros(shape=(N,dim))
    for i in range(N):
        ksi = abs(np.random.normal(size=(n,d)))
        if dim == 2:
            x[:,0] = algo1(delta, c, b, alpha, n, d, ksi, beta, K)
            x[:,1] = algo2(delta, c, b, alpha, n, d, ksi, beta, B)
        elif dim == 5:
            x[:,0] = algo3(delta, c, b, alpha, n, d, ksi, beta, n1, n2, K)
            x[:,1:] = GS(delta, c, b, alpha, n, d, ksi, beta, n1, n2)
        elif dim == 4:
            x[:,0:] = GS(delta, c, b, alpha, n, d, ksi, beta, n1, n2)
        else:
            raise Exception('Incorrect dimension d')
        metric1 = np.dot(c,x)
        ksi_ = abs(np.random.normal(size = (1000000,d)))
        metric2 = (np.sum((np.dot(ksi_,x) <= b).astype(int),0) / 1000000 >= 1 - alpha).astype(int)
        collect1[i] = metric1
        collect2[i] = metric2
    mean_obj, feas_level = np.mean(collect1,0), np.mean(collect2,0)
    return mean_obj, feas_level, collect1, collect2

np.random.seed(2)
d = 10
n = 500
n1, n2 = n // 2, n // 2
data = np.abs(np.random.normal(size=(n, d)))
delta1 = [3 / 2**i for i in range(10)]
delta2 = [delta1[-1] / 1.2**i for i in range(1,16)]
delta = np.array(delta1 + delta2)

alpha = 0.05
beta = 0.1
c = -np.ones(d)
b = .833 * d
N = 100
K = B = 10

dim = 5
column_names = ['Cross Validation','Bootstrapping','Sectioning','NGS','UNGS','NV','UG']
column_names_ = column_names[2:] if dim == 5 else column_names[3:] if dim == 4 else column_names[:2] 
output = 'experiment_results_ro_ellipsoid.csv'

############ RO ###########
delta = 45 * np.arange(1, 26) / 25

for n in [100, 200, 300, 400, 500]:
    print('Data Size:', n)
    for i, n1 in enumerate(range(n // 10, n, n // 10)):
        n2 = n - n1
    # for K in [10]:
    #     B = K

        # delta_SO = np.arange(n1 // 5, n1 + 1, n1 // 5)
        result = perform_stats(c,N,d,n,alpha,beta,K,B,n1,n2,delta=delta,b=b,dim=dim)
        print(result[0], result[1])

        try:
            df1 = pd.read_csv(f'objectives_{n}.csv')
            df2 = pd.read_csv(f'feas_levels_{n}.csv')
        except (pd.errors.EmptyDataError, FileNotFoundError) as e:
            df1, df2 = pd.DataFrame([]), pd.DataFrame([])

        dfr = pd.read_csv(output)
        dfr_ = {'Objective': result[0], 
                'Feasibility': result[1], 
                'Algorithm': column_names_, 
                'Training Percentage': [10 * (i + 1) if dim >= 4 else 100] * dim, 
                'Data Size': [n] * dim, 
                'Folds': [0 if dim >= 4 else K] * dim}
        pd.concat([df1, pd.DataFrame(result[2], columns=column_names_)], axis=1).to_csv(f'objectives_{n}.csv', index=False)
        pd.concat([df2, pd.DataFrame(result[3], columns=column_names_)], axis=1).to_csv(f'feas_levels_{n}.csv', index=False)
        pd.concat([dfr, pd.DataFrame(dfr_)]).to_csv(output, index=False)