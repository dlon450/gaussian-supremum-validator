from rsome import dro
from rsome import E
import rsome as rso
import numpy as np

def mean_risk_portfolio_dro(para, zhat, n=30, d=10, alpha=0.2):

    rho = 10                                  # risk-aversion coefficient
    a1, b1 = -1, rho                          # coefficients of the piecewise expression
    a2, b2 = -1 - rho/alpha, rho - rho/alpha  # coefficients of the piecewise expression

    if isinstance(para, np.floating): para = np.array([para])
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
        fset.probset(pr == 1/n)

        r = z @ x
        model.minsup(E(rso.maxof(a1*r + b1*tau, 
                                a2*r + b2*tau)), fset)
        model.st(x.sum() == 1, x >= 0)
        model.solve()
        solution[:, k] = x.get()

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
            x_ = mean_risk_portfolio_dro(para=delta[i], zhat=ksi_train, n=n-l, d=d, alpha=alpha)
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
    x = mean_risk_portfolio_dro(para=delta_star, zhat=ksi, n=n, d=d, alpha=alpha)
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
            x_ = mean_risk_portfolio_dro(para=delta[i], zhat=ksi_train, n=n, d=d, alpha=alpha)
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
    x = mean_risk_portfolio_dro(para=delta_star, zhat=ksi, n=n, d=d, alpha=alpha)
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
        x_ = mean_risk_portfolio_dro(para=delta[i], zhat=ksi_train, alpha=alpha, n=n1, d=d)
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
    xx = mean_risk_portfolio_dro(para=delta, zhat=ksi_train, n=n1, d=d, alpha=alpha)
    C = np.dot(c,xx)
    res = np.zeros(shape = (d,3))
    indifunc = (np.dot(ksi_v,xx) <= b).astype(int)
    Sigma_hat = np.cov(indifunc.T)
    sigma_hat = np.sqrt(np.diag(Sigma_hat))
    ind_nonzero = sigma_hat > 0
    
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
    
    P_hat = np.mean(indifunc,axis = 0)
    
    feasible_NGS = (P_hat >= 1 - alpha + q_hat * sigma_hat / np.sqrt(n2))
    if np.any(feasible_NGS):
        xx_feasible = xx[:,feasible_NGS]
        C_feasible = C[feasible_NGS]
        ind = np.argmin(C_feasible)
        X_NGS = xx_feasible[:,ind]
    else:
        ind = np.argmax(P_hat - q_hat * sigma_hat/ np.sqrt(n2))
        X_NGS = xx[:,ind]
        
    feasible_NV = (P_hat >= 1 - alpha)
      
    if np.any(feasible_NV):
        xx_feasible = xx[:,feasible_NV]
        C_feasible = C[feasible_NV]
        ind = np.argmin(C_feasible)
        X_NV = xx_feasible[:,ind]
    else:
        ind = np.argmax(P_hat)
        X_NV = xx[:,ind]
        
    feasible_UNGS = (P_hat >= 1 - alpha + qu_hat / np.sqrt(n2))
    if np.any(feasible_UNGS):
        xx_feasible = xx[:,feasible_UNGS]
        C_feasible = C[feasible_UNGS]
        ind = np.argmin(C_feasible)
        X_UNGS = xx_feasible[:,ind]
    else:
        ind = np.argmax(P_hat - qu_hat / np.sqrt(n2))
        X_UNGS = xx[:,ind]
        
    res[:,0],res[:,1],res[:,2] = X_NGS, X_UNGS, X_NV
    return res

def perform_stats(c,N,d,n,alpha,beta,K,B,n1,n2,delta,b,dim=2):
    x = np.zeros(shape=(d,dim))
    collect1 = np.zeros(shape=(N,dim))
    collect2 = np.zeros(shape=(N,dim))
    for i in range(N):
        ksi = abs(np.random.normal(size = (n,d)))
        x[:,0] = algo1(delta, c, b, alpha, n, d, ksi, beta, K)
        x[:,1] = algo2(delta, c, b, alpha, n, d, ksi, beta, B)
        x[:,3] = algo3(delta, c, b, alpha, n, d, ksi, beta, n1, n2, K)
        x[:,4:] = GS(delta, c, b, alpha, n, d, ksi, beta, n1, n2)
        metric1 = np.dot(c,x)
        ksi_ = abs(np.random.normal(size = (1000000,d)))
        metric2 = (np.sum((np.dot(ksi_,x) <= b).astype(int),0) / 1000000 >= 1 - alpha).astype(int)
        collect1[i] = metric1
        collect2[i] = metric2
    mean_obj, feas_level = np.mean(collect1,0), np.mean(collect2,0)
    return mean_obj, feas_level, collect1, collect2

n, d = 100, 10
alpha = 0.2
beta = 0.1
i = np.arange(1, d+1)
np.random.seed(1)
phi = 0.02 * np.random.normal(size=(n, d))
zeta = 0.03*i + 0.025*i*np.random.normal(size=(n, d))
zhat = np.maximum(phi + zeta, -1) # historical data of yearly return
delta = set([b * 10 ** c for b in range(10) for c in [-3, -2, -1]])

dim = 4
column_names = ['Cross Validation','Bootstrapping','Sectioning','NGS','UNGS','NV']
column_names_ = column_names[2:] if dim == 4 else column_names[:2]
output = 'experiment_results_so.csv'

# for n in [1000, 900, 800, 700, 600]:
#     print('Data Size:', n)
#     for i, n1 in enumerate(range(n // 10, n, n // 10)):
n1 = 50
n2 = n - n1

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
        'Training Percentage': [10 * (i + 1) if dim == 4 else 100] * dim, 
        'Data Size': [n] * dim, 
        'Folds': [0 if dim == 4 else K] * dim}
pd.concat([df1, pd.DataFrame(result[2], columns=column_names_)], axis=1).to_csv(f'objectives_{n}.csv', index=False)
pd.concat([df2, pd.DataFrame(result[3], columns=column_names_)], axis=1).to_csv(f'feas_levels_{n}.csv', index=False)
pd.concat([dfr, pd.DataFrame(dfr_)]).to_csv(output, index=False)