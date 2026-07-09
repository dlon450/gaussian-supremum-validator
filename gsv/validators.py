"""Data-driven validation / parameter-selection procedures.

Given a base solver ``solve`` (any of the formulations in :mod:`gsv.formulations`,
all sharing the ``solve(para, c, b, data, alpha, d, n)`` interface) and a grid
of candidate parameter values ``delta``, these routines pick a parameter and
return the corresponding solution. They implement the competing methods
compared in the paper:

    cross_validation    K-fold cross validation                 (algo1)
    bootstrapping       bootstrap resampling                    (algo2)
    sectioning          single train/validate split (sections)  (algo3)
    gaussian_supremum   the Gaussian-supremum validator + variants (GS)

``gaussian_supremum`` returns four candidate solutions, one per selection rule:

    NGS   Normalized Gaussian Supremum   (this paper)
    UNGS  Unnormalized Gaussian Supremum (this paper)
    NV    Naive validation, no margin    (baseline)
    UG    Univariate-Gaussian margin     (baseline)

``run_experiment`` (formerly ``perform_stats``) is the Monte-Carlo driver that
repeats an experiment ``N`` times and reports mean objective and empirical
feasibility for each method.

All random draws happen through NumPy's global RNG, so seed it once in the
calling script (e.g. ``np.random.seed(2)``) for reproducible experiments.
"""

import numpy as np
from scipy.stats import norm


def _select_index(V, C, alpha, beta, folds):
    """Pick the cheapest parameter whose validated feasibility clears the margin.

    A parameter column is "qualified" if at least a ``1 - beta`` fraction of the
    ``folds`` validation splits report feasibility ``>= 1 - alpha``. Among
    qualified columns the cheapest (lowest ``C``) is chosen; if none qualify, the
    column with the highest qualification rate is used instead.
    """
    m = V.shape[1]
    qualified = np.zeros((1, m))
    avg_ind = np.zeros((1, m))
    for i in range(m):
        avg_ind[0, i] = len(V[:, i][V[:, i] >= 1 - alpha]) / folds
        if avg_ind[0, i] >= 1 - beta:
            qualified[0, i] = 1
    if sum(sum(qualified)) == 0:
        min_ind = np.where(avg_ind == np.max(avg_ind))[1][0]
    else:
        min_val = np.inf
        for i in range(m):
            if qualified[0, i] == 1 and C[0, i] < min_val:
                min_val = C[0, i]
                min_ind = i
    return min_ind


def cross_validation(solve, delta, c, b, alpha, n, d, ksi, beta, K, param_grid=None):
    """K-fold cross validation over the parameter grid ``delta``.

    ``param_grid`` (optional) is a callable ``n_eff -> delta`` used by the
    scenario methods (SO/SAA) whose parameter is a discrete constraint count
    that must be regenerated for the ``n - l`` training size and then rescaled
    back to the full sample size ``n``.
    """
    l = int(n / K)
    if param_grid is not None:
        delta = param_grid(n - l)
    m = len(delta)
    V = np.zeros(shape=(K, m))
    C = np.zeros(shape=(1, m))
    ksi_shuffle = np.random.permutation(ksi)
    for i in range(m):
        cost = 0
        for j in range(K):
            ksi_validation = ksi_shuffle[j * l:(j + 1) * l]
            ksi_train = np.vstack((ksi_shuffle[:j * l], ksi_shuffle[(j + 1) * l:]))
            x_ = solve(delta[i], c, b, ksi_train, alpha, d, n - l)
            P = np.dot(ksi_validation, x_) - b
            cost += np.dot(c, x_)
            p = len(P[P < 0]) / l
            V[j, i] = p
        C[0, i] = cost / K
    min_ind = _select_index(V, C, alpha, beta, K)
    if param_grid is not None:
        delta_star = (delta[min_ind] / (n - l) * n).astype(int)
    else:
        delta_star = delta[min_ind]
    return solve(delta_star, c, b, ksi, alpha, d, n)


def bootstrapping(solve, delta, c, b, alpha, n, d, ksi, beta, B, param_grid=None):
    """Bootstrap validation: ``B`` resamples, each validated on its out-of-bag set.

    ``param_grid`` (optional) is a callable ``n -> delta`` used by the scenario
    methods (SO/SAA), whose discrete parameter is defined on the full sample.
    """
    R = np.random.randint(0, n, size=(B, n))
    setup = np.zeros(shape=(B, n))
    setup_c = np.zeros(shape=(B, n))
    if param_grid is not None:
        delta = param_grid(n)
    for i in range(B):
        for j in range(n):
            if j in R[i]:
                setup[i, j] = 1
            else:
                setup_c[i, j] = 1
    m = len(delta)
    denom = np.sum(setup_c, axis=1)
    V = np.zeros(shape=(B, m))
    C = np.zeros(shape=(1, m))
    for i in range(m):
        cost = 0
        for j in range(B):
            ksi_validation = ksi[setup_c[j] == 1]
            ksi_train = ksi[R[j, :], :]
            x_ = solve(delta[i], c, b, ksi_train, alpha, d, n)
            P = np.dot(ksi_validation, x_) - b
            cost += np.dot(c, x_)
            p = len(P[P < 0]) / denom[j]
            V[j, i] = p
        C[0, i] = cost / B
    min_ind = _select_index(V, C, alpha, beta, B)
    delta_star = delta[min_ind]
    return solve(delta_star, c, b, ksi, alpha, d, n)


def sectioning(solve, delta, c, b, alpha, n, d, ksi, beta, n1, n2, K):
    """Sectioning: solve on ``n1`` training points, validate on ``K`` sections of ``n2``."""
    ksi_shuffle = np.random.permutation(ksi)
    ksi_train = ksi_shuffle[:n1]
    ksi_v = ksi_shuffle[n1:]
    m = len(delta)
    l = int(n2 / K)
    V = np.zeros(shape=(K, m))
    C = np.zeros(shape=(1, m))
    xx = np.zeros(shape=(d, m))
    for i in range(m):
        x_ = solve(delta[i], c, b, ksi_train, alpha, d, n1)
        xx[:, i] = x_
        for j in range(K):
            ksi_validation = ksi_v[j * l:(j + 1) * l]
            P = np.dot(ksi_validation, x_) - b
            p = len(P[P < 0]) / l
            V[j, i] = p
        C[0, i] = np.dot(c, x_)
    min_ind = _select_index(V, C, alpha, beta, K)
    return xx[:, min_ind]


def gaussian_supremum(solve, delta, c, b, alpha, n, d, ksi, beta, n1, n2):
    """Gaussian-supremum validator and baselines.

    Solves for every parameter in ``delta`` on ``n1`` training points, then
    validates on ``n2`` held-out points. Returns a ``d x 4`` array whose columns
    are the solutions chosen by the NGS, UNGS, NV and UG rules (in that order).
    """
    ksi_shuffle = np.random.permutation(ksi)
    ksi_train = ksi_shuffle[:n1]
    ksi_v = ksi_shuffle[n1:]
    xx = solve(delta, c, b, ksi_train, alpha, d, n1)
    C = np.dot(c, xx)
    res = np.zeros(shape=(d, 4))

    indifunc = (np.dot(ksi_v, xx) <= b).astype(int)
    Sigma_hat = np.cov(indifunc.T)
    sigma_hat = np.sqrt(np.diag(Sigma_hat))
    ind_nonzero = sigma_hat > 0
    P_hat = np.mean(indifunc, axis=0)

    # Univariate Gaussian (UG) margin.
    feasible_UG = P_hat >= 1 - alpha + norm.ppf(1 - beta) * sigma_hat / np.sqrt(n2)
    if np.sum(feasible_UG) > 0:
        xx_feasible = xx[:, feasible_UG]
        index1 = np.argmin(C[feasible_UG])
        X_UG = xx_feasible[:, index1]
    else:
        ind = np.argmax(P_hat - norm.ppf(1 - beta) * sigma_hat / np.sqrt(n2))
        X_UG = xx[:, ind]

    # Normalized (NGS) and Unnormalized (UNGS) Gaussian Supremum quantiles.
    if np.any(ind_nonzero):
        Sigma_hat = Sigma_hat[ind_nonzero][:, ind_nonzero]
        sim_num = 2000
        l = len(Sigma_hat)
        Z = np.random.multivariate_normal(np.zeros(l), Sigma_hat, sim_num)
        tmat = np.tile(sigma_hat[ind_nonzero], (sim_num, 1))
        q = np.amax(Z / tmat, axis=1)
        qu = np.amax(Z, axis=1)
        q_hat = np.percentile(q, (1 - beta) * 100)
        qu_hat = np.percentile(qu, (1 - beta) * 100)
    else:
        q_hat, qu_hat = 0, 0

    feasible_NGS = (P_hat >= 1 - alpha + q_hat * sigma_hat / np.sqrt(n2))
    if np.any(feasible_NGS):
        xx_feasible = xx[:, feasible_NGS]
        C_feasible = C[feasible_NGS]
        ind = np.argmin(C_feasible)
        X_NGS = xx_feasible[:, ind]
    else:
        ind = np.argmax(P_hat - q_hat * sigma_hat / np.sqrt(n2))
        X_NGS = xx[:, ind]

    feasible_UNGS = (P_hat >= 1 - alpha + qu_hat / np.sqrt(n2))
    if np.any(feasible_UNGS):
        xx_feasible = xx[:, feasible_UNGS]
        C_feasible = C[feasible_UNGS]
        ind = np.argmin(C_feasible)
        X_UNGS = xx_feasible[:, ind]
    else:
        ind = np.argmax(P_hat - qu_hat / np.sqrt(n2))
        X_UNGS = xx[:, ind]

    # Naive validation (NV) baseline: no statistical margin.
    feasible_NV = (P_hat >= 1 - alpha)
    if np.any(feasible_NV):
        xx_feasible = xx[:, feasible_NV]
        C_feasible = C[feasible_NV]
        ind = np.argmin(C_feasible)
        X_NV = xx_feasible[:, ind]
    else:
        ind = np.argmax(P_hat)
        X_NV = xx[:, ind]

    res[:, 0], res[:, 1], res[:, 2], res[:, 3] = X_NGS, X_UNGS, X_NV, X_UG
    return res


def run_experiment(solve, c, N, d, n, alpha, beta, K, B, n1, n2, delta, b,
                   dim=2, cv_grid=None, bs_grid=None):
    """Monte-Carlo comparison of the validation methods.

    Repeats the experiment ``N`` times on freshly sampled data and returns
    ``(mean_obj, feas_level, collect1, collect2)`` where ``collect1``/``collect2``
    hold the per-replication objective and (0/1) feasibility for each method.

    ``dim`` selects which methods are compared and hence how many columns the
    outputs have:

        dim == 2  cross_validation, bootstrapping
        dim == 4  gaussian_supremum only (NGS, UNGS, NV, UG)
        dim == 5  sectioning + gaussian_supremum (NGS, UNGS, NV, UG)

    ``cv_grid`` / ``bs_grid`` are the discrete parameter grids (see
    :func:`cross_validation` / :func:`bootstrapping`) used by SO/SAA at dim 2.
    """
    x = np.zeros(shape=(d, dim))
    collect1 = np.zeros(shape=(N, dim))
    collect2 = np.zeros(shape=(N, dim))
    for i in range(N):
        ksi = abs(np.random.normal(size=(n, d)))
        if dim == 2:
            x[:, 0] = cross_validation(solve, delta, c, b, alpha, n, d, ksi, beta, K, param_grid=cv_grid)
            x[:, 1] = bootstrapping(solve, delta, c, b, alpha, n, d, ksi, beta, B, param_grid=bs_grid)
        elif dim == 5:
            x[:, 0] = sectioning(solve, delta, c, b, alpha, n, d, ksi, beta, n1, n2, K)
            x[:, 1:] = gaussian_supremum(solve, delta, c, b, alpha, n, d, ksi, beta, n1, n2)
        elif dim == 4:
            x[:, 0:] = gaussian_supremum(solve, delta, c, b, alpha, n, d, ksi, beta, n1, n2)
        else:
            raise ValueError('dim must be 2, 4 or 5')
        metric1 = np.dot(c, x)
        ksi_ = abs(np.random.normal(size=(1000000, d)))
        metric2 = (np.sum((np.dot(ksi_, x) <= b).astype(int), 0) / 1000000 >= 1 - alpha).astype(int)
        collect1[i] = metric1
        collect2[i] = metric2
    mean_obj, feas_level = np.mean(collect1, 0), np.mean(collect2, 0)
    return mean_obj, feas_level, collect1, collect2
