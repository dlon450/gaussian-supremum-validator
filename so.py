"""Experiment: scenario optimization (SO).

Compares cross validation against bootstrapping (dim == 2) for the scenario-
optimization solver across a range of data sizes and validation-fold counts.

The SO parameter is a discrete scenario count, so the candidate grid is
regenerated inside the validators via ``cv_grid`` / ``bs_grid`` (the ``delta``
argument passed here is a placeholder and is overridden).

Run with:  python so.py
"""

import numpy as np
import pandas as pd
import time

from gsv.formulations import CCP_SO
from gsv.validators import run_experiment

np.random.seed(2)
d = 10
n = 500
n1, n2 = n // 2, n // 2
data = np.abs(np.random.normal(size=(n, d)))  # kept: advances the RNG stream to match published results
delta1 = [3 / 2**i for i in range(10)]
delta2 = [delta1[-1] / 1.2**i for i in range(1,16)]
delta = np.array(delta1 + delta2)

alpha = 0.05
beta = 0.1
c = -np.ones(d)
b = .833 * d
N = 100
K = B = 10

# Discrete scenario-count grids used by the SO validators (see validators.py).
so_grid = lambda ne: np.arange(1, ne + 1)

dim = 2
column_names = ['Cross Validation','Bootstrapping','Sectioning','NGS','UNGS','NV','UG']
column_names_ = column_names[2:] if dim == 5 else column_names[3:] if dim == 4 else column_names[:2]
output = 'experiment_results_so.csv'

############ RO ###########
# delta = 45 * np.arange(1, 26) / 25

for n in [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]:
    print('Data Size:', n)
    start = time.time()
    for i, n1 in enumerate(range(n // 10, n, n // 10)):
        n2 = n - n1
    for K in [3, 5, 10]:
        B = K

        delta_SO = np.arange(1, n1 + 1, 1)
        result = run_experiment(CCP_SO,c,N,d,n,alpha,beta,K,B,n1,n2,delta=delta_SO,b=b,dim=dim,cv_grid=so_grid,bs_grid=so_grid)
        print(result[0], result[1], time.time() - start)

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
