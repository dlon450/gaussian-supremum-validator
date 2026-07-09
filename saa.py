"""Experiment: sample-average approximation (SAA).

Compares cross validation against bootstrapping (dim == 2) for the SAA solver
across validation-fold counts.

The SAA parameter is a discrete count of satisfied constraints, so the candidate
grid is regenerated inside the validators via ``cv_grid`` / ``bs_grid`` (the
``delta`` argument passed here is a placeholder and is overridden).

Run with:  python saa.py
"""

import numpy as np
import pandas as pd

from gsv.formulations import CCP_SAA
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

# Discrete constraint-count grids used by the SAA validators (see validators.py).
saa_grid = lambda ne: np.arange(np.ceil(ne * (1 - alpha)).astype(int), ne + 1)

dim = 2
column_names = ['Cross Validation','Bootstrapping','Sectioning','NGS','UNGS','NV','UG']
column_names_ = column_names[2:] if dim == 5 else column_names[3:] if dim == 4 else column_names[:2]
output = 'experiment_results_saa.csv'

for n in [500]:
    print('Data Size:', n)
    for i, n1 in enumerate(range(n // 10, n, n // 10)):
        n2 = n - n1
    for K in [3, 5, 10]:
        B = K

        ########### SAA ##########
        delta = np.arange(np.ceil(n1 * (1 - alpha)).astype(int), n1 + 1)
        if len(delta) < 2: delta = np.repeat(delta[0], 2)
        result = run_experiment(CCP_SAA,c,N,d,n,alpha,beta,K,B,n1,n2,delta=delta,b=b,dim=dim,cv_grid=saa_grid,bs_grid=saa_grid)
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
