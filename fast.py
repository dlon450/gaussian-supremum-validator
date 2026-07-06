# so_fast.py
import numpy as np
import pandas as pd
from gurobipy import Model, GRB, quicksum

# ---------------------------
# FAST solver (mirrors MATLAB)
# ---------------------------
def SG_fast(c, b, data, N1, output_flag=0):
    """
    Scenario Optimization (FAST)
      Phase 1:  min c^T x      s.t. data1 x <= b
      Phase 2:  min c^T((1-a)x + a x0)  s.t. data2((1-a)x + a x0) <= b, 0<=a<=1
    Here x0 = 0, so Phase 2 reduces to 1D LP in 'a'.
    """
    c = np.asarray(c, dtype=float).ravel()
    n_obs, d = data.shape
    assert N1 > 0 and N1 < n_obs, "N1 must be in (0, number of rows of data)"
    data1 = np.asarray(data[:N1, :], dtype=float)
    data2 = np.asarray(data[N1:, :], dtype=float)

    # ---------
    # Phase 1
    # ---------
    m1 = Model()
    m1.Params.OutputFlag = output_flag
    x = m1.addVars(d, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="x")
    m1.setObjective(quicksum(float(c[j]) * x[j] for j in range(d)), GRB.MINIMIZE)

    # data1 * x <= b (row-wise)
    for i in range(N1):
        row = data1[i]
        m1.addConstr(quicksum(float(row[j]) * x[j] for j in range(d)) <= float(b))

    m1.optimize()
    if m1.Status != GRB.OPTIMAL:
        return np.full(d, np.nan)

    x_fast = np.array([x[j].X for j in range(d)])
    x0 = np.zeros_like(x_fast)  # feasible anchor

    # ---------
    # Phase 2  (only alpha variable)
    # ---------
    # Constraint: data2 * ((1 - a) * x_fast + a * x0) <= b
    # Since x0=0: (1 - a) * (data2 @ x_fast) <= b
    # Each row i: -(data2_i @ x_fast)*a <= b - (data2_i @ x_fast)
    dx = data2 @ x_fast

    m2 = Model()
    m2.Params.OutputFlag = output_flag
    a = m2.addVar(lb=0.0, ub=1.0, vtype=GRB.CONTINUOUS, name="alpha")

    for i in range(data2.shape[0]):
        coef = -float(dx[i])  # coefficient on 'a'
        rhs = float(b - dx[i])
        m2.addConstr(coef * a <= rhs)

    # Objective: c^T((1 - a) x_fast) = c^T x_fast - a * (c^T x_fast)
    cTx = float(np.dot(c, x_fast))
    m2.setObjective(cTx - a * cTx, GRB.MINIMIZE)
    m2.optimize()

    if m2.Status != GRB.OPTIMAL:
        return np.full(d, np.nan)

    a_star = a.X
    sol = (1.0 - a_star) * x_fast + a_star * x0  # = (1-a)*x_fast since x0=0
    return sol


# ---------------------------
# Eval helper (objective & feasibility on validation)
# ---------------------------
def evaluate_solution(c, b, alpha, x, val_data):
    """Return (objective, P_hat, feasible_boolean) using validation data."""
    obj = float(np.dot(c, x))
    P_hat = float(np.mean((val_data @ x) <= b))
    feas = (P_hat >= 1.0 - alpha)
    return obj, P_hat, feas


# ---------------------------
# Experiment runner (same grid as your original)
# ---------------------------
def run_fast_experiments(seed=2):
    np.random.seed(seed)

    # same base parameters
    d_list = [10]  # can extend if you wish
    n_list = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    alpha = 0.05
    c_base = -1.0
    N_reps = 100
    output = "fast_results.csv"

    all_rows = []

    for d in d_list:
        c = c_base * np.ones(d)
        b = 0.833 * d  # same as your original

        for n in n_list:
            print(f"[FAST] Data Size: {n}")

            # training shares: 10%, 20%, ..., 90%
            for frac in range(1, 10):
                train_pct = 10 * frac
                N1 = n * frac // 10
                N2 = n - N1

                objs, feas_flags = [], []

                for rep in range(N_reps):
                    # generate fresh data each repetition (mirrors your style)
                    data = np.abs(np.random.normal(size=(n, d)))

                    # FAST uses ordered split in MATLAB; but we can also shuffle for robustness.
                    # To mirror the MATLAB code exactly, DO NOT shuffle:
                    # data1 = data[:N1], data2 = data[N1:]
                    sol = SG_fast(c, b, data, N1, output_flag=0)

                    # evaluate on Phase-2 set (validation) as in MATLAB
                    val_data = data[N1:, :]
                    obj, P_hat, feas = evaluate_solution(c, b, alpha, sol, val_data)
                    objs.append(obj)
                    feas_flags.append(feas)

                mean_obj = float(np.mean(objs))
                feas_level = float(np.mean(feas_flags)) * 100.0  # percent feasible

                print(f"  n1/n = {train_pct:>3d}% | mean obj = {mean_obj: .4f} | feasibility = {feas_level:5.1f}%")

                all_rows.append({
                    "Algorithm": "FAST",
                    "Data Size": n,
                    "Training Percentage": train_pct,
                    "Objective": mean_obj,
                    "Feasibility": feas_level,
                })

    df = pd.DataFrame(all_rows)
    df.to_csv(output, index=False)
    print(f"\nSaved results to {output}")
    return df


if __name__ == "__main__":
    run_fast_experiments()
