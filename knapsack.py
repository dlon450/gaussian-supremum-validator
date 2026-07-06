import gurobipy as gp
from gurobipy import GRB
import numpy as np
from array import array


def DRO_solver(para, c, b, alpha, data, kappa=1, warm_start=False):
    # Input processing
    d = len(c)
    n = round(len(data)/d)
    c = list(c)
    data = np.reshape(data, (n, d), order = 'F')
    rM = np.maximum(-data, 0)
    tM = np.maximum(data, 0)
    uM = b - tM.sum(axis=1)
    kappa = min(kappa, n)
    if isinstance(para, float):
        para = array('d', [para])

    # Solve the optimization instances, in the descending order of the radiuses
    ind_para = np.argsort(para)[::-1]
    solution = list()
    feasible_prev = False # An indicator, True if a start point is available from the previous instance, False otherwise
    for idx in ind_para:
        # Create an empty model with lazy constraints
        opt = gp.Model()
        opt.Params.OutputFlag = 0 # Keep the Gurobi log silent
        opt.Params.LazyConstraints = 1 # Use lazy constraints, to be added via the callback function upon violation

        # Create variables
        x = opt.addVars(d, obj = c, vtype = GRB.BINARY)
        w = opt.addVars(d, lb = 0, ub = 1, vtype = GRB.CONTINUOUS)
        z = opt.addVars(n, lb = 0, vtype = GRB.CONTINUOUS)
        r = opt.addVar(lb = 0, vtype = GRB.CONTINUOUS)
        opt.ModelSense = GRB.MINIMIZE

        # Add constraints
        opt.addConstr(gp.quicksum(z) <= alpha*n*r - para[idx]*n)
        opt.addConstrs(x[i] + w[i] == 1 for i in range(d))

        # Use warm start for the next instance
        if warm_start and feasible_prev:
            for i in range(d):
                x[i].Start = x_prev[i]
                w[i].Start = w_prev[i]
            for i in range(n):
                z[i].Start = z_prev[i]
            r.Start = r_prev

        # Pass data to the callback function, and call the solver
        opt._x = x
        opt._w = w
        opt._z = z
        opt._r = r
        opt._rM = rM
        opt._tM = tM
        opt._uM = uM
        opt._kappa = kappa
        opt.optimize(callback)

        # Collect solution
        if opt.Status == GRB.OPTIMAL:
            feasible_prev = True
            x_prev = [x[i].X for i in range(d)]
            w_prev = [w[i].X for i in range(d)]
            z_prev = [z[i].X for i in range(n)]
            r_prev = r.X
            if opt.ObjVal <= 0 or b < 0:
                solution.append(x_prev)
            else:
                solution.append([0]*d)
        elif opt.Status == GRB.INFEASIBLE:
            feasible_prev = False
            if b >= 0:
                solution.append([0]*d)
            else:
                return 'infeasible'
        else:
            return 'unsolved'

    ind_para_rev = np.argsort(ind_para)
    solution_ordered = array('d')
    for idx in ind_para_rev:
        solution_ordered.extend(solution[idx])
    return solution_ordered


def callback(model, where):
    if where == GRB.Callback.MIPSOL:
        # Retrieve the current best solution
        x_hat = model.cbGetSolution(model._x)
        w_hat = model.cbGetSolution(model._w)
        z_hat = model.cbGetSolution(model._z)
        r_hat = model.cbGetSolution(model._r)
        d = len(x_hat)
        x_hat = np.array([x_hat[i] for i in range(d)])
        w_hat = np.array([w_hat[i] for i in range(d)])
        z_hat = np.array([z_hat[i] for i in range(len(z_hat))])

        # Order lazy constraints by amount of violation
        violation = -np.maximum(np.dot(model._rM, x_hat) + np.dot(model._tM, w_hat) + model._uM, 0) - z_hat + r_hat
        idx_constr = violation.argsort()[::-1]

        # Add the most violated lazy constraints
        idx_var = np.argsort(np.concatenate((x_hat, w_hat)))[::-1]
        k = 0
        while k < model._kappa and violation[idx_constr[k]] > 1e-6:
            coeff = np.concatenate((model._rM[idx_constr[k], :], model._tM[idx_constr[k], :]))
            cut = model._uM[idx_constr[k]]
            cut_pos = -max(cut, 0)

            new_constr = gp.LinExpr()
            if cut_pos != 0:
                new_constr.addConstant(cut_pos)
            for idx in idx_var:
                cut += coeff[idx]
                cut_pos_new = -max(cut, 0)
                if cut_pos_new != cut_pos:
                    new_constr.addTerms(cut_pos_new - cut_pos, model._x[idx] if idx < d else model._w[idx-d])
                    cut_pos = cut_pos_new

            model.cbLazy(new_constr <= model._z[idx_constr[k]] - model._r)
            k += 1


d = 20
n = 200
para = [0.001, 0.01, 0.04, 0.1]
c = -np.linspace(start=1, stop=2, num=d)
b = 0.8*d
alpha = 0.05
data = np.abs(np.random.normal(size = d*n))

solution = DRO_solver(para, c, b, alpha, data)
print(np.dot(c, solution[:d]))
print(np.dot(c, solution[d:2*d]))
print(np.dot(c, solution[2*d:3*d]))
print(np.dot(c, solution[3*d:]))



# def callback(model, where):
#     if where == GRB.Callback.MIPNODE:
#         x_hat = model.cbGetNodeRel(model._x)
#         w_hat = model.cbGetNodeRel(model._w)
#         z_hat = model.cbGetNodeRel(model._z)
#         r_hat = model.cbGetNodeRel(model._r)
#     elif where == GRB.Callback.MIPSOL:
#         x_hat = model.cbGetSolution(model._x)
#         w_hat = model.cbGetSolution(model._w)
#         z_hat = model.cbGetSolution(model._z)
#         r_hat = model.cbGetSolution(model._r)
#     else:
#         return

#     d = len(x_hat)
#     x_hat = np.array([x_hat[i] for i in range(d)])
#     w_hat = np.array([w_hat[i] for i in range(d)])
#     z_hat = np.array([z_hat[i] for i in range(len(z_hat))])

#     violation = -np.maximum(np.dot(model._rM, x_hat) + np.dot(model._tM, w_hat) + model._uM, 0) - z_hat + r_hat
#     idx_constr = violation.argsort()[::-1]

#     xw_val = np.concatenate((x_hat, w_hat))
#     idx_var = xw_val.argsort()[::-1]

#     k = 0
#     while k < model._kappa and violation[idx_constr[k]] > 0:
#         coeff = np.concatenate((model._rM[idx_constr[k], :], model._tM[idx_constr[k], :]))
#         cut = model._uM[idx_constr[k]]
#         cut_pos = -max(cut, 0)

#         new_constr = gp.LinExpr()
#         if cut_pos != 0:
#             new_constr.addConstant(cut_pos)
        
#         coeff_list = []
#         var_list = []
#         for var in range(2*d):
#             cut += coeff[idx_var[var]]
#             cut_pos_new = -max(cut, 0)
#             if cut_pos_new != cut_pos:
#                 var_list.append(model._x[idx_var[var]] if idx_var[var] < d else model._w[idx_var[var] - d])
#                 coeff_list.append(cut_pos_new - cut_pos)
#                 cut_pos = cut_pos_new
#         new_constr.addTerms(coeff_list, var_list)

#         model.cbLazy(new_constr <= model._z[idx_constr[k]] - model._r)
#         k += 1
    