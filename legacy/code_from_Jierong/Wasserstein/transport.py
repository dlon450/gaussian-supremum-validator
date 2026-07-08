from ortools.graph import pywrapgraph

def distance(cost, n, precision):
    mean_cost = sum(cost)/len(cost)
    if mean_cost == 0:
        return 0
    coeff = max(1/precision, mean_cost)/mean_cost
    n = round(n)

    assignment = pywrapgraph.LinearSumAssignment()
    for row in range(n):
        for col in range(n):
            assignment.AddArcWithCost(row, col, round(cost[col * n + row]*coeff))

    solve_status = assignment.Solve()
    if solve_status == assignment.OPTIMAL:
        return assignment.OptimalCost()/coeff/n
    return 'unsolved'