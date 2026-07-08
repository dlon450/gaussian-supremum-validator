function x = Optimize(delta, c, b, alpha,n, M, d,ksi)
prob = optimproblem;
y = optimvar('y',n,'Type','integer','LowerBound',0,'UpperBound',1);
v = optimvar('v',1,'LowerBound',0);
z = optimvar('z',n,'LowerBound',0);
r = optimvar('r',1,'LowerBound',0);
s = optimvar('s',n,'LowerBound',0);
x = optimvar('x',d,'LowerBound',0,'UpperBound',1);
prob.Objective = c' * x;
p1 = z + s - r;
p2 = delta * v + sum(z) / n - alpha * r;
p3 = b - s - ksi * x + M * (1 - y);
p4 = M * y - s;
p5 = v - x;
prob.Constraints.const1 = p1 >= 0;
prob.Constraints.const2 = p2 <= 0;
prob.Constraints.const3 = p3 >= 0;
prob.Constraints.const4 = p4 >= 0;
prob.Constraints.const5 = p5 >= 0;
[sol,fval] = solve(prob);
x = sol.x;
end