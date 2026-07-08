function [solution] = CCP_RO_LCX(para, c, b, alpha, data)
% para: a sequence of parameter values specifying the size of the uncertainty set
% c: cost vector in the objective
% b: right-hand threshold in the constraint
% data: data matrix, each row is an observation

d = length(c);
n = size(data, 1);
mesh_size = length(para);
solution = zeros(d, mesh_size);

% cvx_solver gurobi_2
for k = 1:mesh_size
    cvx_begin quiet
    
        variables x_ro(d) f(d) g(d) h(d)
        variable w(n) nonnegative
        variable y(n) nonnegative
        variable z(n) nonnegative
        variable tau nonnegative
        variable theta nonnegative
        variables alfa(3) e(3)
        
        minimize (c'*x_ro)
        
        subject to
        1/alpha*tau - theta + para(k)*sum(e) + 2*para(k)*sum(f) + para(k)*sum(g) <= b
        -e <= alfa <= e
        -f <= h <= f
        -g <= x_ro + h <= g
        -theta + tau + alfa(1) + alfa(2) == sum(w)/n + sum(y)/n + sum(z)/n
        alfa(1) - data*h <= w
        alfa(2) + data*h <= y
        alfa(3) + data*h + data*x_ro <= z
        0 <= x_ro <= 1
        
    cvx_end
    
    solution(:,k) = x_ro;
end

end