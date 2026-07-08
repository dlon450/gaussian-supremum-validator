function solution = CCP_DRO_WassersteinCVaR(para, c, b, alpha, data)
% para: a sequence of delta values
% c: cost vector in the objective
% b: right-hand threshold in the constraint
% alpha: tolerance level in the original chance constraint
% data: data matrix, each row is an observation


[n, d] = size(data);
m = length(para);
solution = zeros(d, m);

% cvx_solver gurobi_2
for k = 1:m
    cvx_begin quiet
        
        variable x(d)
        variable v nonnegative
        variable r nonnegative
        variable z(n) nonnegative

        minimize (c.'*x)

        subject to
        para(k) * v + sum(z) / n <= alpha * r
        b - data * x + z >= r
        v >= x
        0 <= x <= 1

    cvx_end

    solution(:, k) = x;
end

end

