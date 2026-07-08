function [solution] = CCP_SAA(para, c, b, data)
% para: a sequence of parameter values specifying the tolerance level
% c: cost vector in the objective
% b: right-hand threshold in the constraint
% data: data matrix, each row is an observation

d = length(c);
n = size(data, 1);
M = ceil(abs(b) + max(sum(abs(data), 2)));
num_constr = ceil(n*para); % transform the tolerance level into the number of satisfied constraints
mesh_size = length(num_constr);
solution = zeros(d, mesh_size);

% cvx_solver mosek
for k = 1:mesh_size
    if k > 1 && num_constr(k) == num_constr(k - 1)
        solution(:, k) = solution(:, k - 1);
    else
        cvx_begin
    
            variable x_SAA(d)
            variable z(n) binary
        
            minimize (c'*x_SAA)
        
            subject to
            data*x_SAA <= b + M*(1-z)
            sum(z) >= num_constr(k)
            0 <= x_SAA <= 1
        
        cvx_end
    
        solution(:,k) = x_SAA;
    end
end

end