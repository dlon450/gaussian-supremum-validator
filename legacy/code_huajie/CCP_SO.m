function [solution] = CCP_SO(para, c, b, data)
% para: a sequence of integers representing the number of scenarios to use
% c: cost vector in the objective
% b: right-hand threshold in the constraint
% data: data matrix, each row is an observation

d = length(c);
mesh_size = length(para);
solution = zeros(d, mesh_size);

% cvx_solver gurobi_2
for k = 1:mesh_size
    cvx_begin quiet
    
        variable x_gen(d)
        
        minimize(c.'*x_gen)
        
        subject to
        data(1:para(k), :)*x_gen <= b
        0 <= x_gen <= 1
        
    cvx_end
    
    solution(:,k) = x_gen;
end

end