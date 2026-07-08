function [solution] = SG_path(c,b,data1)
% c: cost coefficient
% b: threshold
% data1: phase one data set, each row is an observation of xi
n = length(c);              % dimension
N_data = size(data1,1);
% N_data = 15;

solution = zeros(n,N_data);

for k = 1:N_data    
    cvx_begin quiet
    variable x_gen(n)
    minimize(c.'*x_gen)
    subject to
    data1(1:k,:)*x_gen <= b
    cvx_end
    
    
    if strcmp(cvx_status,'Solved')
        solution(:,k) = x_gen;
    else
        solution(:,k) = NaN(n,1);
    end
end
end