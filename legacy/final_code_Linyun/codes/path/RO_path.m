function [solution] = RO_path(c,b,data1,para)
% c: cost coefficient
% b: threshold
% data1: phase one data set, each row is an observation of xi
% para: a sequence of discrete parameter values

n = length(c);


mu_hat = mean(data1);
sigma_hat = cov(data1);
sigma_hat_rt = sqrtm(sigma_hat);


t_select = para; % the parameter
mesh_size = length(t_select);
solution = zeros(n,mesh_size);


for k = 1:mesh_size
    cvx_begin quiet
    variable x_ro(n)
    minimize(c'*x_ro)
    subject to
    mu_hat*x_ro + sqrt(t_select(k))*norm(sigma_hat_rt*x_ro) - b <= 0
    cvx_end
    
    
    if strcmp(cvx_status,'Solved')
        solution(:,k) = x_ro;
    else
        solution(:,k) = NaN(n,1);
    end
end
end