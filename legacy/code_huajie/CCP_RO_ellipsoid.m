function [solution] = CCP_RO_ellipsoid(para, c, b, data)
% para: a sequence of parameter values specifying the size of the uncertainty set
% c: cost vector in the objective
% b: right-hand threshold in the constraint
% data: data matrix, each row is an observation

d = length(c);

mu_hat = mean(data);
sigma_hat = cov(data);
sigma_hat_rt = sqrtm(sigma_hat);

radius = para;
mesh_size = length(radius);
solution = zeros(d, mesh_size);

% cvx_solver mosek
for k = 1:mesh_size
    cvx_begin quiet
    
        variable x_ro(d)
        
        minimize (c'*x_ro)
        
        subject to
        mu_hat*x_ro + sqrt(radius(k))*norm(sigma_hat_rt*x_ro) - b <= 0
        0 <= x_ro <= 1
        
    cvx_end
    
    solution(:,k) = x_ro;
end

end