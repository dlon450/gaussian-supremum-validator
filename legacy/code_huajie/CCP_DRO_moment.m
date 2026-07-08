function [solution] = CCP_DRO_moment(para, c, b, alpha, data)
% para: a sequence of parameter values specifying the size of the uncertainty set
% c: cost vector in the objective
% b: right-hand threshold in the constraint
% alpha: tolerance level in the original chance constraint
% data: data matrix, each row is an observation


%% construct Jacobian matrix
d = length(c);
d_aug = d * (d + 1) / 2;
mu_hat = mean(data);
grad_A = sqrt(alpha / (1 - alpha)) * speye(d);
grad_C = speye(d_aug);
grad_B = sparse(0, d);

% construct the lower-left block, grad_B, of the Jacobian matrix
for col = 1:d
    block_size = d - col + 1;
    new_block = sparse(1:block_size, col:d, -mu_hat(col), block_size, d) + sparse(1:block_size, col, -mu_hat(col:d), block_size, d);
    grad_B = [grad_B; new_block];
end
grad = [grad_A, sparse(d, d_aug); grad_B, grad_C];


%% compute covariance matrix
tril_ind = tril(true(d));
[row_ind, col_ind] = find(tril_ind);
data_moment = [data, data(:, row_ind) .* data(:, col_ind)];
moment_sigma = cov(data_moment, 1);
V_est = grad * moment_sigma * grad.';


%% solver

% prepare data
Sigma_hat = cov(data, 1);
tilde_c = -sqrt(alpha / (1 - alpha)) * b;
sqrt_cov = sqrtm(V_est);
svec_operator = sqrt(2) * ones(d) + (1-sqrt(2)) * eye(d);
svec_multiplier = svec_operator(tril_ind);
A = sparse(1 : d+d_aug, 1 : d+d_aug, [ones(d,1); svec_multiplier]);

rho = para;
mesh_size = length(rho);
solution = zeros(d, mesh_size);

% cvx_solver mosek
for k = 1:mesh_size
    cvx_begin quiet
    
        variable x_dro(d)
        variable W(d,d)
        variable vec_q(d + d_aug)
        variable eta
        
        minimize (c'*x_dro)
        
        subject to
        sqrt(alpha/(1-alpha))*mu_hat*x_dro + trace(Sigma_hat*W) + rho(k)*norm((A*sqrt_cov).'*vec_q) + tilde_c + eta/4 <= 0
        vec_q == [x_dro; W(tril_ind).*svec_multiplier]
        [W,x_dro; x_dro.',eta] == semidefinite(d+1)
        0 <= x_dro <= 1
        
    cvx_end

    solution(:,k) = x_dro;
end

end