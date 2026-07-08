function [solution] = True_sol(c,b,alpha,mu,sigma)
% c: cost coefficient
% b: threshold
% alpha: safety level
% mu: true mean of xi
% sigma: true covariance matrix of xi

n = length(c);
phi_quantile = norminv(1-alpha, 0, 1);
sigma_rt = sqrtm(sigma);

cvx_begin quiet
    variable x(n)
    minimize(c.'*x)
    subject to
        phi_quantile*norm(sigma_rt*x) + mu(:).'*x - b <= 0
cvx_end


if strcmp(cvx_status,'Solved')
    solution = x;
else
    solution = NaN(n,1);
end
end