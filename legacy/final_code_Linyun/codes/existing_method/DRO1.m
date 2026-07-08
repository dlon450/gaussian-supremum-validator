function [solution] = DRO1(c,b,alpha,mu,sigma)
% c: cost coefficient
% b: threshold
% alpha: safety level
% mu: true mean of xi
% sigma: true covariance matrix of xi
n = length(c);
sigma_rt = sqrtm(sigma);

cvx_begin quiet
    variable x(n)
    minimize(c.'*x)
    subject to
        sqrt((1-alpha)/alpha)*norm(sigma_rt*x) + mu(:).'*x - b <= 0
cvx_end


if strcmp(cvx_status,'Solved')
    solution = x;
else
    solution = NaN(n,1);
end
end