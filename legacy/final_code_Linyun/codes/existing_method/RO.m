function [solution, t_select] = RO(c,b,alpha,beta,data,N1)
% c: cost coefficient
% b: threshold
% alpha: safety level
% beta: confidence level
% data: whole data set, each row is an observation of xi
% N1: sample size for phase one

n = length(c);

data1 = data(1:N1,:);
data2 = data(N1+1:end,:);
N2 = size(data2,1);
order2 = binoinv(1-beta,N2,1-alpha); % estimated quantile
if order2 > N2 - 1
    error("too small N2");
end

% phase one
mu_hat = mean(data1);
sigma_hat = cov(data1);


% phase two
mu_hat_span = repmat(mu_hat,N2,1);
t_value = diag((data2-mu_hat_span)*sigma_hat^-1*(data2-mu_hat_span)');
t_sorted = sort(t_value);
t_select = t_sorted(order2 + 1); % the selected parameter


sigma_hat_rt = sqrtm(sigma_hat);
cvx_begin quiet
variable x_ro(n)
minimize(c'*x_ro)
subject to
mu_hat*x_ro + sqrt(t_select)*norm(sigma_hat_rt*x_ro) - b <= 0
cvx_end


if strcmp(cvx_status,'Solved')
    solution = x_ro;
else
    solution = NaN(n,1);
end

end