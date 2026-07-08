function [solution, t_select] = RO_reconstruct(c,b,alpha,beta,data,N1)
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
order1 = binoinv(1-beta,N1,1-alpha);
order2 = binoinv(1-beta,N2,1-alpha); % estimated quantile
if order1 > N1 - 1
    error("too small N1");
end
if order2 > N2 - 1
    error("too small N2");
end

% estimate mean and covariance
mu_hat = mean(data1);
sigma_hat = cov(data1);


% phase one
mu_hat_span = repmat(mu_hat,N1,1);
t_value = diag((data1-mu_hat_span)*sigma_hat^-1*(data1-mu_hat_span)');
t_sorted = sort(t_value);
t_select = t_sorted(order1 + 1); % the selected parameter


sigma_hat_rt = sqrtm(sigma_hat);
cvx_begin quiet
    variable x_ro1(n)
    minimize(c'*x_ro1)
    subject to
        mu_hat*x_ro1 + sqrt(t_select)*norm(sigma_hat_rt*x_ro1) - b <= 0
cvx_end

% phase two
if strcmp(cvx_status,'Solved')
    t_value2 = data2*x_ro1;
    t_sorted2 = sort(t_value2);
    t_select2 = t_sorted2(order2 + 1);% the selected parameter
    
    
    cvx_begin quiet
        variable p(1)
        minimize((c'*x_ro1)*p)
        subject to
            p*t_select2 <= b
            p >= 0
    cvx_end
    
    
    if strcmp(cvx_status,'Solved')
        solution = p*x_ro1;
    else
        solution = NaN(n,1);
    end
else
    solution = NaN(n,1);
end
end