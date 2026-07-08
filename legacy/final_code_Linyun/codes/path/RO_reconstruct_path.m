function [solution] = RO_reconstruct_path(c,b,alpha,beta,data1,para)
% c: cost coefficient
% b: threshold
% alpha: safety level
% beta: confidence level
% data1: phase one data set, each row is an observation of xi
% para: a sequence of discrete parameter values

n = length(c);
N1 = size(data1,1);

order1 = binoinv(1-beta,N1,1-alpha);
if order1 > N1 - 1
    error("too small N1");
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
    t_select2 = para;% the parameter
    mesh_size = length(t_select2);
    solution = zeros(n,mesh_size);
    
    for k = 1:mesh_size
        cvx_begin quiet
        variable p(1)
        minimize((c'*x_ro1)*p)
        subject to
        p*t_select2(k) <= b
        p >= 0
        cvx_end
        
        
        if strcmp(cvx_status,'Solved')
            solution(:,k) = p*x_ro1;
        else
            solution(:,k) = NaN(n,1);
        end
    end
else
    solution = NaN(n,length(para));
end
end