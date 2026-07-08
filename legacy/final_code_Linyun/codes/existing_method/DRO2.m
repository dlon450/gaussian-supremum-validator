function [solution] = DRO2(c,b,alpha,beta,data)
% c: cost coefficient
% b: threshold
% alpha: safety level
% beta: confidence level
% data: whole data set, each row is an observation of xi
n = length(c);
triu_ind=triu(true(n));

col_num=repmat(1:n,n,1);
row_num=col_num';

row_ind_use=row_num(triu_ind);
col_ind_use=col_num(triu_ind);


n_phi_input=n+n*(n+1)/2;
x_var=sym('x',[n_phi_input,1]);
% syms a b
multi_mat=x_var(1:n)*x_var(1:n).';
phi_origin=[x_var];
phi_origin(1:n)=sqrt(alpha/(1-alpha))*phi_origin(1:n);
phi_origin(n+1:end)=phi_origin(n+1:end)-multi_mat(triu_ind);

jacob_phi_origin=jacobian(phi_origin, x_var);
phi_fun_temp(x_var)=jacob_phi_origin;
phi_fun=matlabFunction(phi_fun_temp,'Vars',{x_var.'});



N_data = size(data,1);
dataset_for_theta=[data, data(:,row_ind_use).*data(:,col_ind_use)];

theta_estimate=mean(dataset_for_theta);
new_sigma=cov(dataset_for_theta,1);

value_jacob= phi_fun(theta_estimate);

V_est=value_jacob*new_sigma*value_jacob';
% V_est_inv=pinv(V_est);

mu_hat=mean(data);
Sigma_hat=cov(data,1);
% R_inv=inv(sqrtm(V_est_inv));
R_inv=sqrtm(V_est);
rho = sqrt(chi2inv(1-beta,n_phi_input)/N_data);
% rho = 0.5;
tilde_c=-b;

svec_operator=sqrt(2)*ones(n)-sqrt(2)*eye(n)+eye(n);
svec_multiplier=svec_operator(triu_ind);
A=[eye(n),zeros(n,n_phi_input-n);
    zeros(n_phi_input-n,n),diag(svec_multiplier)];

cvx_begin quiet
    variable x_data_dro(n)
    variable W(n,n)
    variable vec_q(n_phi_input)
    variable eta
    minimize (c'*x_data_dro)
    subject to
        sqrt(alpha/(1-alpha))*mu_hat*x_data_dro + trace (Sigma_hat*W) + rho *norm((A*R_inv)'*(vec_q)) +sqrt(alpha/(1-alpha))* tilde_c+eta/4  <= 0
        vec_q == [x_data_dro;W(triu_ind).*svec_multiplier]
        [W,x_data_dro;x_data_dro',eta] == semidefinite(n+1)
cvx_end

if strcmp(cvx_status,'Solved')
    solution = x_data_dro;
else
    solution = NaN(n,1);
end
end