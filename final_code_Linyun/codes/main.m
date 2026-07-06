%% settings
cd F:\GoogleDrive\OhWork\Prof.Lam\NewFrameWorkforChanceConstraintProgram\simulations\codes
addpath(genpath(pwd));
load('setting');

save('setting', 'xi_mean', 'xi_cov', 'c_exp', 'b_exp', 'alpha', 'sol_true',...
    'beta_exp', 'objv_true', 'n1', 'N', 'num_exp', 'sim_num');

%% Non-Data-Driven
xi_mean = zeros(10, 1);
xi_cov = rand(10);
xi_cov = xi_cov' * xi_cov;
c_exp = rand(10, 1);
b_exp = 10 * rand(1);
alpha = 0.1;
sol_true = True_sol(c_exp, b_exp, alpha, xi_mean, xi_cov);
objv_true = c_exp' * sol_true;

%Safe Convex Approximation
sol_SCA = SCA(c_exp, b_exp, alpha, xi_mean, xi_cov);
objv_SCA = c_exp' * sol_SCA;
%Distributionally Robust Optimization
sol_DRO = DRO1(c_exp, b_exp, alpha, xi_mean, xi_cov);
objv_DRO = c_exp' * sol_DRO;

save('nondtdriven', 'sol_SCA', 'sol_DRO', 'objv_SCA', 'objv_DRO');
load('nondtdriven');
%% Data-Driven
beta_exp = 0.05;
n1 = 100;
N = 200;
num_exp = 1000;
sim_num = 2000;

%DRO
result_DRO = expDRO(c_exp, b_exp, alpha, beta_exp, xi_mean, xi_cov, n1, N, num_exp, sim_num);
result_DRO2 = expDRO(c_exp, b_exp, alpha, beta_exp, xi_mean, xi_cov, 250, 500, num_exp, sim_num);

%RO
result_RO = expRO(c_exp, b_exp, alpha, beta_exp, xi_mean, xi_cov, n1, N, num_exp, sim_num);
result_RO2 = expRO(c_exp, b_exp, alpha, beta_exp, xi_mean, xi_cov, 250, 500, num_exp, sim_num);

%RO_reconstruction
result_RO_reconstruct = expRO_reconstruct(c_exp, b_exp, alpha, beta_exp, xi_mean, xi_cov, n1, N, num_exp, sim_num);
result_RO_reconstruct2 = expRO_reconstruct(c_exp, b_exp, alpha, beta_exp, xi_mean, xi_cov, 250, 500, num_exp, sim_num);

%SG
result_SG = expSG(c_exp, b_exp, alpha, beta_exp, xi_mean, xi_cov, n1, N, num_exp, sim_num);
result_SG2 = expSG(c_exp, b_exp, alpha, beta_exp, xi_mean, xi_cov, 250, 500, num_exp, sim_num);%new on Apr.22

%SG_fast
result_SG_fast = expSG_fast(c_exp, b_exp, alpha, beta_exp, xi_mean, xi_cov, n1, N, num_exp, sim_num);
result_SG_fast2 = expSG_fast(c_exp, b_exp, alpha, beta_exp, xi_mean, xi_cov, 250, 500, num_exp, sim_num);

save('dtdriven', 'result_RO', 'result_RO_reconstruct', 'result_SG', 'result_DRO', 'result_SG_fast');
load('dtdriven');