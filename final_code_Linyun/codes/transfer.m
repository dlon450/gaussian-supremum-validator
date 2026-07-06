%% add one extra method and adjustments
cd F:\GoogleDrive\OhWork\Prof.Lam\NewFrameWorkforChanceConstraintProgram\simulations\codes

%% setting 1
addpath(genpath('F:\GoogleDrive\OhWork\Prof.Lam\NewFrameWorkforChanceConstraintProgram\reports\version4_May-08-2019\exp_results\setting1'))
load('setting');

% N = 200
%DRO
result_DRO = cellbuild(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num, 'DRO_data_09-Mar-2019', 'DRO_exp_09-Mar-2019', 'DRO');
%RO
result_RO = cellbuild(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num, 'RO_data_07-Mar-2019', 'RO_exp_08-Mar-2019', 'RO');
%RO_reconstruction
result_RO_reconstruct = cellbuild(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num, 'RO_rec_data_04-Mar-2019', 'RO_rec_exp_08-Mar-2019','RO_rec');
%SG
result_SG = cellbuild2(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num, 'SG_data_05-Mar-2019', 'SG_exp_08-Mar-2019', 'SG');
%extra
result_SG2 = cellbuild2(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num, 'SG_data_15-May-2019', 'SG_exp_16-May-2019', 'SG');
%SG_fast
result_SG_fast = cellbuild(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num, 'SG_fast_data_04-Mar-2019', 'SG_fast_exp_07-Mar-2019','SG_fast');

save('dtdriven', 'result_RO', 'result_RO_reconstruct', 'result_SG', 'result_DRO', 'result_SG_fast');

% N = 500
%DRO
result_DRO = cellbuild(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num, 'DRO_data26-Apr-2019', 'DRO_exp_26-Apr-2019', 'DRO');
%RO
result_RO = cellbuild(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num, 'RO_data_26-Apr-2019', 'RO_exp_26-Apr-2019', 'RO');
%RO_reconstruction
result_RO_reconstruct = cellbuild(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num, 'RO_rec_data_26-Apr-2019', 'RO_rec_exp_27-Apr-2019','RO_rec');
%SG
result_SG = cellbuild2(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num, 'SG_data_23-Apr-2019', 'SG_exp_24-Apr-2019', 'SG');
%SG_fast
result_SG_fast = cellbuild(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num, 'SG_fast_data_27-Apr-2019', 'SG_fast_exp_27-Apr-2019','SG_fast');

save('dtdriven', 'result_RO', 'result_RO_reconstruct', 'result_SG', 'result_DRO', 'result_SG_fast');

%% setting 2
addpath(genpath('F:\GoogleDrive\OhWork\Prof.Lam\NewFrameWorkforChanceConstraintProgram\reports\version4_May-08-2019\exp_results\setting2'))
load('setting');

% N = 200
%DRO
result_DRO = cellbuild(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num, 'DRO_data_06-Mar-2019', 'DRO_exp_09-Mar-2019', 'DRO');
%RO
result_RO = cellbuild(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num, 'RO_data_08-Mar-2019', 'RO_exp_08-Mar-2019', 'RO');
%RO_reconstruction
result_RO_reconstruct = cellbuild(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num, 'RO_rec_data_05-Mar-2019', 'RO_rec_exp_08-Mar-2019','RO_rec');
%SG
result_SG = cellbuild2(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num, 'SG_data_05-Mar-2019', 'SG_exp_09-Mar-2019', 'SG');
%extra
result_SG2 = cellbuild2(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num, 'SG_data_15-May-2019', 'SG_exp_16-May-2019', 'SG');
%SG_fast
result_SG_fast = cellbuild(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num, 'SG_fast_data_05-Mar-2019', 'SG_fast_exp_07-Mar-2019','SG_fast');

save('dtdriven', 'result_RO', 'result_RO_reconstruct', 'result_SG', 'result_DRO', 'result_SG_fast');

% N = 500
%DRO
result_DRO = cellbuild(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num, 'DRO_data26-Apr-2019', 'DRO_exp_26-Apr-2019', 'DRO');
%RO
result_RO = cellbuild(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num, 'RO_data_26-Apr-2019', 'RO_exp_26-Apr-2019', 'RO');
%RO_reconstruction
result_RO_reconstruct = cellbuild(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num, 'RO_rec_data_26-Apr-2019', 'RO_rec_exp_27-Apr-2019','RO_rec');
%SG
result_SG = cellbuild2(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num, 'SG_data_23-Apr-2019', 'SG_exp_25-Apr-2019', 'SG');
%SG_fast
result_SG_fast = cellbuild(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num, 'SG_fast_data_27-Apr-2019', 'SG_fast_exp_27-Apr-2019','SG_fast');

save('dtdriven', 'result_RO', 'result_RO_reconstruct', 'result_SG', 'result_DRO', 'result_SG_fast');

%% setting 3
%DRO
result_DRO = 0;
%RO
result_RO = cellbuild(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num, 'RO_data_08-May-2019', 'RO_exp_08-May-2019', 'RO');
%RO_reconstruction
result_RO_reconstruct = cellbuild(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num, 'RO_rec_data_08-May-2019', 'RO_rec_exp_08-May-2019','RO_rec');
%SG
result_SG = cellbuild2(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num,'SG_data_10-May-2019', 'SG_exp_11-May-2019', 'SG');
%SG_fast
result_SG_fast = cellbuild(b_exp, c_exp, xi_cov, alpha, beta_exp, sim_num, 'SG_fast_data_08-May-2019', 'SG_fast_exp_08-May-2019','SG_fast');

save('dtdriven', 'result_RO', 'result_RO_reconstruct', 'result_SG', 'result_DRO', 'result_SG_fast');
load('dtdriven');