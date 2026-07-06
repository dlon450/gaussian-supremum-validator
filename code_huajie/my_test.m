% clear all;
close all;
clc;


% clear classes;
% py.importlib.reload(py.importlib.import_module('knapsack'));

d = 10;
% c = -rand(d, 1);
c = -ones(d, 1);
% c = [-7;-6;-7;-2;-9;-4;-6;-10;-4;-2];
b = 8.33;
alpha = 0.05;
% beta = 0.05;


n = 1000;
data = rng_data(n, d);
% delta = [0.001, 0.005, 0.01, 0.05, 0.1];
% delta = 0.001;
delta1 = 6./2.^(1:10);
delta2 = delta1(end)./1.3.^(1:15);
delta = [delta1, delta2];
test_data = rng_data(1e6, d);


% cvx_solver gurobi_2
% tic;
% x = CCP_DRO_Wasserstein(delta, c, b, alpha, data);
% constraint1 = mean(test_data*x <= b, 1);
% toc;

kappa = 5;
tic;
x = py.knapsack.DRO_solver(delta, c.', b, alpha, data(:).', kappa);
x = numeric(x);
x = reshape(x, d, []);
cost = c.'*x;
constraint = mean(test_data*x <= b, 1);
toc;


% tic;
% x1 = CCP_DRO_WassersteinGRB(delta, c, b, alpha, data);
% toc;

% tic;
% x = CCP_DRO_Wasserstein(delta, c, b, alpha, data);
% toc;

% test_size = 1e5;
% data_test = rng_data(test_size, d);
% feas_test = mean(data_test * x <= b, 1);
% feas_test2 = mean(data_test * x <= b, 1);


% tic;
% x = DRO_moment(delta, c, b, alpha, ksi);
% cost = c.'*x;
% toc;

% tic;
% x = RO_ellipsoid(delta, c, b, ksi);
% cost = c.'*x;
% toc;

% tic;
% x = SO_FAST(delta, c, b, ksi);
% cost = c.'*x;
% toc;

% tic;
% x = SO(1:n, c, b, ksi);
% cost = c.'*x;
% toc;

% tic;
% x = RO_LCX(delta, c, b, alpha, ksi);
% cost = c.'*x;
% toc;

% cvx_solver mosek
% tic;
% x = CCP_SAA(0.95, c, b, data);
% cost = c.'*x;
% toc;

% cvx_solver mosek
% tic;
% x2 = new_optimize(delta, c, b, alpha,n,d,ksi);
% toc;

% cvx_solver gurobi_2
% tic;
% x3 = new_optimize(delta, c, b, alpha,n,d,ksi);
% toc;

% cvx_solver gurobi_2
% tic;
% [x4, x5, x6] = GS(delta, c, b, alpha,n,d,ksi,beta,n1,n2);
% toc;