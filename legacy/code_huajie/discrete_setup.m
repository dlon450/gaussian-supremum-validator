close all;
clc;

d = 10;
c = -ones(d, 1);
b = 8.33;
alpha = 0.05;
beta = 0.05;


n = 500;
data = abs(randn(n, d));


delta1 = 6./1.8.^(1:10); % large delta values
delta2 = delta1(end)./1.3.^(1:10); % small delta values



delta = [delta1, delta2];


kappa = 5;
tic;
x = py.knapsack.DRO_solver(delta, c.', b, alpha, data(:).', kappa);
x = numeric(x);
x = reshape(x, d, []);
cost = c.'*x;
test_data = abs(randn(1e6, d));
constraint = mean(test_data*x <= b, 1);
toc;
