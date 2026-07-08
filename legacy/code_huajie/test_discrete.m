% clear all;
clc;

% clear classes;
% py.importlib.reload(py.importlib.import_module('knapsack'));


d = 20;
c = -rand(d, 1);
b = 0.8*d;
alpha = 0.05;
% beta = 0.05;
n = 1000;
data = abs(randn(n, d));
test_data = abs(randn(1e6, d));

% delta = 6./1.5.^(1:30);
delta = 0.01;


% Use the default kappa=1, i.e., in each round, one violated constraint is added to the model
tic;
x1 = py.knapsack.DRO_solver(delta, c.', b, alpha, data(:).');
x1 = numeric(x1);
x1 = reshape(x1, d, []);
cost1 = c.'*x1;
constraint1 = mean(test_data*x1 <= b, 1);
toc;

% Use a different kapppa
kappa = 5;
tic;
x2 = py.knapsack.DRO_solver(delta, c.', b, alpha, data(:).', kappa);
x2 = numeric(x2);
x2 = reshape(x2, d, []);
cost2 = c.'*x2;
constraint2 = mean(test_data*x2 <= b, 1);
toc;