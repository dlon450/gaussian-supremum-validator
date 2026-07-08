% Uncomment the following two lines if changes are made to transport.py
% clear classes;
% py.importlib.reload(py.importlib.import_module('transport'));


n = 1000;
beta = 0.05;
tic;
max_delta = get_max_delta(n, beta, @rng_ksi);
toc;