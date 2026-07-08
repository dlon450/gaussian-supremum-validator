%test the selection of range in S (discretized uncertainty set)
load('setting');
b = b_exp;
c = c_exp;
beta = beta_exp;

i = 1;
data = gendata(((i-1)*N + 1):(i*N), :); %(n1+n2)-by-d
[~, rho_est] = DRO2(c, b, alpha, beta, data);
    
[~, rho_est] = DRO2(c, b, alpha, beta, data(1:n1, :));

sol_star = DRO2_path(c, b, alpha, data(1:n1, :), linspace(0,(rho_est + 20),50));
sol_star = DRO2_path(c, b, alpha, data(1:n1, :), linspace(0,(rho_est*1.5),50));
sol_star = DRO2_path(c, b, alpha, data(1:n1, :), linspace(0,rho_est,50));
sol_star = DRO2_path(c, b, alpha, data(1:n1, :), linspace(0,0.8 * rho_est,50));
sol_star = DRO2_path(c, b, alpha, data(1:n1, :), linspace(0,0.5 * rho_est,50));

sol_star = sol_star(:, ~isnan(sol_star(1,:)));
data2 = data((n1+1):N, :);

[n2, n] = size(data2);

obj_star = data2 * sol_star; %n2-by-m
indifunc = obj_star <= b;
P_hat = mean(indifunc, 1); % 1-by-m
sum(P_hat)
sigma_hat_sqr = P_hat .* (1 - P_hat); %1-by-m
sigma_hat = sqrt(sigma_hat_sqr); 

%all results the same: P_hat = ones(50, 1)