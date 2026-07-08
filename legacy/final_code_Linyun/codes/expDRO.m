function [out] = expDRO(c, b, alpha, beta, mu, sigma, n1, N, num_exp, sim_num)
%expDRO: experiment on RO method
% 
% c: cost coefficient
% b: threshold
% alpha: safety level
% beta: confidence level
% mu: true mean of xi
% sigma: true covariance matrix of xi
% n1: num of observations to enter phase one
% N: total num of observations

n = length(c); % dimension
sol_DRO = zeros(n, num_exp);
sol_PCDRO1 = zeros(n, num_exp);
sol_PCDRO2 = zeros(n, num_exp);
sol_PCDRO3 = zeros(n, num_exp);
sol_PCDRO4 = zeros(n, num_exp);
sol_phase1 = cell(1, num_exp);
rho_est = sqrt(chi2inv(1-beta, (n+n*(n+1)/2))/n1);

gendata = mvnrnd(mu, sigma, N * num_exp);
file_name = ['DRO_data', date()];
testdata = mvnrnd(mu, sigma, 10000);


for i = 1:num_exp
    data = gendata(((i-1)*N + 1):(i*N), :); %(n1+n2)-by-d
    %sol_DRO(:, i) = DRO2(c, b, alpha, beta, data);
    
    %phase one
    sol_star = DRO2_path(c, b, alpha, data(1:n1, :), linspace(0,(rho_est * 1.5),50));
    sol_star = sol_star(:, ~isnan(sol_star(1,:))); %sol_star is d-by-m
    sol_phase1{1, i} = sol_star;
    
    %phase two
    [sol_PCDRO1(:, i), sol_PCDRO2(:, i), sol_PCDRO3(:, i), sol_PCDRO4(:, i)] = phase2_1(data((n1+1):N, :),...
        sol_star, b, c, alpha, beta, sim_num);
    
    disp(i);
end

save(file_name, 'gendata', 'testdata', 'sol_phase1', 'rho_est');%save phase one solutions and the material data

%output: including saving the solutions
out = output(b, c, sigma, alpha, sol_DRO, sol_PCDRO1, sol_PCDRO2, sol_PCDRO3, sol_PCDRO4, testdata, 'DRO');
end


