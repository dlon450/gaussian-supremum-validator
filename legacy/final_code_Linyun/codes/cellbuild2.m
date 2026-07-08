function [out] = cellbuild2(b, c, sigma, alpha, beta, sim_num, name1, name2, name3)
% transfer from 3 methods to 4
% name1: raw data
% name2: outcomes
% name3: for output

load(name1);
load(name2);

[~, num_exp] = size(sol_phase1);
[k, n] = size(gendata);
N = k / num_exp;
n1 = N / 2;

sol_1 = zeros(n, num_exp);
sol_2 = zeros(n, num_exp);
sol_3 = zeros(n, num_exp);
sol_4 = zeros(n, num_exp);

for i = 1:num_exp
    data = gendata(((i-1)*N + 1):(i*N), :); %(n1+n2)-by-d
    [sol_1(:, i), sol_2(:, i), sol_3(:, i), sol_4(:, i)] = phase2_2(data((n1+1):N, :),...
        sol_phase1{1, i}, b, c, alpha, beta, sim_num);%avoid too many NaNs
end

out = output(b, c, sigma, alpha, sol_ori_omit, sol_1, sol_2, sol_3, sol_4, testdata, name3);
end