function [x_star1, x_star2, x_star3, x_star4] = phase2_1(data2, sol_star, b, c, alpha, beta, sim_num)
%phase2_1: contains 4 methods. phase one should be accomplished in PCxx.m
% 
% data2: data used to feed into phase two, n2-by-n, n is the dimension
% sol_star: solutions obtained from phase one n-by-m
% note: adjustment to phase2_2 as to avoid to many NaNs

[n2, n] = size(data2);

obj_star = data2 * sol_star; %n2-by-m
indifunc = obj_star <= b;
P_hat = mean(indifunc, 1); % 1-by-m
sigma_hat_sqr = P_hat .* (1 - P_hat); %1-by-m
sigma_hat = sqrt(sigma_hat_sqr); %1-by-m

% method one
ind1 = (P_hat >= (1 - alpha + norminv(1-beta, 0, 1) * sigma_hat / sqrt(n2)));
if sum(ind1) >0
    sol_feasible1 = sol_star(:, ind1);
    [~, index1] = min(c' * sol_feasible1);
    x_star1 = sol_feasible1(:, index1);
else
    x_star1 = NaN(n, 1);
end

% method two and three
Sigma_hat = cov(indifunc); % m-by-m
Z = mvnrnd(zeros(size(sol_star, 2), 1), Sigma_hat, sim_num);%m2
Z2 = mvnrnd(zeros(size(sol_star, 2), 1), Sigma_hat, sim_num);%m3
ind2 = sigma_hat > 0;
if sum(ind2) > 0
    Z = Z(:, ind2);
    Z2 = Z2(:, ind2);
    tmat = diag(1./sigma_hat(ind2));
    q = max(Z * tmat, [], 2); %m2, by row
    q2 = max(Z2, [], 2);%m3
    q_hat = prctile(q, (1-beta) * 100); 
    qu_hat = prctile(q2, (1-beta) * 100); 
else
    q_hat = 0;
    qu_hat = 0;
end

%m2
ind3 = (P_hat >= (1 - alpha + q_hat * sigma_hat / sqrt(n2)));
if sum(ind3) > 0
    sol_feasible2 = sol_star(:, ind3);
    [~, index2] = min(c' * sol_feasible2);
    x_star2 = sol_feasible2(:, index2);
else
    x_star2 = NaN(n, 1);
end

%m3
ind4 = (P_hat >= (1 - alpha + qu_hat / sqrt(n2)));
if sum(ind4) > 0
    sol_feasible3 = sol_star(:, ind4);
    [~, index3] = min(c' * sol_feasible3);
    x_star3 = sol_feasible3(:, index3);
else
    x_star3 = NaN(n, 1);
end

% method 4: a benchmark
ind5 = (P_hat >= (1 - alpha));
if sum(ind5) > 0
    sol_feasible4 = sol_star(:, ind5);
    [~, index4] = min(c' * sol_feasible4);
    x_star4 = sol_feasible4(:, index4);
else
    x_star4 = NaN(n, 1);
end 
end