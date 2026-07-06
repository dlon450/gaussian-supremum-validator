function sample = sample_KDE(data, m)
% data: the data matrix, each row is an observation
% m: sample size

bw = bw_KDE(data); % calculate bandwidth
n = size(data, 1);
mu = data(randsample(n, m, true), :);
sample = mvnrnd(mu, bw.^2);

end