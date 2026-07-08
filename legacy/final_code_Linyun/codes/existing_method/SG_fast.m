function [solution] = SG_fast(c,b,data,N1)
% c: cost coefficient
% b: threshold
% data: whole data set, each row is an observation of xi
% N1: sample size for phase one
n = length(c);    % dimension
x_fast0=zeros(n,1);


data1 = data(1:N1,:);
data2 = data(N1+1:end,:);


% phase one
cvx_begin quiet
    variable x_fast(n)
    minimize(c.'*x_fast)
    subject to
        data1*x_fast <= b
cvx_end


% phase two
if strcmp(cvx_status,'Solved')
    cvx_begin quiet
    variable alpha_fast
    minimize(c.'*((1-alpha_fast)*x_fast + alpha_fast*x_fast0))
    subject to
    data2*((1-alpha_fast)*x_fast + alpha_fast*x_fast0) <= b
    alpha_fast >= 0
    alpha_fast <= 1
    cvx_end
    
    solution = (1 - alpha_fast)*x_fast + alpha_fast*x_fast0;
else
    solution = NaN(n,1);
end
end