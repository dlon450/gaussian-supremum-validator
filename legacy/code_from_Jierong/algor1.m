function x = algor1(delta, c, b, alpha,n, M, d,ksi,beta,K)
m = length(delta);
l = n / K;
V = zeros(K,m);
C = zeros(1,m);
ksi_shuffle = ksi(randperm(size(ksi,1))',:);
for i = 1:m
    cost = 0;
    for j = 1:K
        ksi_validation = ksi_shuffle((j - 1) * l + 1:j * l,:);
        ksi_train = ksi_shuffle([1:(j - 1) * l,j * l + 1 : n],:);
        x_ = Optimize(delta(i), c, b, alpha,n - l, M, d,ksi_train);
        P = ksi_validation * x_ - b;
        cnt = 0;
        cost = cost + c' * x_;
        for k = 1:length(P)
            if P(k,1) <= 0
                cnt = cnt + 1;
            end
        end
        p = cnt / l;
        V(j,i) = p;
    end
    C(1,i) = cost / K;
end
qualified = zeros(1,m);
for i = 1:m
    cnt = 0;
    for j = 1: K
        if V(j,i) >= 1 - alpha
            cnt = cnt + 1;
        end
    end
    if cnt / K >= 1 - beta
        qualified(1,i) = 1;
    end
end
min_val = inf;
for i = 1:m
    if qualified(1,i) == 1
        if C(1,i) < min_val
            min_val = C(1,i);
            min_ind = i;
        end
    end
end
delta_star = delta(min_ind);
x = Optimize(delta_star, c, b, alpha,n, M, d,ksi);

end