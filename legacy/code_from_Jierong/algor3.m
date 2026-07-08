function x = algor3(delta, c, b, alpha,n, M, d,ksi,beta,n1,n2,K)
ksi_shuffle = ksi(randperm(size(ksi,1))',:);
ksi_train = ksi_shuffle(1:n1,:);
ksi_v = ksi_shuffle(n1 + 1:n,:);
m = length(delta);
l = n2 / K;
V = zeros(K,m);
C = zeros(1,m);
xx = zeros(d,m);
for i = 1:m
    x_ = Optimize(delta(i), c, b, alpha,n1, M, d,ksi_train);
    xx(:,i) = x_;
 
    for j = 1:K
        ksi_validation = ksi_v((j - 1) * l + 1:j * l,:);
        P = ksi_validation * x_ - b;
        cnt = 0;
        for k = 1:length(P)
            if P(k,1) <= 0
                cnt = cnt + 1;
            end
        end
        p = cnt / l;
        V(j,i) = p;
    end
    C(1,i) = c' * x_;
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
% delta_star = delta(min_ind);
% x = Optimize(delta_star, c, b, alpha,n, M, d,ksi);
x = xx(:,min_ind);
end