function x = algor2(delta, c, b, alpha,n, M, d,ksi,beta,B)
R = randi(n,B,n);
set = zeros(B,n);
set_comp = zeros(B,n);
for i = 1:B
    for j = 1:n
        if ismember(j,R(i,:))
            set(i,j) = 1;
        else
            set_comp(i,j) = 1;
        end
    end
end
m = length(delta);
denom = sum(set_comp,2);
V = zeros(B,m);
C = zeros(1,m);
for i = 1:m
    cost = 0;
    for j = 1:B
        ksi_validation = ksi(set_comp(j,:) == 1,:);
        ksi_train = ksi(R(j,:),:);
        x_ = Optimize(delta(i), c, b, alpha,n , M, d,ksi_train);
        P = ksi_validation * x_ - b;
        cnt = 0;
        cost = cost + c' * x_;
        for k = 1:length(P)
            if P(k,1) <= 0
                cnt = cnt + 1;
            end
        end
        p = cnt / denom(j);
        V(j,i) = p;
    end
    C(1,i) = cost / B;
end

qualified = zeros(1,m);
for i = 1:m
    cnt = 0;
    for j = 1:B
        if V(j,i) >= 1 - alpha
            cnt = cnt + 1;
        end
    end
    if cnt / B >= 1 - beta
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
            
    