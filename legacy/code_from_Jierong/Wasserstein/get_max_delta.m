function delta = get_max_delta(n, beta, rng_ksi)
% n: sample size
% beta: 1 - beta is the confidence level
% rng_ksi: a function that can generate a sample of ksi of a given size

B = 200;
dist_set = zeros(B, 1);
for b = 1:B
    ksi_1 = rng_ksi(n);
    ksi_2 = rng_ksi(n);

    cost = zeros(n, n);
    for i = 1:n
        for j = 1:n
            cost(i,j) = norm(ksi_1(i,:) - ksi_2(j,:), 1);
        end
    end
    cost = cost(:).';
    precision = 0.01;
    
    distance = py.transport.distance(cost, n, precision);
    if distance == 'unsolved'
        error('Wasserstein distance calculation failed.')
    end
    dist_set(b) = distance;
    
    
%     A1 = sparse(n, n^2);
%     for i = 1:n
%         A1(i, (i-1)*n + 1:i*n) = 1;
%     end 
%     A2 = repmat(sparse(eye(n)), 1, n);
%     
%     Aeq = [A1;A2];
%     beq = ones(2*n, 1);
%     LB = zeros(n^2, 1);
%     UB = ones(n^2, 1);
    
%     tic;
%     [~, fval] = linprog(f, [], [], Aeq, beq, LB);
%     [~, fval] = intlinprog(f, 1:n^2, [], [], Aeq, beq, LB, UB);
%     dist_lin(b) = fval/n;
%     toc;
end

dist_set = sort(dist_set);
delta = dist_set(ceil(B*(1-beta)));
end