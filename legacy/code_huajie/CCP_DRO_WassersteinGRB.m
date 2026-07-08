function solution = CCP_DRO_WassersteinGRB(para, c, b, alpha, data)
% para: a sequence of parameter values specifying the KL distance
% c: cost vector in the objective
% b: right-hand threshold in the constraint
% alpha: tolerance level in the original chance constraint
% data: data matrix, each row is an observation


[n, d] = size(data);
M = ceil(abs(b) + max(sum(abs(data),2)));
m = length(para);
solution = zeros(d, m);


% A_f_dro_1=reshape(Samples_for_SAA',n,N_data_solve*m)';
%     b_f_dro_1=repmat(b1,N_data_solve,1);
%
%     A_f_dro=[A_f_dro_1];
%     b_f_dro=[b_f_dro_1];

%objective
%     model_obj=[c;zeros(N_data_solve,1)];
d_opt = d+3*n+2;
obj = [c;zeros(d_opt-d, 1)];
A1 = [sparse(n, d+n), -speye(n), -speye(n), sparse(1:n, 2, 1, n, 2)];
rhs1 = zeros(n, 1);
A2 = [sparse(1, d+n), ones(1, n), sparse(1, n+1), -n*alpha];
rhs2 = 0;
A3 = [data, M*speye(n), sparse(n, n), speye(n), sparse(n, 2)];
rhs3 = (b+M)*ones(n, 1);
A4 = [sparse(n, d), -M*speye(n), sparse(n, n), speye(n), sparse(n, 2)];
rhs4 = zeros(n, 1);
A5 = [speye(d, d), sparse(d, 3*n), sparse(1:d, 1, -1, d, 2)];
rhs5 = zeros(d, 1);

A = sparse([A1;A2;A3;A4;A5]);
rhs = [rhs1;rhs2;rhs3;rhs4;rhs5];
sense = repmat('<', 3*n+d+1, 1);

lb = zeros(d_opt, 1);
ub = [ones(d+n, 1); inf(d_opt-d-n, 1)];
type = [repmat('C',d,1);repmat('B',n,1);repmat('C',d_opt-d-n,1)];

%     % setup constraints
%     chance_lhs=[zeros(n,1);ones(N_data_solve,1)]';
%     chance_rhs=p;
%     chance_sns='<';

%     cut_region_lhs=zeros(N_data_solve,N_data_solve+n);
%     for kk=1:N_data_solve
%         for_z=zeros(1,N_data_solve);
%         for_z(kk)=-big_M;
%         cut_region_lhs(kk,:)=[A_f_dro(kk,:),for_z];
%     end
%     cut_region_rhs=repmat(b1,N_data_solve,1);
%     cut_region_sns=repmat('<',N_data_solve,1);
%
%     % X,Z deterministic
%     lb_dec=[repmat(-inf,n,1);zeros(N_data_solve,1)];
%     ub_dec=[repmat(inf,n,1);ones(N_data_solve,1)];


%     constraint_lhs=[chance_lhs;cut_region_lhs];
%     constraint_rhs=[chance_rhs;cut_region_rhs];
%     constraint_sense=[chance_sns;cut_region_sns];



% params.outputflag = 0;

% cvx_solver gurobi_2
for k = 1:m
    
    A(n+1, end-1) = n*para(k);
    
    clear model;
    model.A = A;
    model.obj = obj;
    model.modelsense = 'min';
    model.rhs = rhs;
    model.sense = sense;
    model.lb = lb;
    model.ub = ub;
    model.vtype = type;
    
%     results=gurobi(model, params);
    results=gurobi(model);
    
    
    %     cvx_begin quiet
    %
    %         variable x(d)
    %         variable y(n) binary
    %         variable v nonnegative
    %         variable r nonnegative
    %         variable z(n) nonnegative
    %         variable s(n) nonnegative
    %
    %         minimize (c.'*x)
    %
    %         subject to
    %         z + s >= r
    %         para(k) * v + sum(z) / n <= alpha * r
    %         b - data * x + M * (1 - y) >= s
    %         M * y >= s
    %         v >= x
    %         0 <= x <= 1
    %
    %     cvx_end
    
    solution(:, k) = results.x(1:d);
end


end

