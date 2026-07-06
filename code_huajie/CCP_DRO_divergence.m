function [solution] = CCP_DRO_divergence(para, c, b, alpha, data)
% para: a sequence of parameter values specifying the KL distance
% c: cost vector in the objective
% b: right-hand threshold in the constraint
% alpha: tolerance level in the original chance constraint
% data: data matrix, each row is an observation


%% calculate the perturbed tolerance level
m = length(para);
tolerance = zeros(m, 1);
precision = 1e-6;
for k = 1:m
    x_l=0;
    x_r=1;
    while (x_r-x_l)/x_r > precision
        x_m = (x_l + x_r)/2;
        if 1 - alpha*exp(-para(k))*x_m^(1-alpha) - (1-alpha)*exp(-para(k))*x_m^(-alpha) > 0
            x_r = x_m;
        else
            x_l = x_m;
        end
    end
    
    tolerance(k) = 1 - (exp(-para(k))*((x_l+x_u)/2)^(1-alpha) - 1) / (x-1);
    tolerance(k) = 1 - max(tolerance(k), 0);
end


%% solve the perturbed chance constrained problem using sample average approximation
SAA_size = 500;
SAA_data = sample_KDE(data, SAA_size);
solution = CCP_SAA(tolerance, c, b, SAA_data);


end