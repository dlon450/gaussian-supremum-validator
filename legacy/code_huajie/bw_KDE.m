function bw = bw_KDE(data)
% data: the data matrix, each row is an observation

[n, d] = size(data);
bw = std(data)/n^(1/(d+4)); % Using Scott's rule to set bandwidth
end

