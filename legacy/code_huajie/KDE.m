function density = KDE(data, pts)
% data: the data matrix, each row is an observation
% pts: matrix of points where the density is to be evaluated

bw = bw_KDE(data); % calculate bandwidth
m = size(pts, 1);
density = zeros(m, 1);

for k = 1:m
    density(k) = mean(mvnpdf(pts(k, :), data, bw.^2));
end

end