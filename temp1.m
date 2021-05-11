clear;
M = Mesh3D('vase.off');
N = 200;

F_Curv = abs(M.IV2F * M.MeanCurvature);
[density,I] = sort(F_Curv);
int_density = cumsum(density);
rand_results = randsrc(N, 1, [I, density / sum(density)]');

temp = zeros(M.NF, 1);
temp(rand_results) = 1;
while( sum(temp) ~= N ) % Make sure we have N unique centers
    rand_results = randsrc(N - sum(temp), 1, [I, density / sum(density)]');
    temp(rand_results) = 1;
end
triangle_centroids = (M.V(M.F(:,1),:) + M.V(M.F(:,2),:) + M.V(M.F(:,3),:)) / 3;
centers = triangle_centroids(logical(temp), :);

prev_cells = FindIsotropicVoronoiCells(M, centers);
[centers, metrics] = FindCentroidAndMetric(M, prev_cells, N);
next_cells = FindVoronoiCells(M, centers, metrics, prev_cells);

while( sum(prev_cells ~= next_cells) > 30 )
    prev_cells = next_cells;
    [centers, metrics] = FindCentroidAndMetric(M, prev_cells, N);
    next_cells = FindVoronoiCells(M, centers, metrics, prev_cells);
end

figure;
M.ShowMesh;
M.PlotFunction(next_cells, false);
hold on;
scatter3(centers(:,1), centers(:,2), centers(:,3),20,'filled','r');

[U, S, ~] = svd(metrics{1});
[x,y,z] = ellipsoid(centers(1,1), centers(1,2), centers(1,3), S(1,1), S(2,2), S(3,3));
r = (U * [x(:) y(:) z(:)]')';
X = reshape(r(:,1),21,21);
Y = reshape(r(:,2),21,21);
Z = reshape(r(:,3),21,21);
surf(X, Y, Z);


