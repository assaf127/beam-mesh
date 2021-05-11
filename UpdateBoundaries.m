function [boundaries_mat, cells_in_triangles] = UpdateBoundaries( ...
    cells_in_triangles, boundaries_mat, triangle_centroids, metrics, ...
    cell_centers, F_adj)
%UpdateBoundaries Updates the Voronoi cells
%   Updates the Voronoi cells in the previous boundaries, and produces the
%   new boundaries to be updated in the next iteration.
% Input:
%   cells_in_triangles  - F x 1, contains the Voronoi cell of each triangle
%   boundaries_mat      - F x N, contains the cells that have to be checked
%                                for every triangle on the boundaries
%   triangle_centroids  - F x 3, the (x,y,z) locations of each triangle
%   metrics             - N x 9, every line is the metric matrix (after
%                                reshaping to 1 x 9) of the cell
%   cell_centers        - N x 3, the (x,y,z) locations of cell centroids
%   F_adj               - F x F, the faces adjacency matrix
% Output:
%   boundaries_mat      - F x N, contains the *updated* cells that will
%                                have to be checked for every triangle on
%                                the boundaries
%   cells_in_triangles  - F x 1, contains the *updated* Voronoi cell of
%                                each triangle
N = size(boundaries_mat, 2);
F = size(boundaries_mat, 1);
dist_mat = double(boundaries_mat);  % preallocate

for i=1:N
    indices = boundaries_mat(:, i) > 0;
    dist_mat(indices, i) = ...
        1 ./FindAnisotropicDistance(triangle_centroids(indices, :) - ...
                            cell_centers(i, :), reshape(metrics(i, :),3,3));
end

[update_rows, ~] = find(dist_mat);
rows = unique(update_rows);
[~, rowmax] = max(dist_mat(rows, :),[], 2);

cells_in_triangles(rows) = rowmax;
[triangles,~,cells]  = find(cells_in_triangles);
triangles_in_cells_mat = sparse(triangles, cells, ones(size(triangles)), F, N);
boundaries_mat = GetBoundaryTrianglesMatrix(triangles_in_cells_mat, F_adj);
end

% delta: Nx3 matrix, metric: 3x3 matrix
function [dist] = FindAnisotropicDistance(delta, metric)
    dist = sum((delta * metric) .* delta, 2);
end



