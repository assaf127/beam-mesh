function [metrics, cell_centers_positions, cell_center_triangles] = ...
    UpdateCentersAndMetrics( ...
    STLFormat, triangle_centroids, cells_in_triangles, num_centers, triangles_area, ...
    cell_center_triangles)
%UpdateCentersAndMetrics Updates the Voronoi generators and metrics
%   Computes the new Voronoi cell centroids and metrics to be used in the
%   next iteration
% Input:
%   STLFormat               - F x 9, every row contains the positions of
%                             a triangle vertices in the following format:
%                             [v1x v2x v3x v1y v2y v3y v1z v2z v3z]
%   triangle_centroids      - F x 3, contains the (x,y,z) locations of each
%                             triangle centroid
%   cells_in_triangles      - F x 1, the Voronoi cell of each triangle
%   num_centers             - The number of Voronoi cells
%   triangles_area          - F x 1, the area of the triangles
%   cell_center_triangles   - N x 1, the index of the triangle in the
%                             center of each Voronoi cell
% Output:
%   metrics                 - N x 9, every line is the metric matrix (after
%                             reshaping to 1 x 9) of the cell
%   cell_centers_positions  - N x 3, contains the *updated* centroids of
%                             each Voronoi cell
%   cell_center_triangles   - N x 1, the *updated* index of the triangle in
%                             the center of each Voronoi cell
num_triangles = size(cells_in_triangles,1);

[cell_triangles,~,triangle_cells]  = find(cells_in_triangles);
cell_triangles_centroids = triangle_centroids(cell_triangles, :);
cols_arr = repmat([1 2 3], size(cell_triangles));

cell_centers_positions = accumarray([repmat(triangle_cells,3,1), cols_arr(:)], ...
            cell_triangles_centroids(:), [num_centers 3], @mean);

metrics = zeros(num_centers, 9);        % preallocate
covariance_template = [2 1 1;
                       1 2 1;
                       1 1 2] / 12;

[triangles,~,cells]  = find(cells_in_triangles);
dist_mat = sparse(triangles, cells, ones(size(triangles)), num_triangles, num_centers);  % preallocate
for i=1:num_centers     % parfor
    temp_metric = zeros(3);
    temp_centroid = cell_centers_positions(i, :); % might be different
    curr_triangles = find(dist_mat(:, i));
    if(isempty(curr_triangles))
        continue;
    end
    [~, min_ind] = min(FindIsotropicDistance(triangle_centroids(curr_triangles, :) - ...
                            temp_centroid ));

    cell_center_triangles(i) = curr_triangles(min_ind);

    cell_vertices_mats = reshape(STLFormat(curr_triangles,:), 3*size(curr_triangles,1), 3) - temp_centroid;
    cell_vertices_mats = reshape(cell_vertices_mats', 3*size(curr_triangles,1), 3)';
    cell_triangles_area = triangles_area(curr_triangles);
    
    for t=1:size(curr_triangles,1)       % t=curr_triangles'
        %vertices = reshape(STLFormat(t,:),3,3) - temp_centroid;
        %temp_metric = temp_metric + triangles_area(t) * ...
        %                vertices' * covariance_template * vertices;
        vertices = cell_vertices_mats(:,(3*t-2):(3*t));
        temp_metric = temp_metric + cell_triangles_area(t) * ...
                        vertices' * covariance_template * vertices;
    end
    [U, S, V] = svd(temp_metric);
    S(1,1) = S(3,3) / S(1,1);
    S(2,2) = S(3,3) / S(2,2);
    S(3,3) = 1;
    metrics(i, :) = reshape(U * S * V', 1, 9);
end

cell_centers_positions = triangle_centroids(cell_center_triangles, :);
end

function [dist] = FindIsotropicDistance(delta)
    dist = dot(delta, delta, 2);
end
