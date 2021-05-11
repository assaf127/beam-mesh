function [new_V, new_F, mesh_vertices, edges] = ...
    CoarseRemeshMahalanobisCVD(M, N, showProgress)
%CoarseRemeshMahalanobisCVD Computes a coarse dual triangular mesh
%   Computes a coarse dual triangular mesh with N faces for mesh M using
%   Lloyd's algorithm for finding CVD with the Mahalanobis metric.
% Input:
%   M           - Mesh3D object
%   N           - The number of generators (and therefor faces) for the CVD
%   showProgress- Pass true to view the algorithm progression in situ
% Output:
%   new_V   - V x 3 matrix with the (x,y,z) locations of the vertices of
%             the resulting coarse mesh
%   new_F   - F x max_polygon_degree matrix with the faces of the resulting
%             mesh. Unused areas are filled with NaNs
%   mesh_vertices   - The indices of the resulting vertices in the original
%                     mesh
%   edges   - E x 2 matrix where each row contains 2 vertex indices
%             representing an edge in the coarse mesh
if(~exist('M', 'var'))
    M = Mesh3D('cat_s3.off');
end
if(~exist('N', 'var'))
    N = 100;
end
if(~exist('showProgress', 'var'))
    showProgress = true;
end

F_Curv = abs(M.IV2F * M.MeanCurvature);
[density,I] = sort(F_Curv);
%int_density = cumsum(density);
rand_results = randsrc(N, 1, [I, density / sum(density)]');

temp = zeros(M.NF, 1);
temp(rand_results) = 1;
while( sum(temp) ~= N ) % Make sure we have N unique centers
    rand_results = randsrc(N - sum(temp), 1, [I, density / sum(density)]');
    temp(rand_results) = 1;
end
triangle_centroids = (M.V(M.F(:,1),:) + M.V(M.F(:,2),:) + M.V(M.F(:,3),:)) / 3;
prev_center_triangles = find(temp == 1);
temp(prev_center_triangles) = (1:N)';
prev_centers = triangle_centroids(prev_center_triangles, :);
prev_faces_in_cells = sparse(prev_center_triangles, (1:N)', ones(N,1), M.NF, N);
F_adj = M.F_adj + speye(M.NF);
F_Area = M.F_Area;
STLFormat = reshape(M.V(M.F,:), M.NF, 9);

prev_boundaries = GetBoundaryTrianglesMatrix(prev_faces_in_cells, F_adj);
prev_metrics = repmat([1 0 0 0 1 0 0 0 1], N, 1);


[next_boundaries, next_triangle_cells] = ...
    UpdateBoundaries( temp, prev_boundaries, ...
    triangle_centroids, prev_metrics, prev_centers, F_adj);

[next_metrics, next_centers, next_center_triangles] = UpdateCentersAndMetrics( ...
    STLFormat, triangle_centroids, next_triangle_cells, N, F_Area, prev_center_triangles);
prev_triangle_cells = temp;

iteration = 1;

if(showProgress)
    figure;
end
%% just initializations: sum(next_triangle_cells == 0) > 0
while( sum(next_triangle_cells == 0) > 0  )
    if(showProgress && mod(iteration, 11) == 0 )
        M.PlotRemeshing(next_triangle_cells, false);
        drawnow;
    end
    prev_triangle_cells = next_triangle_cells;
    prev_boundaries = next_boundaries;
    [next_boundaries, next_triangle_cells] = ...
            UpdateBoundaries( prev_triangle_cells, prev_boundaries, ...
                triangle_centroids, prev_metrics, prev_centers, F_adj);
    prev_prev_centers = prev_centers;
    prev_centers = next_centers;
    prev_metrics = next_metrics;
    prev_center_triangles = next_center_triangles;
    [next_metrics, next_centers, next_center_triangles] = UpdateCentersAndMetrics( ...
            STLFormat, triangle_centroids, next_triangle_cells, N, F_Area, prev_center_triangles);
    
    iteration = iteration + 1;
end


%% Waiting for convergence
while( sum(prev_triangle_cells ~= next_triangle_cells) > 10 && iteration < 500  ) % parameter!
    if(showProgress && mod(iteration, 11) == 0)
        M.PlotRemeshing(next_triangle_cells, false);
        drawnow;
    end
    prev_triangle_cells = next_triangle_cells;
    prev_boundaries = next_boundaries;
    [next_boundaries, next_triangle_cells] = ...
            UpdateBoundaries( prev_triangle_cells, prev_boundaries, ...
                triangle_centroids, prev_metrics, prev_centers, F_adj);
    prev_prev_centers = prev_centers;
    prev_centers = next_centers;
    prev_metrics = next_metrics;
    prev_center_triangles = next_center_triangles;
    [next_metrics, next_centers, next_center_triangles] = UpdateCentersAndMetrics( ...
            STLFormat, triangle_centroids, next_triangle_cells, N, F_Area, prev_center_triangles);
    
    if(all(prev_prev_centers(:) == next_centers(:))) % 2-state oscillations
        break;
    end
    iteration = iteration + 1;
end

if(showProgress)
    M.PlotRemeshing(next_triangle_cells, false); % Update to a better view
    drawnow;
end

%% Removing holes and repeating
while( (sum(prev_triangle_cells ~= next_triangle_cells) > 0 && ...
        iteration < 1000) || ...
        sum(next_triangle_cells == 0) > 0)
    if(showProgress && mod(iteration, 11) == 0)
        M.PlotRemeshing(next_triangle_cells, false);
        drawnow;
    end
    prev_triangle_cells = next_triangle_cells;
    next_triangle_cells = RemoveVoronoiUnconnectedCompounds(next_center_triangles, next_triangle_cells, F_adj);
    prev_boundaries = next_boundaries;
    [triangles,~,cells]  = find(next_triangle_cells);
    faces_in_cells = sparse(triangles, cells, ones(size(triangles)), M.NF, N);
    next_boundaries = GetBoundaryTrianglesMatrix(faces_in_cells, F_adj);
    [next_boundaries, next_triangle_cells] = ...
            UpdateBoundaries( next_triangle_cells, next_boundaries, ...
                triangle_centroids, prev_metrics, prev_centers, F_adj);
    prev_prev_centers = prev_centers;
    prev_centers = next_centers;
    prev_metrics = next_metrics;
    prev_center_triangles = next_center_triangles;
    [next_metrics, next_centers, next_center_triangles] = UpdateCentersAndMetrics( ...
            STLFormat, triangle_centroids, next_triangle_cells, N, F_Area, prev_center_triangles);
    
    if(all(prev_prev_centers(:) == next_centers(:))) % 2-state oscillations
        break;
    end
    iteration = iteration + 1;
end

%% Creating the coarse mesh
vertices_in_cells = double((M.VF_adj * faces_in_cells) > 0);
mesh_vertices = find(sum(vertices_in_cells, 2) >= 3); % Hopefully exactly 3
mesh_vertices_in_cells = vertices_in_cells(mesh_vertices, :);   % coarse mesh VF_adj matrix

% Setting the faces matrix, including the correct directions
[row, col] = find(mesh_vertices_in_cells);
max_polygon_size = max(accumarray(col,1));
new_F = NaN(N, max_polygon_size);
num_edges = size(mesh_vertices, 1) + N - M.Chi;
edges = zeros(num_edges + 1, 2);
edge_index = 1;
for i=1:N
    vertices = row(col == i);
    if isempty(vertices)
        continue
    end
    org_vertices = mesh_vertices(vertices);
    total_curr_triangles_indices = find(faces_in_cells(:, i));
    curr_V_adj = sparse(M.F(total_curr_triangles_indices, [1 2 3]), ...
        M.F(total_curr_triangles_indices, [2 3 1]), ...
        ones(size(total_curr_triangles_indices, 1), 3), M.NV, M.NV);
    curr_V_adj = double((curr_V_adj - curr_V_adj') > 0);
    start_v = sparse(org_vertices(1), 1, 1, M.NV, 1);
    ind = 1;
    edges(edge_index, 1) = vertices(1);
    flag = false;
    while( ind <= size(vertices, 1) )
        pos = find(org_vertices == find(start_v));
        if(~isempty(pos))
            new_F(i, ind) = vertices(pos);
            if(flag)
                edges(edge_index, 2) = vertices(pos);
                edge_index = edge_index + 1;
                edges(edge_index, 1) = vertices(pos);
            else
                flag = true;
            end
            ind = ind + 1;
        end
        start_v = curr_V_adj * start_v;
    end
    edges(edge_index, 2) = vertices(1);
    edge_index = edge_index + 1;
end
new_V = M.V(mesh_vertices, :);
edges = unique(sort(edges,2), 'rows'); % Removing edge repetitions

% Drawing the result
if(showProgress)
    hold off;
    patch('Faces',new_F,'Vertices',new_V,'FaceColor','red', 'FaceAlpha', 0.3);
    axis('equal');
    cameratoolbar;
end


%edge_length = vecnorm(new_V(edges(:,1),:) - new_V(edges(:,2),:), 2, 2);

end







