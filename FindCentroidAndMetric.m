function [centers, metrics] = FindCentroidAndMetric(M, cells, cell_num)
%FINDCENTROIDANDMETRIC Returns the centroid and Mahalanobis metric for the
%specified voronoi cells
%   Detailed explanation goes here
% Input:
%   M - Mesh3D object
%   cells - The voronoi cell indices of each triangle (M.NF x 3)
%   cell_num - The number of voronoi cells
% Output:
%   centers - A column vector with size M.NF x 3 where each element is the
%             center of the voronoi cell
%   metrics - A cell array with the Mahalanobis distance metric for each
%             voronoi cell

    metrics = cell(cell_num, 1);
    centers = zeros(cell_num, 3);
    F_Area = M.F_Area;
    covariance_template = [2 1 1;
                           1 2 1;
                           1 1 2] / 12;
    triangle_centroids = (M.V(M.F(:,1),:) + M.V(M.F(:,2),:) + M.V(M.F(:,3),:)) / 3;
    
    for i=1:cell_num
        metrics{i} = zeros(3);
        centers(i,:) = mean(M.V(M.F((cells == i),:),:)); % Maybe should be updated
    end

    for i=1:M.NF
        cell_index = cells(i);
        triangle_vertices = M.V(M.F(i,:),:) - centers(cell_index,:);
        metrics{cell_index} = metrics{cell_index} + F_Area(i) * ...
                        triangle_vertices' * covariance_template * triangle_vertices;
    end

    for i=1:cell_num
        [U, S, V] = svd(metrics{i});
        S(1,1) = S(3,3) / S(1,1);
        S(2,2) = S(3,3) / S(2,2);
        S(3,3) = 1;
        metrics{i} = U * S * V';
        try
        indices = find(cells == i);
        cell_centers = triangle_centroids(indices, :);
        [~, centerIndex] = min(FindIsotropicDistance(cell_centers, centers(i,:)));
        centers(i,:) = triangle_centroids(indices(centerIndex), :);
        catch ME
            5;
        end
    end
end


function [dist] = FindIsotropicDistance(p1, p2)
    dist = dot(p1 - p2, p1 - p2,2);
end
