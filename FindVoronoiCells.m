function [cells] = FindVoronoiCells(M, centers, metrics, previous_cells)
%FINDVORONOICELLS Summary of this function goes here
%   Detailed explanation goes here
% Input:
%   M - Mesh3D object
%   centers - The centers of the voronoi cells
%   Metric - 3x3 distance metric
%   previous_cells - A column vector with size M.NF x 1 where each element
%                    is the index of the cell for the previous traingle
% Output:
%   cells - A column vector with size M.NF x 1 where each element is the
%           index of the cell for the traingle

cells = zeros(M.NF,1);
triangle_centroids = (M.V(M.F(:,1),:) + M.V(M.F(:,2),:) + M.V(M.F(:,3),:)) / 3;


parfor i=1:M.NF
    [~, center] = min(FindIsotropicDistance(centers, triangle_centroids(i,:), metrics));
    cells(i,1) = center;
end

end

function [dist] = FindIsotropicDistance(p1, p2, metrics)
    dist = zeros(size(p1,1), 1);
    for i=1:size(p1,1)
        vec = p1(i,:) - p2;
        dist(i,1) = vec * metrics{i} * vec';
    end
end

