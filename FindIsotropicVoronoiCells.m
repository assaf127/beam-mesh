function [cells] = FindIsotropicVoronoiCells(M, centers)
%FINDISOTROPICVORONOICELLS Finds the isotropic (Metric = Identity) voronoi
%cells
%   Detailed explanation goes here
% Input:
%   M - Mesh3D object
%   centers - The centers of the voronoi cells (M.NF x 3)
% Output:
%   cells - A column vector with size M.NF x 1 where each element is the
%           index of the cell for the traingle

cells = zeros(M.NF,1);
triangle_centroids = (M.V(M.F(:,1),:) + M.V(M.F(:,2),:) + M.V(M.F(:,3),:)) / 3;


for i=1:M.NF
    [~, center] = min(FindIsotropicDistance(centers, triangle_centroids(i,:)));
    cells(i,1) = center;
end

end

function [dist] = FindIsotropicDistance(p1, p2)
    dist = dot(p2 - p1, p2 - p1, 2);
end

