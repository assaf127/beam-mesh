function [triangle_cells] = RemoveVoronoiUnconnectedCompounds( ...
    cell_centers, triangle_cells, F_adj)
%RemoveVoronoiUnconnectedCompounds Removes unconnected components of
% Voronoi cells
%   Removes any unconnected components of the Voronoi cells while keeping
%   the main connected component (which contains the centroid)
% Input:
%   cell_centers    - N x 1, the traingle index of the centroid of each
%                            Voronoi cell
%   triangle_cells  - F x 1, contains the cell number of every triangle
%   F_adj           - F x F, the faces adjacency matrix
% Output:
%   triangle_cells  - F x 1, contains the *updated* cell number of every
%                            triangle. Can be zero if it needs to be
%                            recalculated due to its being in an
%                            unconnected component of a Voronoi cell

for i=1:size(cell_centers, 1)
    curr_triangles = find(triangle_cells == i);
    center = find(curr_triangles == cell_centers(i));
    curr_F_graph = graph(F_adj(curr_triangles, curr_triangles));
    connected_compounds = conncomp(curr_F_graph);
    triangle_cells( ...
        curr_triangles(connected_compounds ~= connected_compounds(center)) ...
        ) = 0;
end

end

