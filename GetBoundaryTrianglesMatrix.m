function [faces_next_to_cells] = GetBoundaryTrianglesMatrix(faces_in_cells, F_adj)
%GetBoundaryTrianglesMatrix Returns a matrix representing the boundary
%triangles for each cell
%   Returns a M.NF x (NumOfCells) matrix where each column represents the
%   triangles that are in the boundary of a specific cell (inside and
%   outside triangles). For each triangle, the distance calculation will
%   only have to take into account the cells that are in his row.
% Input:
%   faces_in_cells  - M.NF x (NumOfCells) matrix representing the triangles
%                     that are in each cell
%   F_adj           - The faces adjacency matrix


%faces_next_to_cells = ismember(F_adj * faces_in_cells, 1:3);
faces_next_to_cells = F_adj * faces_in_cells;
faces_next_to_cells = (faces_next_to_cells == 1) | (faces_next_to_cells == 2) ...
    | (faces_next_to_cells == 3);
end

