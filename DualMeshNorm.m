function [norm_V] = DualMeshNorm(V,F)
%DUALMESHNORM Summary of this function goes here
%   Detailed explanation goes here
NV = size(V, 1);
NF = size(F, 1);
norm_V = zeros(NV,3);
for i=1:NF
    current_edges_indices = find(~isnan(F(i,:)));
    if isempty(current_edges_indices)
        end_index = size(F, 2);
    else
        end_index = current_edges_indices(end);
    end
    current_edges_indices = F(i,1:end_index);
    norm_V(current_edges_indices(1), :) = ...
        norm_V(current_edges_indices(1), :) + ...
        cross(V(current_edges_indices(end),:) - V(current_edges_indices(1),:), ...
        V(current_edges_indices(2),:) - V(current_edges_indices(1),:), 2);
    for j=2:end_index-1
        norm_V(current_edges_indices(j), :) = ...
        norm_V(current_edges_indices(j), :) + ...
        cross(V(current_edges_indices(j-1),:) - V(current_edges_indices(j),:), ...
        V(current_edges_indices(j+1),:) - V(current_edges_indices(j),:), 2);
    end
    norm_V(current_edges_indices(end), :) = ...
    norm_V(current_edges_indices(end), :) + ...
    cross(V(current_edges_indices(end-1),:) - V(current_edges_indices(end),:), ...
    V(current_edges_indices(1),:) - V(current_edges_indices(end),:), 2);
end
norm_V = norm_V ./ vecnorm(norm_V, 2, 2);
end

