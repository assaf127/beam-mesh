function [v_res, affected_v, prev_error] = ProjectionMinLength(v, edges, min_length)
%PROJECTIONMINLENGTH Summary of this function goes here
%   Detailed explanation goes here
v_res = v;
edge_length = vecnorm(v(edges(:,1), :) - v(edges(:,2), :), 2, 2);
short_edges = find(edge_length < min_length);

edited_edges = edges(short_edges, :);
affected_v = zeros(size(v,1),1);
affected_v(edited_edges(:)) = 1;
new_edge_halves = 0.5 * min_length * ...
                (v(edited_edges(:,1), :) - v(edited_edges(:,2), :)) ./ ...
                edge_length(short_edges);
edge_centers = (v(edited_edges(:,1), :) + v(edited_edges(:,2), :)) ./ 2;

new_vs = [edge_centers + new_edge_halves; edge_centers - new_edge_halves];

v_counts = accumarray(edited_edges(:), 1, [size(v,1) 1]);
v_sums = sparse([edited_edges(:) edited_edges(:) edited_edges(:)], ...
    repmat(1:3,size(edited_edges(:), 1),1), ...
    new_vs, size(v,1), 3);
updated_locations = v_counts > 0;
v_res(updated_locations, :) = v_sums(updated_locations, :) ./ v_counts(updated_locations);

prev_error = norm(v - v_res, 'fro') .^ 2;
end

