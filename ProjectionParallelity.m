function [v_in_res, v_out_res, prev_error] = ...
    ProjectionParallelity(v_in, v_out, edges)
%PROJECTIONPARALLELITY Summary of this function goes here
%   Detailed explanation goes here
V = size(v_in, 1);

v1_centers = (v_in(edges(:,1),:) + v_out(edges(:,1),:)) ./ 2;
v2_centers = (v_in(edges(:,2),:) + v_out(edges(:,2),:)) ./ 2;
edge_norm = v1_centers - v2_centers;
edge_norm = edge_norm ./ vecnorm(edge_norm,2,2);
edge_in_centers = (v_in(edges(:,1),:) + v_in(edges(:,2),:)) ./ 2;
edge_out_centers = (v_out(edges(:,1),:) + v_out(edges(:,2),:)) ./ 2;

half_new_in_edge = dot( (v_in(edges(:,1),:) - v_in(edges(:,2),:)) ./ 2, edge_norm, 2) ...
    .* edge_norm;
half_new_out_edge = dot( (v_out(edges(:,1),:) - v_out(edges(:,2),:)) ./ 2, edge_norm, 2) ...
    .* edge_norm;

projected_v_in = [edge_in_centers + half_new_in_edge; edge_in_centers - half_new_in_edge];
projected_v_out = [edge_out_centers + half_new_out_edge; edge_out_centers - half_new_out_edge];

count_v = accumarray(edges(:), 1, [V 1]);
sum_v_in = full(sparse([edges(:), edges(:), edges(:)], repmat(1:3,size(edges(:),1),1), ...
    projected_v_in, V, 3));
v_in_res = sum_v_in ./ count_v;

sum_v_out = full(sparse([edges(:), edges(:), edges(:)], repmat(1:3,size(edges(:),1),1), ...
    projected_v_out, V, 3));
v_out_res = sum_v_out ./ count_v;

prev_error = norm([v_out;v_in] - [v_out_res;v_in_res], 'fro') .^ 2;

end

