function [v_in_res, v_out_res, prev_error] = ...
    ProjectionPlanarity(v_in, v_out, edges, edge_length)
%PROJECTIONPLANARITY Summary of this function goes here
%   Detailed explanation goes here
V = size(v_in, 1);
edge_norm = v_in(edges(:,1),:) + v_out(edges(:,1),:) - v_in(edges(:,2),:) - v_out(edges(:,2),:);
edge_norm = edge_norm ./ vecnorm(edge_norm,2,2);
quad_centers = (v_in(edges(:,1),:) + v_out(edges(:,1),:) + ...
    v_in(edges(:,2),:) + v_out(edges(:,2),:)) ./ 4;

quad_norms = cross(v_out(edges(:,1),:) + v_out(edges(:,2),:) ...
    - v_in(edges(:,1),:) - v_in(edges(:,2),:), ...
    edge_norm, 2);
quad_norms = quad_norms ./ vecnorm(quad_norms, 2, 2);

quad_centers = [quad_centers; quad_centers];
quad_norms = [quad_norms; quad_norms];

projected_v_in = v_in(edges(:),:) - ...
    dot(v_in(edges(:),:) - quad_centers, quad_norms, 2) .* quad_norms;
projected_v_out = v_out(edges(:),:) - ...
    dot(v_out(edges(:),:) - quad_centers, quad_norms, 2) .* quad_norms;

% count_v = accumarray(edges(:), 1, [V 1]);
% sum_v_in = full(sparse([edges(:), edges(:), edges(:)], repmat(1:3,size(edges(:),1),1), ...
%     projected_v_in, V, 3));
% v_in_res = sum_v_in ./ count_v;
% 
% sum_v_out = full(sparse([edges(:), edges(:), edges(:)], repmat(1:3,size(edges(:),1),1), ...
%     projected_v_out, V, 3));
% v_out_res = sum_v_out ./ count_v;

count_v = full(sparse(edges, 1, repmat(1 ./ edge_length, 1, 2), V, 1));
sum_v_in = full(sparse([edges(:), edges(:), edges(:)], repmat(1:3,size(edges(:),1),1), ...
    projected_v_in ./ repmat(edge_length, 2, 1), V, 3));
v_in_res = sum_v_in ./ count_v;

sum_v_out = full(sparse([edges(:), edges(:), edges(:)], repmat(1:3,size(edges(:),1),1), ...
    projected_v_out ./ repmat(edge_length, 2, 1), V, 3));
v_out_res = sum_v_out ./ count_v;


prev_error = norm([v_out;v_in] - [v_out_res;v_in_res], 'fro') .^ 2;

end

