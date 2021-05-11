function [v_in_res, v_out_res, affected_v, prev_error] = ...
    ProjectionHeight(v_in,v_out, edges, min_height, max_height)
%PROJECTIONHEIGHT Summary of this function goes here
%   Detailed explanation goes here
V = size(v_in, 1);
v_in_res = v_in;
v_out_res = v_out;
v_vectors = v_out - v_in;
edge_norm = v_in(edges(:,1),:) + v_out(edges(:,1),:) - v_in(edges(:,2),:) - v_out(edges(:,2),:);
edge_norm = edge_norm ./ vecnorm(edge_norm,2,2);
edges = edges(:);
height_vec = cross(cross([-edge_norm; edge_norm], v_vectors(edges, :), 2), [-edge_norm; edge_norm], 2);
current_height = vecnorm(height_vec, 2, 2);

short_edges = find(current_height < min_height);
tall_edges = find(current_height > max_height);

short_v = edges(short_edges);
tall_v = edges(tall_edges);
affected_v = zeros(V,1);
affected_v([short_v;tall_v]) = 1;


% 
% vertices_vectors = 0.5 .* (v_out - v_in);
% vertices_centers = (v_out + v_in) ./ 2;
% 
% delta_height_vector = (min_height ./ current_height(short_edges)) .* vertices_vectors(short_v, :);
% new_short_edges_out = vertices_centers(short_v, :) + delta_height_vector;
% new_short_edges_in = vertices_centers(short_v, :) - delta_height_vector;
% delta_height_vector = (max_height ./ current_height(tall_edges)) .* vertices_vectors(tall_v, :);
% new_tall_edges_out = vertices_centers(tall_v, :) + delta_height_vector;
% new_tall_edges_in = vertices_centers(tall_v, :) - delta_height_vector;
% 
% count_v = accumarray([edges(short_edges);edges(tall_edges)], 1, [V 1]);
% sum_v_in = full(sparse([repmat(short_v,1,3);repmat(tall_v,1,3)], ...
%     repmat(1:3,size(short_v,1)+size(tall_v,1),1), ...
%     [new_short_edges_in;new_tall_edges_in], V, 3));
% sum_v_out = full(sparse([repmat(short_v,1,3);repmat(tall_v,1,3)], ...
%     repmat(1:3,size(short_v,1)+size(tall_v,1),1), ...
%     [new_short_edges_out;new_tall_edges_out], V, 3));
% update_v = count_v > 0;
% v_in_res(update_v, :) = sum_v_in(update_v, :) ./ count_v(update_v);
% v_out_res(update_v, :) = sum_v_out(update_v, :) ./ count_v(update_v);


height_norm = height_vec ./ vecnorm(height_vec, 2, 2);
delta_height_vector = 0.5 * (min_height - current_height(short_edges)) .* height_norm(short_edges, :);
new_short_edges_out = v_out(edges(short_edges), :) + delta_height_vector;
new_short_edges_in = v_in(edges(short_edges), :) - delta_height_vector;
delta_height_vector = 0.5 * (max_height - current_height(tall_edges)) .* height_norm(tall_edges, :);
new_tall_edges_out = v_out(edges(tall_edges), :) + delta_height_vector;
new_tall_edges_in = v_in(edges(tall_edges), :) - delta_height_vector;

counts = zeros(V, 1);
sums_in = zeros(V, 3);
sums_out = zeros(V, 3);
for i=1:size(short_edges,1)
    index = edges(short_edges(i));
    counts(index) = counts(index) + 1;
    sums_out(index, :) = sums_out(index, :) + new_short_edges_out(i, :);
    sums_in(index, :) = sums_in(index, :) + new_short_edges_in(i, :);
end
for i=1:size(tall_edges,1)
    index = edges(tall_edges(i));
    counts(index) = counts(index) + 1;
    sums_out(index, :) = sums_out(index, :) + new_tall_edges_out(i, :);
    sums_in(index, :) = sums_in(index, :) + new_tall_edges_in(i, :);
end

changed_vs = counts > 0;
v_in_res(changed_vs, :) = sums_in(changed_vs, :) ./ counts(changed_vs);
v_out_res(changed_vs, :) = sums_out(changed_vs, :) ./ counts(changed_vs);

prev_error = norm([v_out;v_in] - [v_out_res;v_in_res], 'fro') .^ 2;

end

