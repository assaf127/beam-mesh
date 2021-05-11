function [beam_quads_in1, beam_quads_in2, beam_quads_out2, beam_quads_out1, ...
    hole1, hole2, text1, text2] = ...
    CreateBeams(v_in, v_out, edges, hole_radius)
%CREATEBEAMS Summary of this function goes here
%   Detailed explanation goes here
v = (v_out + v_in) ./ 2;
E = size(edges, 1);
quad_centers = (v_in(edges(:,1),:) + v_out(edges(:,1),:) + ...
    v_in(edges(:,2),:) + v_out(edges(:,2),:)) ./ 4;
edges_vec = v(edges(:,1),:) -  v(edges(:,2),:);
edges_length = vecnorm(edges_vec,2,2);
edge_norm = edges_vec ./ edges_length;

quad_norms = cross(v_out(edges(:,1),:) + v_out(edges(:,2),:) ...
    - v_in(edges(:,1),:) - v_in(edges(:,2),:), ...
    edge_norm, 2);
quad_norms = quad_norms ./ vecnorm(quad_norms, 2, 2);

rotation_vec = cross(quad_norms, repmat([0 0 1],E,1), 2);
rotation_vec = rotation_vec ./ vecnorm(rotation_vec, 2, 2);
rotation_angle = GetVectorsAngle(quad_norms, repmat([0 0 1],E,1));

% pre-allocation
beam_quads_in1 = zeros(2,E);
beam_quads_in2 = zeros(2,E);
beam_quads_out2 = zeros(2,E);
beam_quads_out1 = zeros(2,E);
hole1 = zeros(2,E);
hole2 = zeros(2,E);
text1 = zeros(2,E);
text2 = zeros(2,E);

for i=1:E
    curr_rotation = vrrotvec2mat([rotation_vec(i,:), rotation_angle(i)]);
    curr_rotation = curr_rotation([1 2],:);
    beam_quads_in1(:,i) = curr_rotation * (v_in(edges(i,1),:) - quad_centers(i,:))';
    beam_quads_in2(:,i) = curr_rotation * (v_in(edges(i,2),:) - quad_centers(i,:))';
    beam_quads_out2(:,i) = curr_rotation * (v_out(edges(i,2),:) - quad_centers(i,:))';
    beam_quads_out1(:,i) = curr_rotation * (v_out(edges(i,1),:) - quad_centers(i,:))';
    
    edge_vec = (beam_quads_in1(:,i) + beam_quads_out1(:,i) - ...
        beam_quads_in2(:,i) - beam_quads_out2(:,i)) ./ 2;
    edge_norm = edge_vec ./ vecnorm(edge_vec);
    
    curr_rotation = vrrotvec2mat([0,0,1,-atan2(edge_norm(2), edge_norm(1))]);
    curr_rotation = curr_rotation(1:2, 1:2);
    beam_quads_in1(:,i) = curr_rotation * beam_quads_in1(:,i);
    beam_quads_in2(:,i) = curr_rotation * beam_quads_in2(:,i);
    beam_quads_out2(:,i) = curr_rotation * beam_quads_out2(:,i);
    beam_quads_out1(:,i) = curr_rotation * beam_quads_out1(:,i);
    hole1(1,i) = min([beam_quads_in1(1,i) beam_quads_out1(1,i)]) - 2 * hole_radius;
    hole2(1,i) = max([beam_quads_in2(1,i) beam_quads_out2(1,i)]) + 2 * hole_radius;
    text1(1,i) = hole2(1,i) - hole_radius;
    text2(1,i) = text1(1,i);
    text1(2,i) = -1.2 * hole_radius;
    text2(2,i) = min([beam_quads_in1(2,i) beam_quads_in2(2,i)]) - hole_radius;
end


end

