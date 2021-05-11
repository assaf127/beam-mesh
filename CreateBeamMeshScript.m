
clear;
M = Mesh3D('bunny2.off');
N = 50;
max_dimension_length_cm = 30;
hole_radius = 0.1;
sheet_width = 50;
output_filename = 'result';
%% Performing Lloyd algorithm for obtaining the CVD
[first_V, first_F, ~, ~] = ...
    CoarseRemeshMahalanobisCVD(M, N, true);

bounding_box_size = max(first_V) - min(first_V);
stretch_factor = max_dimension_length_cm / max(bounding_box_size);
first_V = stretch_factor * first_V;

optimization_parameters.min_length = 2;
optimization_parameters.min_height = 1.0;
optimization_parameters.max_height = 1.2;
optimization_parameters.min_length_weight = 10;
optimization_parameters.normal_direction_weight = 0.2;
optimization_parameters.height_weight = 20;
optimization_parameters.planarity_weight = 200;
optimization_parameters.parallelity_weight = 0.4;


[V, F, v_in, v_out, edges, error] = ...
    BeamMeshOptimization(first_V, first_F, optimization_parameters);

PlotBeamMesh(v_out, v_in, edges);

% %% Beam mesh evaluation
% edges_vec = V(edges(:,1),:) -  V(edges(:,2),:);
% edges_length = vecnorm(edges_vec,2,2);
% edge_norm = edges_vec ./ edges_length;
% v_vectors = v_out - v_in;
% norm_V = v_vectors ./ vecnorm(v_vectors,2,2);
% cross1 = cross(edges_vec, norm_V(edges(:,1),:),2);
% cross2 = cross(edges_vec, norm_V(edges(:,2),:),2);
% torsion_angles = GetVectorsAngleD(cross1, cross2);
% height_vec = cross(cross([-edge_norm; edge_norm], v_vectors(edges, :), 2), [-edge_norm; edge_norm], 2);
% edges_height = vecnorm(height_vec, 2, 2);

[beam_quads_in1, beam_quads_in2, beam_quads_out2, beam_quads_out1, ...
    hole1, hole2, text1, text2] = ...
    CreateBeams(v_in, v_out, edges, hole_radius);

[quad_pos, sheet_height, order] = ...
    LayoutQuads(beam_quads_in1, beam_quads_in2, beam_quads_out2, beam_quads_out1, ...
    hole_radius, sheet_width);

new_beam_quads_in1 = beam_quads_in1 + quad_pos;
new_beam_quads_in2 = beam_quads_in2 + quad_pos;
new_beam_quads_out2 = beam_quads_out2 + quad_pos;
new_beam_quads_out1 = beam_quads_out1 + quad_pos;
new_hole1 = hole1 + quad_pos;
new_hole2 = hole2 + quad_pos;
new_text1 = text1 + quad_pos;
new_text2 = text2 + quad_pos;


theta = linspace(0,2*pi,100)';
circleX = hole_radius * cos(theta);
circleY = hole_radius * sin(theta);
figure;
X = [new_beam_quads_in1(1,:);new_beam_quads_in2(1,:); ...
     new_beam_quads_out2(1,:);new_beam_quads_out1(1,:);new_beam_quads_in1(1,:)];
Y = [new_beam_quads_in1(2,:);new_beam_quads_in2(2,:); ...
     new_beam_quads_out2(2,:);new_beam_quads_out1(2,:);new_beam_quads_in1(2,:)];
plot(X,Y);
hold on;
plot(repmat(circleX,1,size(new_hole1,2)) + repmat(new_hole1(1,:),size(circleX,1),1), ...
    repmat(circleY,1,size(new_hole1,2)) + repmat(new_hole1(2,:),size(circleY,1),1));
plot(repmat(circleX,1,size(new_hole2,2)) + repmat(new_hole2(1,:),size(circleX,1),1), ...
    repmat(circleY,1,size(new_hole2,2)) + repmat(new_hole2(2,:),size(circleY,1),1));

% plot(repmat(circleX,1,size(new_text1,2)) + repmat(new_text1(1,:),size(circleX,1),1), ...
%     repmat(circleY,1,size(new_text1,2)) + repmat(new_text1(2,:),size(circleY,1),1));
% plot(repmat(circleX,1,size(new_text2,2)) + repmat(new_text2(1,:),size(circleX,1),1), ...
%     repmat(circleY,1,size(new_text2,2)) + repmat(new_text2(2,:),size(circleY,1),1));
axis equal;



reference_location = [50*sheet_width; 100*sheet_height + 100];

SaveLayoutAsSVG(output_filename, 100*sheet_height + 300, 100*sheet_width, 100*new_beam_quads_in1(:,order), ...
    100*new_beam_quads_in2(:,order), 100*new_beam_quads_out2(:,order), ...
    100*new_beam_quads_out1(:,order), 100*new_hole1(:,order), 100*new_hole2(:,order), ...
    100*new_text1(:,order), 100*new_text2(:,order), 100*hole_radius, edges(order,:), ...
    reference_location);

%% Creating a construction list for the vertices conections
edges_sheet = zeros(size(F,1) * (size(F,2)),3);
curr_index = 1;
for i=1:size(F,1)
    curr_row = F(i, ~isnan(F(i,:)));
    curr_size = size(curr_row,2);
    for j=1:(curr_size-2)
        edges_sheet(curr_index, :) = F(i, j:(j+2));
        curr_index = curr_index + 1;
    end
    edges_sheet(curr_index, :) = F(i, [curr_size-1,curr_size,1]);
    curr_index = curr_index + 1;
    edges_sheet(curr_index, :) = F(i, [curr_size,1,2]);
    curr_index = curr_index + 1;
end
edges_sheet(curr_index:end, :) = [];

[~,permutation] = sort(edges_sheet(:,2));
edges_sheet = edges_sheet(permutation, :);
vertex_layout = zeros(size(V,1),3);
for i=1:(size(edges_sheet,1)/3)
    curr_mat = sortrows(edges_sheet((1+3*(i-1)):(3*i), [1 3]));
    vertex_layout(i,1:2) = curr_mat(1,:);
    vertex_layout(i,3) = curr_mat(curr_mat(:,1) == curr_mat(1,2),2);
end
vertex_layout = sortrows(vertex_layout);


PlotBeamMesh(v_out, v_in, edges);
hold on;
scatter3(V(1,1),V(1,2),V(1,3),'k');
scatter3(V(8,1),V(8,2),V(8,3),'r');
scatter3(V(36,1),V(36,2),V(36,3),'g');
scatter3(V(52,1),V(52,2),V(52,3),'b');

