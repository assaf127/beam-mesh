function [new_V, new_F, v_in, v_out, edges, error] = ...
    BeamMeshOptimization(V, F, optimization_parameters)
%BeamMeshOptimization Optimizes a dual mesh for beam extraction
%   Optimizes a dual triangle mesh for extracting beams for production.
%   The optimization is done based on ShapeUp with projections for minimal
%   beam length, optimal beam height, vertices normal direction and beam
%   planarity, with specific weights for each projection.
% Input:
%   V           - Vertices (V x 3) matrix with (x,y,z) coordinates
%   F           - Faces (F x ?) of the mesh
%   optimization_parameters: struct containing the following fields:
%       min_length
%       min_height
%       max_height
%       min_length_weight
%       normal_direction_weight
%       height_weight
%       planarity_weight
%       parallelity_weight
% Output:
%   new_V   - Optimized V matrix (size might change from the input)
%   new_F   - Optimized V matrix (size might change from the input)
%             mesh. Unused areas are filled with NaNs
%   v_in    - Locations of the inside vertices of the beam mesh
%   v_out   - Locations of the outside vertices of the beam mesh
%   edges   - E x 2 matrix where each row contains 2 vertex indices
%             representing an edge in the coarse mesh
%   error   - The sum of squared differences from the projections

% clear;
% M = Mesh3D('cat_s3.off');
% N = 120;
% %% Performing Lloyd algorithm for obtaining the CVD
% [new_V, new_F, ~, ~] = ...
%     CoarseRemeshMahalanobisCVD(M, N, true);

[new_V,new_F, edges] = CleanDualMesh(V, F);

%% Setting up constants and parameters
NV = size(new_V,1);
NF = size(new_F,1);
edge_length = vecnorm(new_V(edges(:,1),:) - new_V(edges(:,2),:), 2, 2);
norm_edges = (new_V(edges(:,2),:) - new_V(edges(:,1),:)) ./ edge_length;

norm_V = DualMeshNorm(new_V, new_F);
% Constraint parameters
min_length  = optimization_parameters.min_length;
min_height  = optimization_parameters.min_height;
max_height  = optimization_parameters.max_height;
height      = (min_height + max_height) / 2;
% Constraint weight parameters
min_length_weight       = optimization_parameters.min_length_weight;
normal_direction_weight = optimization_parameters.normal_direction_weight;
height_weight           = optimization_parameters.height_weight;
planarity_weight        = optimization_parameters.planarity_weight;
parallelity_weight      = optimization_parameters.parallelity_weight;

total_weight = min_length_weight + normal_direction_weight + ...
    height_weight + planarity_weight + parallelity_weight;

%% Starting guess
v_out = new_V + (0.5 * height) .* norm_V;
v_in = new_V - (0.5 * height) .* norm_V;

% PlotBeamMesh(v_out, v_in, edges);

weights_mat = zeros(NV, 5);
weights_mat(:,1) = min_length_weight;
weights_mat(:,2) = normal_direction_weight;
weights_mat(:,3) = height_weight;
weights_mat(:,4) = planarity_weight;
weights_mat(:,5) = parallelity_weight;

num_iterations = 2000;
last_errors.min_length = zeros(num_iterations,1);
last_errors.normal_direction = zeros(num_iterations,1);
last_errors.height = zeros(num_iterations,1);
last_errors.planarity = zeros(num_iterations,1);
last_errors.parallelity = zeros(num_iterations,1);

for i=1:num_iterations
    % Min length projection
    [v_in2, affected_v_in_min_length, prev_min_length_in_error] = ...
        ProjectionMinLength(v_in, edges, min_length);
    [v_out2, affected_v_out_min_length, prev_min_length_out_error] = ...
        ProjectionMinLength(v_out, edges, min_length);
    prev_min_length_error = prev_min_length_in_error + prev_min_length_out_error;

    % Normal direction projection
    [v_in3, v_out3, prev_normal_direction_error] = ...
        ProjectionNormalDirection(v_in, v_out, norm_V);

    % Height projection
    [v_in4, v_out4, affected_v_height, prev_height_error] = ...
        ProjectionHeight(v_in,v_out, edges, min_height, max_height);

    % Planarity projection
    [v_in5, v_out5, prev_planarity_error] = ...
        ProjectionPlanarity(v_in, v_out, edges, edge_length);

    % Parallelity projection
    [v_in6, v_out6, prev_parallelity_error] = ...
        ProjectionParallelity(v_in, v_out, edges);

    % Average the projections
    current_weights = weights_mat;
    current_weights(~affected_v_height, 3) = 0;
    current_weights(~affected_v_in_min_length, 1) = 0;
    
    v_in = (current_weights(:, 1) .* v_in2 + current_weights(:, 2) .* v_in3 + ...
        current_weights(:, 3) .* v_in4 + current_weights(:, 4) .* v_in5 + ...
        current_weights(:, 5) .* v_in6) ./ sum(current_weights,2);
    
    current_weights(:,1) = min_length_weight;
    current_weights(~affected_v_out_min_length, 1) = 0;
    v_out = (current_weights(:, 1) .* v_out2 + current_weights(:, 2) .* v_out3 + ...
        current_weights(:, 3) .* v_out4 + current_weights(:, 4) .* v_out5 + ...
        current_weights(:, 5) .* v_out6) ./ sum(current_weights,2);

    last_errors.min_length(i) = prev_min_length_error;
    last_errors.normal_direction(i) = prev_normal_direction_error;
    last_errors.height(i) = prev_height_error;
    last_errors.planarity(i) = prev_planarity_error;
    last_errors.parallelity(i) = prev_parallelity_error;
    
    % PlotBeamMesh(v_out, v_in, edges);
    % drawnow;
end

% PlotBeamMesh(v_out, v_in, edges);
% drawnow;

total_error = (last_errors.min_length * min_length_weight + ...
    last_errors.normal_direction * normal_direction_weight + ...
    last_errors.height * height_weight + ...
    last_errors.planarity * planarity_weight + ...
    last_errors.parallelity * parallelity_weight) ./ total_weight;
% figure;
% plot(1:num_iterations, last_errors.min_length, ...
%     1:num_iterations, last_errors.normal_direction, ...
%     1:num_iterations, last_errors.height, ...
%     1:num_iterations, last_errors.planarity, ...
%     1:num_iterations, last_errors.parallelity, ...
%     1:num_iterations, total_error);
% legend('min length','normal direction', 'height', 'planarity', 'parallelity', 'total error');
% title('Errors vs. iterations');
% xlabel('iterations');
% ylabel('Error (MSE)');

error = total_error(end);

new_V = (v_in + v_out) ./ 2;

end

