function [v_in_res, v_out_res, prev_error] = ProjectionNormalDirection(v_in, v_out, normal_direction)
%PROJECTIONNORMALDIRECTION Summary of this function goes here
%   Detailed explanation goes here
centers = (v_in + v_out) ./ 2;

v_out_res = centers - dot(centers - v_out, normal_direction, 2) .* normal_direction;
v_in_res = centers - dot(centers - v_in, normal_direction, 2) .* normal_direction;

prev_error = norm([v_out;v_in] - [v_out_res;v_in_res], 'fro') .^ 2;
end

