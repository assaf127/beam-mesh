function [quad_pos, sheet_height, order] = ...
    LayoutQuads(beam_quads_in1, beam_quads_in2, beam_quads_out2, beam_quads_out1, ...
    space_size, sheet_width)
%LAYOUTQUADS Summary of this function goes here
%   Detailed explanation goes here
[quad_x_sizes, order] = ...
    sort(max([beam_quads_in1(1,:);beam_quads_out1(1,:)]) - min([beam_quads_in2(1,:);beam_quads_out2(1,:)]), 'descend');
quad_y_sizes = max([beam_quads_in1(2,:);beam_quads_in2(2,:)]) - min([beam_quads_out1(2,:);beam_quads_out2(2,:)]);
quad_y_sizes = quad_y_sizes(order);
quad_bounding_center_x = (max([beam_quads_in1(1,:);beam_quads_out1(1,:)]) + ...
    min([beam_quads_in2(1,:);beam_quads_out2(1,:)])) ./ 2;
quad_bounding_center_y = (max([beam_quads_in1(2,:);beam_quads_in2(2,:)]) + ...
    min([beam_quads_out1(2,:);beam_quads_out2(2,:)])) ./ 2;
quad_bounding_center = [quad_bounding_center_x;quad_bounding_center_y];
quad_bounding_center = quad_bounding_center(:,order);

quad_pos = zeros(size(beam_quads_in1));
curr_x = 0;
curr_y = 0;
max_row_height = 0;
for i=1:size(beam_quads_in1,2)
    if curr_x + quad_x_sizes(i) > sheet_width
        curr_x = 0;
        curr_y = curr_y + max_row_height + space_size;
        max_row_height = 0;
    end
    curr_height = quad_y_sizes(i);
    if curr_height > max_row_height
        max_row_height = curr_height;
    end
    quad_pos(:, order(i)) = [curr_x;curr_y] + ([quad_x_sizes(i);quad_y_sizes(i)] ./ 2) ...
        - quad_bounding_center(:,i);
    curr_x = curr_x + quad_x_sizes(i) + space_size;
end

sheet_height = curr_y + max_row_height;

end

