function [h] = PlotBeamMesh(v_out, v_in, edges)
%PlotBeamMesh Summary of this function goes here
%   Detailed explanation goes here
V = size(v_in, 1);
quads = [edges, edges(:,[2 1])+V];
vertices = [v_in; v_out];

h = patch('Faces',quads,'Vertices',vertices,'FaceColor',[0.3 0.3 0.3]);
axis('equal');
cameratoolbar;

end

