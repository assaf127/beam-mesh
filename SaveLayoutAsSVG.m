function SaveLayoutAsSVG(filename, height, width, beam_quads_in1, beam_quads_in2, ...
    beam_quads_out2, beam_quads_out1, hole1, hole2, text1, text2, hole_radius, edges, reference_location)
%SAVELAYOUTASSVG Summary of this function goes here
%   Detailed explanation goes here

%% Save main file
fid = fopen([filename '.svg'],'w', 'n', 'UTF-8');
if fid < 0
    error('Can''t open the file for writing');
end

% Header
fprintf(fid, ['<?xml version="1.0" encoding="utf-8"?>\n' ...
'<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">\n' ...
'<svg version="1.0" id="Layer_1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" x="0px" y="0px" ' ...
'width="%fpx" height="%fpx" viewBox="0 0 %f %f" enable-background="new 0 0 %f %f" ' ...
'xml:space="preserve">\n'], width, height, width, height, width, height);

% Beams
fprintf(fid, '<line x1="%f" y1="%f" x2="%f" y2="%f" stroke="black" stroke-width="1" fill="none" />\n', ...
    [beam_quads_in2(1,:); beam_quads_in2(2,:); beam_quads_in1(1,:); beam_quads_in1(2,:)]);
fprintf(fid, '<line x1="%f" y1="%f" x2="%f" y2="%f" stroke="black" stroke-width="1" fill="none" />\n', ...
    [beam_quads_out2(1,:); beam_quads_out2(2,:); beam_quads_out1(1,:); beam_quads_out1(2,:)]);
fprintf(fid, '<line x1="%f" y1="%f" x2="%f" y2="%f" stroke="black" stroke-width="1" fill="none" />\n', ...
    [beam_quads_out2(1,:); beam_quads_out2(2,:); beam_quads_in2(1,:); beam_quads_in2(2,:)]);
fprintf(fid, '<line x1="%f" y1="%f" x2="%f" y2="%f" stroke="black" stroke-width="1" fill="none" />\n', ...
    [beam_quads_out1(1,:); beam_quads_out1(2,:); beam_quads_in1(1,:); beam_quads_in1(2,:)]);

% Holes
fprintf(fid, ['<circle cx="%f" cy="%f" r="' num2str(hole_radius) '" stroke="black" stroke-width="1" fill="none" />\n' ...
    '<circle cx="%f" cy="%f" r="' num2str(hole_radius) '" stroke="black" stroke-width="1" fill="none" />\n'], ...
    [hole1(1,:); hole1(2,:); hole2(1,:); hole2(2,:)]);

% Labels
fprintf(fid, ['<text x="%f" y="%f" stroke="red" ' ...
    'lengthAdjust="spacingAndGlyphs" stroke-width="1" font-family="''MyriadPro-Regular''" fill="none" ' ...
    'font-size="' num2str(4*hole_radius) '">%d</text>\n' ...
    '<text x="%f" y="%f" stroke="red" '...
    'lengthAdjust="spacingAndGlyphs" stroke-width="1" font-family="''MyriadPro-Regular''" fill="none" '...
    'font-size="' num2str(4*hole_radius) '">%d</text>\n'], ...
    [text2(1,:); text2(2,:); edges(:,2)'; text1(1,:); text1(2,:); edges(:,1)']);

% Reference
deltaX = 50;
deltaY = 50;
delta_mat = ...
    [deltaX,  -deltaX, -deltaX, deltaX; ...
     deltaY,  deltaY,  -deltaY, -deltaY; ...
     -deltaX, -deltaX, deltaX,  deltaX; ...
     deltaY,  -deltaY, -deltaY, deltaY];
delta_mat = delta_mat + repmat(reference_location,2,4);
fprintf(fid, '<line x1="%f" y1="%f" x2="%f" y2="%f" stroke="blue" stroke-width="1" fill="none" />\n', ...
    delta_mat);

% Finish
fprintf(fid, '</svg>');
fclose(fid);

%% Save cutting only file
fid = fopen([filename '_cuts.svg'],'w', 'n', 'UTF-8');
if fid < 0
    error('Can''t open the file for writing');
end

% Header
fprintf(fid, ['<?xml version="1.0" encoding="utf-8"?>\n' ...
'<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">\n' ...
'<svg version="1.0" id="Layer_1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" x="0px" y="0px" ' ...
'width="%fpx" height="%fpx" viewBox="0 0 %f %f" enable-background="new 0 0 %f %f" ' ...
'xml:space="preserve">\n'], width, height, width, height, width, height);

% Beams
fprintf(fid, '<line x1="%f" y1="%f" x2="%f" y2="%f" stroke="black" stroke-width="1" fill="none" />\n', ...
    [beam_quads_in2(1,:); beam_quads_in2(2,:); beam_quads_in1(1,:); beam_quads_in1(2,:)]);
fprintf(fid, '<line x1="%f" y1="%f" x2="%f" y2="%f" stroke="black" stroke-width="1" fill="none" />\n', ...
    [beam_quads_out2(1,:); beam_quads_out2(2,:); beam_quads_out1(1,:); beam_quads_out1(2,:)]);
fprintf(fid, '<line x1="%f" y1="%f" x2="%f" y2="%f" stroke="black" stroke-width="1" fill="none" />\n', ...
    [beam_quads_out2(1,:); beam_quads_out2(2,:); beam_quads_in2(1,:); beam_quads_in2(2,:)]);
fprintf(fid, '<line x1="%f" y1="%f" x2="%f" y2="%f" stroke="black" stroke-width="1" fill="none" />\n', ...
    [beam_quads_out1(1,:); beam_quads_out1(2,:); beam_quads_in1(1,:); beam_quads_in1(2,:)]);

% Holes
fprintf(fid, ['<circle cx="%f" cy="%f" r="' num2str(hole_radius) '" stroke="black" stroke-width="1" fill="none" />\n' ...
    '<circle cx="%f" cy="%f" r="' num2str(hole_radius) '" stroke="black" stroke-width="1" fill="none" />\n'], ...
    [hole1(1,:); hole1(2,:); hole2(1,:); hole2(2,:)]);

% Reference
deltaX = 50;
deltaY = 50;
delta_mat = ...
    [deltaX,  -deltaX, -deltaX, deltaX; ...
     deltaY,  deltaY,  -deltaY, -deltaY; ...
     -deltaX, -deltaX, deltaX,  deltaX; ...
     deltaY,  -deltaY, -deltaY, deltaY];
delta_mat = delta_mat + repmat(reference_location,2,4);
fprintf(fid, '<line x1="%f" y1="%f" x2="%f" y2="%f" stroke="blue" stroke-width="1" fill="none" />\n', ...
    delta_mat);

% Finish
fprintf(fid, '</svg>');
fclose(fid);


%% Save labels only file
fid = fopen([filename '_labels.svg'],'w', 'n', 'UTF-8');
if fid < 0
    error('Can''t open the file for writing');
end

% Header
fprintf(fid, ['<?xml version="1.0" encoding="utf-8"?>\n' ...
'<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">\n' ...
'<svg version="1.0" id="Layer_1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" x="0px" y="0px" ' ...
'width="%fpx" height="%fpx" viewBox="0 0 %f %f" enable-background="new 0 0 %f %f" ' ...
'xml:space="preserve">\n'], width, height, width, height, width, height);

% Labels
fprintf(fid, ['<text x="%f" y="%f" stroke="red" ' ...
    'lengthAdjust="spacingAndGlyphs" stroke-width="1" font-family="''MyriadPro-Regular''" fill="none" ' ...
    'font-size="' num2str(4*hole_radius) '">%d</text>\n' ...
    '<text x="%f" y="%f" stroke="red" '...
    'lengthAdjust="spacingAndGlyphs" stroke-width="1" font-family="''MyriadPro-Regular''" fill="none" '...
    'font-size="' num2str(4*hole_radius) '">%d</text>\n'], ...
    [text2(1,:); text2(2,:); edges(:,2)'; text1(1,:); text1(2,:); edges(:,1)']);

% Reference
deltaX = 50;
deltaY = 50;
delta_mat = ...
    [deltaX,  -deltaX, -deltaX, deltaX; ...
     deltaY,  deltaY,  -deltaY, -deltaY; ...
     -deltaX, -deltaX, deltaX,  deltaX; ...
     deltaY,  -deltaY, -deltaY, deltaY];
delta_mat = delta_mat + repmat(reference_location,2,4);
fprintf(fid, '<line x1="%f" y1="%f" x2="%f" y2="%f" stroke="blue" stroke-width="1" fill="none" />\n', ...
    delta_mat);

% Finish
fprintf(fid, '</svg>');
fclose(fid);


end

