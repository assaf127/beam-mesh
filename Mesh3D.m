classdef Mesh3D
    %3DMESH Represents a 3D triagonal mesh
    %   This class contains the (x,y,z) locations of the mesh vertices, as
    %   well as the mesh's faces.
    %   Provided functions:
    %   Reading and writing of OFF mesh files
    %   Visualization of the 3D mesh via the function ShowMesh
    %   Computation of Vertex-Vertex and Vertex-Face adjancy matrices (sparse)
    %   Computation of vertex/face areas
    %   Visualization of a function over the vertices/faces
    %   Interpolation of functions from faces to vertices and vice-versa
    %   Findingthe number of boundary edges and boundary components
    %   Calculation of topological measures: Euler characteristic (Chi) and
    %   the genus
    
    properties
        V       % |V|x3 matrix containing the (x,y,z) coordinates of the vertices
        F       % |F|x3 matrix containing the faces' triangles in a CCW manner
    end
    
    properties (Dependent)
        NV          % Number of vertices
        NF          % Number of faces
        NE          % Number of edges
        V_adj       % Vertex adjacency matrix (sparse)
        VF_adj      % Face-Vertex connectivity matrix (sparse)
        V_Area      % Returns a NV x 1 vector with the vertices area
        F_Area      % Returns a NF x 1 vector with the faces area
        IF2V        % Returns a matrix that interpolates from faces to
                    % vertices, i.e, multiplying a vector containing values
                    % for the faces A like IF2V*A will result in a
                    % vector containing values for the vertices
        IV2F        % Performs the reverse operation of IF2V
        Genus       % Returns the genus of the mesh
        Chi         % Returns the Euler characteristic of the mesh
        NB          % Returns the number of boundary edges in the mesh
        NBComps     % Returns the number of boundary components in the mesh
        Boundaries  % Returns the boundary edges in the mesh as a NB x 2
                    % matrix where each row contains the indices of the
                    % vertices of the boundary edge
                    
        F_norm      % Returns a NF x 3 vector containing the normal to the 
                    % faces (normalized)
        V_norm      % Returns a NV x 3 vector containing the normal to the 
                    % vertices (normalized)
        
                    
        VerticesValence         % Returns a column vector containing the 
                                % degree of each vertex
        GaussianCurvature       % Returns the gaussian curvature on the
                                % vertices as a column vector
                                
        PrependicularEdges      % Returns 3NF x 3 vector with rotated
                                % triangle endges. First NF rows are for
                                % the rotated edge opposite of vertex 1,
                                % and so on. The vectors are directed into
                                % the triangle
        Gradient    % Returns a matrix representing the gradient operator
                    % on a function defined on the vertices. The matrix
                    % size is 3NF x NV, and the result should be reshaped
                    % ( reshape(result,size(A.F)) ) in order to obtain the
                    % gradient vector field
        Divergence  % 
        Laplacian   % 
        
        MeanCurvature           % 
        CotangentWeights        %
        
        F_adj       % Faces adjacency matrix
        
    end
    
    methods
        %%%%% Constructors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = Mesh3D(arg1, arg2)
            if (nargin == 2 && ismatrix(arg1) && ismatrix(arg2)) % Mesh3D(V, F)
                obj.V = arg1;
                if (size(obj.V,2) ~= 3)
                    error('Bad vertices format');
                end
                obj.F = arg2;
                if (size(obj.F,2) ~= 3)
                    error('Bad faces format');
                end
                %obj.NV = size(obj.V,1);
                %obj.NF = size(obj.F,1);
                if (max(obj.F(:)) > obj.NV)
                    error('Faces contain vertices that don''t exist');
                end
                return;
            elseif (nargin == 1)  % Mesh3D(Filename)
                if(isfile(arg1))
                    filetype = lower(arg1((find(arg1 == '.', 1, 'last')+1):end));
                    switch (filetype)
                        case 'off'
                            obj = Mesh3D.ReadOFFFile(arg1);
                        % TODO: Add here more filetypes loading support
                    end
                    return;
                end
                error(['Unable to open file ' arg1]);
            end
        end
        
        %%%%% Getters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function V_adj = get.V_adj(obj)
            if (isnan(obj.V))
                error('Mesh3D object is not initialized with a mesh');
            end
            V_adj = sparse(obj.F,obj.F(:,[2,3,1]),ones(obj.NF,3),obj.NV,obj.NV);
            V_adj = V_adj + speye(size(V_adj));
            V_adj = double((V_adj + V_adj') > 0);
        end
        function VF_adj = get.VF_adj(obj)
            if (isnan(obj.V))
                error('Mesh3D object is not initialized with a mesh');
            end
            VF_adj = sparse(obj.F,[(1:obj.NF)' (1:obj.NF)' (1:obj.NF)'],ones(obj.NF,3),obj.NV,obj.NF);
        end
        function NV = get.NV(obj)
            if (isnan(obj.V))
                error('Mesh3D object is not initialized with a mesh');
            end
            NV = size(obj.V,1);
        end
        function NF = get.NF(obj)
            if (isnan(obj.F))
                error('Mesh3D object is not initialized with a mesh');
            end
            NF = size(obj.F,1);
        end
        function F_Area = get.F_Area(obj)
            if (isnan(obj.F))
                error('Mesh3D object is not initialized with a mesh');
            end
            edge_vec12 = obj.V(obj.F(:,2),:) - obj.V(obj.F(:,1),:);
            edge_vec13 = obj.V(obj.F(:,3),:) - obj.V(obj.F(:,1),:);
            F_Area = 0.5*vecnorm(cross(edge_vec12,edge_vec13)')';   % Half parallelogram area
        end
        function V_Area = get.V_Area(obj)
            if (isnan(obj.V))
                error('Mesh3D object is not initialized with a mesh');
            end
            V_Area = obj.VF_adj*obj.F_Area/3;
        end
        function IF2V = get.IF2V(obj)
            if (isnan(obj.V))
                error('Mesh3D object is not initialized with a mesh');
            end
            IF2V = sparse(obj.F,[(1:obj.NF)' (1:obj.NF)' (1:obj.NF)'], ...
                repmat(obj.F_Area,1,3)./(3*obj.V_Area(obj.F)),obj.NV,obj.NF);
        end
        function IV2F = get.IV2F(obj)
            if (isnan(obj.V))
                error('Mesh3D object is not initialized with a mesh');
            end
            IV2F = sparse([(1:obj.NF)' (1:obj.NF)' (1:obj.NF)'],obj.F, ...
                ones(size(obj.F))/3,obj.NF,obj.NV);
        end
        function NB = get.NB(obj)
            if (isnan(obj.V))
                error('Mesh3D object is not initialized with a mesh');
            end
            V_adj_ex = sparse(obj.F,obj.F(:,[2,3,1]),ones(obj.NF,3),obj.NV,obj.NV);
            V_adj_ex = V_adj_ex + V_adj_ex';
            NB = full(sum(V_adj_ex==1,'all')/2);
        end
        function Boundaries = get.Boundaries(obj)
            if (isnan(obj.V))
                error('Mesh3D object is not initialized with a mesh');
            end
            V_adj_ex = sparse(obj.F,obj.F(:,[2,3,1]),ones(obj.NF,3),obj.NV,obj.NV);
            V_adj_ex_sym = V_adj_ex + V_adj_ex';
            
            B_adj = (V_adj_ex_sym == 1) & (V_adj_ex == 1);
            [row, ~] = find(B_adj);
            if(isempty(row))    % No boundaries
                Boundaries = {};
                return;
            end
            i = 1;
            current_v = row(1);
            e = sparse(obj.NV,1);
            e(current_v) = 1;
            B_adj(:,current_v) = 0;
            Boundaries{i,1} = current_v;
            while (1)
                e = B_adj' * e;
                current_v = find(e);
                if(isempty(current_v)) % finished the boundary
                    [remaining_rows,~] = find(B_adj);
                    if(isempty(remaining_rows))
                        break;
                    else
                        current_v = remaining_rows(1);
                        e(current_v) = 1;
                        i = i + 1;
                        Boundaries{i,1} = [];
                    end
                end
                Boundaries{i,1} = [Boundaries{i,1}, current_v];
                B_adj(:,current_v) = 0;
            end
        end
        function NBComps = get.NBComps(obj)
            if (isnan(obj.V))
                error('Mesh3D object is not initialized with a mesh');
            end
            NBComps = size(obj.Boundaries,1);
        end
        function NE = get.NE(obj)
            if (isnan(obj.V))
                error('Mesh3D object is not initialized with a mesh');
            end
            NE = (obj.NF*3 + obj.NB)/2;
        end
        function Chi = get.Chi(obj)
            if (isnan(obj.V))
                error('Mesh3D object is not initialized with a mesh');
            end
            Chi = obj.NV - obj.NE + obj.NF;
        end
        function Genus = get.Genus(obj)
            if (isnan(obj.V))
                error('Mesh3D object is not initialized with a mesh');
            end
            Genus = 1 - (obj.Chi + obj.NBComps)/2;
        end
        function VerticesValence = get.VerticesValence(obj)
            if (isnan(obj.V))
                error('Mesh3D object is not initialized with a mesh');
            end
            VerticesValence = sum(obj.V_adj,2)-1;
        end
        function GaussianCurvature = get.GaussianCurvature(obj)
            if (isnan(obj.V))
                error('Mesh3D object is not initialized with a mesh');
            end
            V1 = obj.V(obj.F(:,1),:);   % 1st vertices of faces
            V2 = obj.V(obj.F(:,2),:);   % 2nd vertices of faces
            V3 = obj.V(obj.F(:,3),:);   % 3rd vertices of faces

            % Head angles near the vertices
            V1Angles = GetVectorsAngle(V2 - V1, V3 - V1);
            V2Angles = GetVectorsAngle(V1 - V2, V3 - V2);
            V3Angles = GetVectorsAngle(V2 - V3, V1 - V3);

            % Calculating the sum of the angles around the vertices
            sum_angles=accumarray(obj.F(:), [V1Angles;V2Angles;V3Angles]);
            % Finding if the point is on the boundary
            V_adj_ex = sparse(obj.F,obj.F(:,[2,3,1]),ones(obj.NF,3),obj.NV,obj.NV);
            V_adj_ex_sym = V_adj_ex + V_adj_ex';
            B_adj = (V_adj_ex_sym == 1) & (V_adj_ex == 1);
            is_boundary = double(sum(B_adj)>0)';
            GaussianCurvature = ((2-is_boundary)*pi - sum_angles)./obj.V_Area;
        end
        function F_norm = get.F_norm(obj)
            if (isnan(obj.V))
                error('Mesh3D object is not initialized with a mesh');
            end
            % Positions of the vertices in each triangle
            F1 = obj.V(obj.F(:,1),:);
            F2 = obj.V(obj.F(:,2),:);
            F3 = obj.V(obj.F(:,3),:);
            
            F_norm = cross(F3-F2,F1-F3,2);
            F_norm = F_norm ./ repmat(vecnorm(F_norm,2,2), 1, 3);
        end
        function V_norm = get.V_norm(obj)
            if (isnan(obj.V))
                error('Mesh3D object is not initialized with a mesh');
            end
            X = obj.VF_adj;
            Y = obj.F_norm .* repmat(obj.F_Area,1,3);

            V_norm = X * Y;
            V_norm = V_norm ./ repmat(vecnorm(V_norm,2,2),1,3);
        end
        function PrependicularEdges = get.PrependicularEdges(obj)
            if (isnan(obj.V))
                error('Mesh3D object is not initialized with a mesh');
            end
            % Positions of the vertices in each triangle
            F1 = obj.V(obj.F(:,1),:);
            F2 = obj.V(obj.F(:,2),:);
            F3 = obj.V(obj.F(:,3),:);
            % The edge vectors - Edge 1 is the edge opposite of the first
            % vertex in the triangle
            Edge1 = F3-F2;
            Edge2 = F1-F3;
            Edge3 = F2-F1;
            % Findong the normal to the triangles
            Fnorm = cross(F3-F2,F1-F3,2);
            Fnorm = Fnorm ./ repmat(vecnorm(Fnorm,2,2), 1, 3);
            % Turning the edge vectors 90 degrees into the triangle while
            % keeping their length
            PrependicularEdges = cross(repmat(Fnorm,3,1), [Edge1; Edge2; Edge3], 2);
        end
        function Gradient = get.Gradient(obj)
            if (isnan(obj.V))
                error('Mesh3D object is not initialized with a mesh');
            end
            RotatedEdges = obj.PrependicularEdges;
            % Dividing by 2*Triangle_area
            RotatedEdges = RotatedEdges ./ (2 * repmat(obj.F_Area,3,3));

            % Creating the sparse gradient matrix
            row_mat = repmat(reshape((1:(3*obj.NF))',size(obj.F)),3,1);
            col_mat = repmat(reshape(obj.F,[3*obj.NF 1]),1,3);
            Gradient = sparse(row_mat, col_mat, RotatedEdges, 3*obj.NF, obj.NV);
        end
        function Divergence = get.Divergence(obj)
            if (isnan(obj.V))
                error('Mesh3D object is not initialized with a mesh');
            end
            inv_Gv = sparse(1:obj.NV, 1:obj.NV, 1 ./ obj.V_Area, obj.NV, obj.NV);
            Gf = sparse(1:(3*obj.NF), 1:(3*obj.NF), repmat(obj.F_Area,3,1)', 3*obj.NF, 3*obj.NF);

            Divergence = -inv_Gv*(obj.Gradient)'*Gf;
        end
        function Laplacian = get.Laplacian(obj)
            if (isnan(obj.V))
                error('Mesh3D object is not initialized with a mesh');
            end
            inv_Gv = sparse(1:obj.NV, 1:obj.NV, 1 ./ obj.V_Area, obj.NV, obj.NV);
            Gf = sparse(1:(3*obj.NF), 1:(3*obj.NF), repmat(obj.F_Area,3,1)', 3*obj.NF, 3*obj.NF);
            grad = obj.Gradient;
            Laplacian = inv_Gv*grad'*Gf*grad;
        end
        function MeanCurvature = get.MeanCurvature(obj)
            if (isnan(obj.V))
                error('Mesh3D object is not initialized with a mesh');
            end
            MeanCurvature = dot(obj.Laplacian*obj.V,obj.V_norm,2) / 2;
        end
        function CotangentWeights = get.CotangentWeights(obj)
            if (isnan(obj.V))
                error('Mesh3D object is not initialized with a mesh');
            end
            row_mat = repmat(reshape((1:(3*obj.NF))',size(obj.F)),3,1);
            col_mat = repmat(reshape(obj.F,[3*obj.NF 1]),1,3);
            E = sparse(row_mat, col_mat, obj.PrependicularEdges, 3*obj.NF, obj.NV);
            inv_Gf = sparse(1:(3*obj.NF), 1:(3*obj.NF), repmat(1 ./ obj.F_Area,3,1)', 3*obj.NF, 3*obj.NF);
            CotangentWeights = E' * inv_Gf * E / 4;
            %CotangentWeights(1:(3*obj.NF+1):end) = 0;
        end
        function F_adj = get.F_adj(obj)
            if (isnan(obj.V))
                error('Mesh3D object is not initialized with a mesh');
            end
            vf_adj = obj.VF_adj;
            F_adj = double((vf_adj' * vf_adj) == 2);
        end
        
        
        %%%%% Saving %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Save(obj, filename, filetype)
            if (nargin < 1)
                error('Please input a Mesh3D object');
            elseif (nargin == 1)
                error('Please input a filename');
            elseif (nargin == 2)
                filetype = lower(filename((find(filename == '.', 1, 'last')+1):end));
            end
            switch lower(filetype)
                case 'off'
                    obj.SaveOFFFile(filename);
                otherwise
                    error('Please input a filename with a supported extension (.off)');
            end
                    
        end
        
        %%%%% Viewing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [h] = ShowMesh(obj, color)
            if (~exist('obj','var'))
                error('Please input a Mesh3D object');
            end
            if (~exist('color','var'))
                color = 'blue';
            end
            h=patch('Faces',obj.F,'Vertices',obj.V,'FaceColor',color);
            axis('equal');
            cameratoolbar;
        end
        function [h] = PlotFunction(obj, values, onVertices)
            if (~exist('obj','var'))
                error('Please input a Mesh3D object');
            end
            if (~exist('values','var'))
                error('Please input function values, or use ShowMesh to display the mesh');
            end
            if (~ismatrix(values) || size(values,2) ~= 1)
                error('Please input function values as NF x 1 or as NV x 1 column vector');
            end
            if (~exist('onVertices','var'))
                onVertices = false;
            end
            
            if (~onVertices && size(values,1) == obj.NF)
                h=patch('Faces',obj.F,'Vertices',obj.V,'FaceVertexCData',values,'FaceColor','flat');
            elseif (onVertices && size(values,1) == obj.NV)
                h=patch('Faces',obj.F,'Vertices',obj.V,'FaceVertexCData',values,'FaceColor','interp');
            else
                error('Please input function values as NF x 1 or as NV x 1 column vector');
            end
            axis('equal');
            cameratoolbar;
            colorbar;
        end
        function [h] = PlotRemeshing(obj, values, onVertices)
            if (~exist('obj','var'))
                error('Please input a Mesh3D object');
            end
            if (~exist('values','var'))
                error('Please input function values, or use ShowMesh to display the mesh');
            end
            if (~ismatrix(values) || size(values,2) ~= 1)
                error('Please input function values as NF x 1 or as NV x 1 column vector');
            end
            if (~exist('onVertices','var'))
                onVertices = false;
            end
            
            colormap hsv;
            if (~onVertices && size(values,1) == obj.NF)
                h=patch('Faces',obj.F,'Vertices',obj.V,'FaceVertexCData',values, ...
                    'FaceColor','flat', 'LineStyle', 'none');
            elseif (onVertices && size(values,1) == obj.NV)
                h=patch('Faces',obj.F,'Vertices',obj.V,'FaceVertexCData',values, ...
                    'FaceColor','interp', 'LineStyle', 'none');
            else
                error('Please input function values as NF x 1 or as NV x 1 column vector');
            end
            axis('equal');
            cameratoolbar;
            colorbar;
        end
        function [h] = PlotVectorField(obj, vectors, onVertices)
            if (~exist('obj','var'))
                error('Please input a Mesh3D object');
            end
            if (~exist('vectors','var'))
                error('Please input the vectors to be shown');
            end
            if (~exist('onVertices','var'))
                onVertices = false;
            end
            if (~ismatrix(vectors) || size(vectors,2) ~= 3)
                if(all(size(vectors) == [obj.NF*3 1]) && onVertices == false)
                    vectors = reshape(vectors,size(obj.F));
                else
                    error('Please input the vectors as NF x 3 or as NV x 3 column vector');
                end
            end
            
            
            if (~onVertices && size(vectors,1) == obj.NF)
                FX = mean(reshape(obj.V(obj.F,1),size(obj.F)),2);
                FY = mean(reshape(obj.V(obj.F,2),size(obj.F)),2);
                FZ = mean(reshape(obj.V(obj.F,3),size(obj.F)),2);
                h=quiver3(FX,FY,FZ,vectors(:,1),vectors(:,2),vectors(:,3),'r');
            elseif (onVertices && size(vectors,1) == obj.NV)
                h=quiver3(obj.V(:,1),obj.V(:,2),obj.V(:,3),vectors(:,1),vectors(:,2),vectors(:,3),'r');
            else
                error('Please input vector values as NF x 3 or as NV x 3 column vector');
            end
            axis('equal');
            cameratoolbar;
            colorbar;
        end
        
        %%%%% Manipulating %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [mesh] = ExtraTriangles(obj)
        % This function creates a new mesh that has the same structure as
        % the received mesh, but cuts every edge in the middle to add new
        % triangles. The number of triangles will increase to 4 times its
        % previous number.
            F1 = obj.F(:,1);
            F2 = obj.F(:,2);
            F3 = obj.F(:,3);

            total_edges = [F1 F2; F2 F3; F3 F1];
            total_edges = sort(total_edges,2);
            edges = unique(total_edges,'rows');
            new_vertices = (obj.V(edges(:,1),:) + obj.V(edges(:,2),:))/2;

            spedges = sparse([edges(:,1);edges(:,2)], ...
                [edges(:,2);edges(:,1)], ...
                [1:size(edges,1),1:size(edges,1)]',obj.NV,obj.NV);

            V12 = full(spedges(sub2ind(size(spedges),F1,F2))) + obj.NV;
            V23 = full(spedges(sub2ind(size(spedges),F2,F3))) + obj.NV;
            V31 = full(spedges(sub2ind(size(spedges),F3,F1))) + obj.NV;

            new_faces = [F1 V12 V31; ...
                         V12 F2 V23; ...
                         V31 V23 F3; ...
                         V31 V12 V23];

            mesh = Mesh3D([obj.V;new_vertices],new_faces);
        end
    end
    
    methods (Access = private)
        function SaveOFFFile(obj, filename) % Assume filename includes OFF extension
            fileID = fopen(filename,'w');
            if (fileID < 0)
                error(['Could not open/create file ' filename]);
            end
            fprintf(fileID,'OFF\n');
            fprintf(fileID,'%d %d %d\n', [obj.NV obj.NF 0]);
            fprintf(fileID,'%f %f %f\n', obj.V');
            fprintf(fileID,'3 %d %d %d\n', obj.F' - 1);
            fclose(fileID);
        end
    end
    
    %%%%%%%% Static Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static)
        function obj = ReadOFFFile(filename)    % Load an OFF File
            if (nargin < 1)
                error('Please input a filename');
            end
            obj = Mesh3D;
            fileID = fopen(filename,'r');
            if (fileID < 0)
                error(['File ' filename ' not found!']);
            end
            tempLine = fgetl(fileID);
            lineNumber = 1;
            if (strcmpi(tempLine,'OFF') == 0)   % Header
                error(['Bad file format for ' filename]);
            end
            while ~feof(fileID) % Ignore empty lines and comments that start with '#'
                lineNumber = lineNumber + 1;
                tempLine = fgetl(fileID);
                tempLine = strtrim(tempLine);
                if (isempty(tempLine) || tempLine(1) == '#')
                    continue;
                end
                break;
            end
            sizes = sscanf(tempLine, '%d'); % Read NV, NF, ignore NE
            if (~ismatrix(sizes) || size(sizes,1)<2)
                error(['Bad file format for ' filename]);
            end
            nv = sizes(1);
            nf = sizes(2);
%             lineNumber = lineNumber + 1;
%             fclose(fileID);
%             
%             % Read the vertices' positions
%             obj.V = readmatrix(filename,'FileType','text', ...
%                 'Range',[num2str(lineNumber) ':' num2str(lineNumber + nv - 1)]);
%             obj.V = obj.V(:,1:3);   % Remove any irrelevant column
%             % Read the faces
%             obj.F = readmatrix(filename,'FileType','text', ...
%                 'Range',[num2str(lineNumber + nv - 1) ':' num2str(lineNumber + nv + nf - 1)]);
            obj.V = fscanf(fileID, '%f%f%f\n', [3 nv]);
            obj.V = obj.V';
            obj.F = fscanf(fileID, '%f%f%f\n', [4 nf]);
            obj.F = obj.F';
            %obj.F = obj.F(2:end,:); % Required because of wierd bug in readmatrix
            if (norm(obj.F(:,1)-3) ~= 0)    % Allow only triangular meshes
                error(['File ' filename ' is not a triangular mesh']);
            end
            obj.F = obj.F(:,2:4) + 1; % Indices should start from 1 in MATLAB
        end
    end
    
    
    
end

