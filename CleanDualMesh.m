function [V, F, edges] = CleanDualMesh(V, F)
%CLEANDUALMESH Summary of this function goes here
%   Detailed explanation goes here
[bad_faces_indices,~] = find(isnan(F(:,3)));
bad_faces = sortrows(sort(F(bad_faces_indices, 1:2),2));
new_F = F;
deleted_indices = false(size(V,1), 1);

% Mark some vertices to be deleted and compute updated vertices
for i=1:size(bad_faces,1)
    [row, ~] = find(F == bad_faces(i,1));
    [row2, ~] = find(F == bad_faces(i,2));
    
    if (size(row,1) > 3 && size(row2,1) > 3) || ...
            isnan(bad_faces(i,2))
        new_F(bad_faces_indices(i), :) = NaN; % Only delete the face
    elseif size(row,1) > 3
        V(bad_faces(i,1),:) = (V(bad_faces(i,1),:) + V(bad_faces(i,2),:)) ./ 2;
        new_F(F == bad_faces(i,2)) = NaN;
        deleted_indices(bad_faces(i,2)) = true;
    elseif size(row2,1) > 3
        V(bad_faces(i,2),:) = (V(bad_faces(i,1),:) + V(bad_faces(i,2),:)) ./ 2;
        new_F(F == bad_faces(i,1)) = NaN;
        deleted_indices(bad_faces(i,1)) = true;
    else
        new_F(F == bad_faces(i,1) | F == bad_faces(i,2)) = NaN;
        deleted_indices(bad_faces(i,1)) = true;
        deleted_indices(bad_faces(i,2)) = true;
    end
end
% "Compress" new_F
keep_faces = true(size(F,1), 1);
edges_num = 0;
for i=1:size(new_F,1)
    curr_row = new_F(i, ~isnan(new_F(i,:)));
    curr_size = size(curr_row,2);
    edges_num = edges_num + curr_size;
    if curr_size == 0
        keep_faces(i) = false;
        continue
    end
    new_F(i,1:size(curr_row,2)) = curr_row;
    new_F(i,(size(curr_row,2)+1):end) = NaN;
end
F = new_F(keep_faces,:);

% Update vertices
for i=flip(find(deleted_indices)')
    F(F>i) = F(F>i) - 1;
end
V(deleted_indices, :) = [];

% Create edges matrix
edges = zeros(edges_num, 2);
edges_num = 1;
flag = false;
for i=1:size(F,1)
    curr_row = F(i, ~isnan(F(i,:)));
    curr_size = size(curr_row,2);
    if curr_size < 3
        flag = true;
        break
    end
    edges(edges_num:(edges_num+curr_size-1), 1) = curr_row';
    edges(edges_num:(edges_num+curr_size-1), 2) = curr_row([2:curr_size 1])';
    edges_num = edges_num + curr_size;
end

if flag
    [V, F, edges] = CleanDualMesh(V, F);
    return
end

edges = unique(sort(edges,2), 'rows'); % Removing edge repetitions

end

