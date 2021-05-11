function [theta] = GetVectorsAngle(u, v)
% theta = GETVECTORSANGLE(U, V)
% Returns the angle (in radians) between two vectors.
%
% U and V are N x 3 matrices whose rows represent vectors
    theta = atan2(vecnorm(cross(u,v),2,2),dot(u,v,2));
end