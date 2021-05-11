function [theta] = GetVectorsAngleD(u, v)
% theta = GETVECTORSANGLED(U, V)
% Returns the angle (in degrees) between two vectors.
%
% U and V are N x 3 matrices whose rows represent vectors
    theta = atan2d(vecnorm(cross(u,v),2,2),dot(u,v,2));
end