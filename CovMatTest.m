clear;
A = [2 1 1;1 2 1;1 1 2] / 12;
V1 = [0 0 0];
V2 = [1 0 0];
V3 = [0 1 0];
S1 = norm(cross(V1-V2, V3-V2)) / 2;
V4 = [1 0 0];
V5 = [1 1 0.5];
V6 = [0 1 0];
S2 = norm(cross(V4-V5, V6-V5)) / 2;

T1 = [V1;V2;V3];
T2 = [V4;V5;V6];
center = mean([T1;T2]);
T1_u = T1 - center;
T2_u = T2 - center;

cov = S1 * T1_u' * A * T1_u + S2 * T2_u' * A * T2_u;

[x,y,z] = sphere;
R = 5 * cov * [x(:), y(:), z(:)]';
X = reshape(R(1,:),21,21) + center(1);
Y = reshape(R(2,:),21,21) + center(2);
Z = reshape(R(3,:),21,21) + center(3);

figure;
surf(X,Y,Z);
xlabel('x');
ylabel('y');
zlabel('z');
hold on;
T = [T1;T2];
scatter3(T(:,1),T(:,2),T(:,3),'filled','r');


[U, S, V] = svd(cov);
newS = S;
newS(1,1) = S(3,3) / S(1,1);
newS(2,2) = S(3,3) / S(2,2);
newS(3,3) = 1;
M = U * newS * V';

[x,y,z] = sphere;
R = 0.05 * metrics{i} * [x(:), y(:), z(:)]';
X = reshape(R(1,:),21,21) + center(1);
Y = reshape(R(2,:),21,21) + center(2);
Z = reshape(R(3,:),21,21) + center(3);

figure;
surf(X,Y,Z);
xlabel('x');
ylabel('y');
zlabel('z');
hold on;
T = [T1;T2];
scatter3(T(:,1),T(:,2),T(:,3),'filled','r');







