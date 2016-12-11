function [Emk] = graphene_bandstructure()
%Calculates and plots band structure for graphene
%returns bands along kx direction

Ep_C = 1.2057;
Ed_C = 24.1657;
pp_pi =  -3.2600;
pd_pi =  2.4000;
dd_pi =  3.6000;
dd_delta =  -7.4000;

r_1 = [3/2  sqrt(3)/2];
r_2 = [3/2 -sqrt(3)/2];
a = 3/2;
b = sqrt(3)/2;
kmax_Y = 2*pi/(3*b);
kmax_X = pi/a;
K_point_y = pi / (3*b);
K_point_x = pi/a;

H_onsite = [Ep_C/2 0      0      pp_pi  pd_pi  0;
            0      Ed_C/2 0      pd_pi  dd_pi  0;
            0      0      Ed_C/2 0      0      dd_delta;
            0      0      0      Ep_C/2 0      0;
            0      0      0      0      Ed_C/2 0;
            0      0      0      0      0      Ed_C/2;];
H_onsite = H_onsite + H_onsite';
% 
% H_off_nonzero = [pp_pi/2 l*(1-2*m*m)*pd_pi          m*(1-2*l*l)*pd_pi;
%                  0       (l*l*dd_pi+m*m*dd_delta)/2 l*m*dd_pi+l*m*dd_delta;
%                  0       0                          (m*m*dd_pi+l*l*dd_delta)/2;];

zero_block = zeros(3,3);

% Off diagonal elements for positive r1
theta = pi/3;
l = cos(theta);
m = sin(theta);
% H_off_pos_r1 = [pp_pi/2 l*pd_pi                    m*pd_pi;
%                 0       (l*l*dd_pi+m*m*dd_delta)/2 l*m*dd_pi - l*m*dd_delta;
%                 0       0                          (m*m*dd_pi+l*l*dd_delta)/2];
% H_off_pos_r1 = H_off_pos_r1 + H_off_pos_r1';
H_off_pos_r1 = [pp_pi    l*pd_pi                  m*pd_pi;
                -l*pd_pi (l*l*dd_pi+m*m*dd_delta) l*m*dd_pi - l*m*dd_delta;
                -m*pd_pi  l*m*dd_pi - l*m*dd_delta (m*m*dd_pi+l*l*dd_delta)];

H_off_pos_r1 = [zero_block H_off_pos_r1; zero_block zero_block];

% Off diagonal elements for positive r2
theta = -pi/3;
l = cos(theta);
m = sin(theta);
H_off_pos_r2 = [pp_pi    l*pd_pi                  m*pd_pi;
                -l*pd_pi (l*l*dd_pi+m*m*dd_delta) l*m*dd_pi - l*m*dd_delta;
                -m*pd_pi  l*m*dd_pi - l*m*dd_delta (m*m*dd_pi+l*l*dd_delta)];
% H_off_pos_r2 = H_off_pos_r2 + H_off_pos_r2';
H_off_pos_r2 = [zero_block H_off_pos_r2; zero_block zero_block];
% H_off_pos_r2 = H_off_pos_r1;

% Off diagonal elements for negative r1
% theta = 2*pi/3;
% l = cos(theta);
% m = sin(theta);
% H_off_neg_r1 = [pp_pi/2 l*pd_pi                    m*pd_pi;
%                 0       (l*l*dd_pi+m*m*dd_delta)/2 l*m*dd_pi + l*m*dd_delta;
%                 0       0                          (m*m*dd_pi+l*l*dd_delta)/2];
% H_off_neg_r1 = H_off_neg_r1 + H_off_neg_r1';
% H_off_neg_r1 = [zero_block zero_block; H_off_neg_r1' zero_block];
H_off_neg_r1 = H_off_pos_r1';

% Off diagonal elements for negative r2
% theta = 4*pi/3;
% l = cos(theta);
% m = sin(theta);
% H_off_neg_r2 = [pp_pi/2 l*pd_pi                    m*pd_pi;
%                 0       (l*l*dd_pi+m*m*dd_delta)/2 l*m*dd_pi + l*m*dd_delta;
%                 0       0                          (m*m*dd_pi+l*l*dd_delta)/2];
% H_off_neg_r2 = H_off_neg_r2 + H_off_neg_r2';
% H_off_neg_r2 = [zero_block zero_block; H_off_neg_r2' zero_block];
H_off_neg_r2 = H_off_pos_r2';

% H_off = [zero_block H_off_nonzero; H_off_nonzero' zero_block];

% H_off_neg_r1 = H_off_pos_r1';
% H_off_neg_r2 = H_off_pos_r2';

Nt = 1000;
E = zeros(3 * Nt, 6);
Emk = zeros(Nt, 6);
index = 0;

for Nk = 1:Nt
    k = [kmax_X 0] * (Nk - 1)/(Nt - 1);
%     k = [kmax_Y-K_y/2 -K_x/2] * (Nk - 1)/(Nt - 1);
    p1 = exp( 1i * dot(k,r_1))*H_off_pos_r1;
    p2 = exp(-1i * dot(k,r_1))*H_off_neg_r1;
    p3 = exp( 1i * dot(k,r_2))*H_off_pos_r2;
    p4 = exp(-1i * dot(k,r_2))*H_off_neg_r2;
    H = H_onsite + p1 + p2 + p3 + p4;
    [V,eigst] = eig(H, 'vector');
    E(index + (Nt - Nk) + 1,:) = eigst;
    Emk(index + Nk + 1, :) = eigst;
end

index = Nt;
for Nk = 1:Nt
    k = [K_point_x K_point_y] * (Nk - 1)/(Nt - 1);
    p1 = exp( 1i * dot(k,r_1))*H_off_pos_r1;
    p2 = exp(-1i * dot(k,r_1))*H_off_neg_r1;
    p3 = exp( 1i * dot(k,r_2))*H_off_pos_r2;
    p4 = exp(-1i * dot(k,r_2))*H_off_neg_r2;
    H = H_onsite + p1 + p2 + p3 + p4;
    [V,eigst] = eig(H, 'vector');
    E(index + Nk,:) = eigst;
end

index = 2*Nt;

for Nk = 1:Nt
    k = ([0  K_point_y] * (Nk - 1)/(Nt - 1)) + [K_point_x 0];
    %disp(k)
    p1 = exp( 1i * dot(k,r_1))*H_off_pos_r1;
    p2 = exp(-1i * dot(k,r_1))*H_off_neg_r1;
    p3 = exp( 1i * dot(k,r_2))*H_off_pos_r2;
    p4 = exp(-1i * dot(k,r_2))*H_off_neg_r2;
    H = H_onsite + p1 + p2 + p3 + p4;
    [V,eigst] = eig(H, 'vector');
    E(index + (Nt - Nk) + 1,:) = eigst;
end

h = plot(E,'b');
% axis([0.5 2 -3 6])
end