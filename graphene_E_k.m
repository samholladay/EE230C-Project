function [E] = graphene_E_k(k_x, k_y)

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
H_off_pos_r2 = [zero_block H_off_pos_r2; zero_block zero_block];
% H_off_pos_r2 = H_off_pos_r1;

% Off diagonal elements for negative r1
H_off_neg_r1 = H_off_pos_r1';

% Off diagonal elements for negative r2
H_off_neg_r2 = H_off_pos_r2';

k = [k_x k_y];
%     k = [kmax_Y-K_y/2 -K_x/2] * (Nk - 1)/(Nt - 1);
p1 = exp( 1i * dot(k,r_1))*H_off_pos_r1;
p2 = exp(-1i * dot(k,r_1))*H_off_neg_r1;
p3 = exp( 1i * dot(k,r_2))*H_off_pos_r2;
p4 = exp(-1i * dot(k,r_2))*H_off_neg_r2;
H = H_onsite + p1 + p2 + p3 + p4;
[V,eigst] = eig(H, 'vector');
E = eigst;

end