function [X,E] = graphene_bandstructure()
clear all;

Ep_C = 1.2057;
Ed_C = 24.1657;
pp_pi =  -3.2600;
pd_pi =  2.4000;
dd_pi =  3.6000;
dd_delta =  -7.4000;

theta = pi/2;
phi = pi/3;
l = sin(theta)*cos(phi);
m = sin(theta)*sin(phi);

H_onsite = [Ep_C/2 0 0 pp_pi pd_pi 0;
    0 Ed_C/2 0  pd_pi dd_pi 0;
    0 0 Ed_C/2 0 0 dd_delta;
    0 0 0 Ep_C/2 0 0;
    0 0 0 0 Ed_C/2 0;
    0 0 0 0 0 Ed_C/2;];

H_onsite = H_onsite + H_onsite';

H_off_nonzero = [pp_pi/2 l*(1-2*m*m)*pd_pi m*(1-2*l*l)*pd_pi;
    0 (l*l*dd_pi+m*m*dd_delta)/2 l*m*dd_pi+l*m*dd_delta;
    0 0 (m*m*dd_pi+l*l*dd_delta)/2;];

H_off_nonzero = H_off_nonzero + H_off_nonzero';
zero_block = zeros(3,3);
H_off = [zero_block H_off_nonzero; H_off_nonzero' zero_block];

r_1 = [3/2 sqrt(3)/2];
r_2 = [3/2 -sqrt(3)/2];

% This is for kx, the x part of the Gamma L-X plot
kmax = 2*pi;
Nt = 100;

for Nk = 1:Nt
    k = [0 1] * kmax * (Nk - 1)/(Nt - 1);
    p1 = exp(i * sum(k.*r_1));
    p2 = exp(-i * sum(k.*r_1));
    p3 = exp(i * sum(k.*r_2));
    p4 = exp(-i * sum(k.*r_2));
    H = H_onsite + (p1 + p2 + p3 + p4)*H_off;
    [V,D] = eig(H);
    eigst = sum(D);
    E(Nk,:) = sort(real(eigst));
    X(Nk) = kmax * (Nk-1)/(Nt-1);
end
dirac_point=(2*pi)/3;
h = plot(X,E,'b');
xlabel('k')
ylabel('E_k')
end


