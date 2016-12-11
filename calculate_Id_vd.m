%UNTITLED Summary of this function goes here
%   Detailed explanation goes hered
kbT    = 0.026;
a_0    = 1.42e-10; %Graphene lattice constant
h_bar  = (6.626e-34)/(2*pi);
q      = 1.6e-19;
w      = 1e-6; % How wide is the transistor?
Vd     = linspace(0,0.5,500);
mu     = 0;

a      = 3/2;
b      = sqrt(3)/2;
kmax_x = pi/(a);
kmax_y = 2*pi/(3*b);
kmin_y = pi / (3*b);
Ek = graphene_bandstructure();
length_of_Ek = size(Ek,1);
k_x = linspace(0, pi/(3/2), length_of_Ek);
k_x = k_x(2:length_of_Ek);
k_y_limit = linspace(kmax_y, kmin_y, length_of_Ek);
k_y_limit = k_y_limit(2:length_of_Ek);
k_step = pi/(a_0*size(Ek,1));

y_resolution = 30;

% Recall v=(1/h_bar)dE/dk
grad = diff(Ek)/k_step;
v = (q/h_bar).*grad;
Ek_cut = Ek(2:length_of_Ek,:);
x = linspace(0,1e-2,size(Ek_cut,1));
num_of_x_indicies = length_of_Ek - 1;
% TODO: Remove the hard-coded 6 in the future.
fermi_across_y = zeros(num_of_x_indicies, 6);

for x_index=1:num_of_x_indicies
    % TODO: Remove the hard-coded 6 in the future.
    Ek_y = zeros(y_resolution, 6);
    k_y = linspace(k_y_limit(x_index), -k_y_limit(x_index), y_resolution);
    for y_index=1:y_resolution
        Ek_y(y_index,:) = graphene_E_k(k_x(x_index), k_y(y_index));
    end
    fermi = 1./(1+exp((Ek_y-mu)./kbT)); % What is mu?
    fermi_across_y(x_index, :) = trapz(k_y, fermi);
end

integrand = (1/(4.*pi^2)).*fermi_across_y.*v;

% MUST USE X TO DEFINE LIMITS OF INTEGRATION
fk_vk = trapz(x,integrand); % Integrate to infinity
fk = trapz(x,fermi_across_y);
vth = nansum(fk_vk./fk);
disp(vth)
%Id=-q.*w.*fk_vk;

v_inj = vth.*(1-exp(-Vd./kbT))./(1+exp(-Vd./kbT));
Id=-q.*w.*v_inj;

plot(Vd,-Id);