%UNTITLED Summary of this function goes here
%   Detailed explanation goes hered
kbT    = 0.026; %eV
a_0    = 1.42e-10; %Graphene lattice constant
% h_bar  = (6.626e-34)/(2*pi);
h_bar  = 6.582e-16;
% q      = 1.6e-19;
q      = 1;
w      = 1e-6; % How wide is the transistor?
Vd     = linspace(0,0.5,500); %Volts
mu     = 0;

a      = 3/2;
b      = sqrt(3)/2;
kmax_x = pi/(a); 
kmax_y = 2*pi/(3*b);
kmin_y = pi / (3*b);
%Ek = graphene_bandstructure();
%length_of_Ek = size(Ek,1);
%k_x = linspace(0, kmax_x, length_of_Ek);
k_x = linspace(0, kmax_x, 1000);
k_y = linspace(kmax_y, kmin_y, 1000);
k_y_limit = k_y_limit(3:length_of_Ek);
k_step = pi/(a_0*1000);
Ek_x = zeros(1000,1);
for i =1:1000:
    Ek_x(i) = graphene_E_k
%k_x = k_x(3:length_of_Ek);
k_y_limit = linspace(kmax_y, kmin_y, length_of_Ek);
k_y_limit = k_y_limit(3:length_of_Ek);
k_step = pi/(a_0*size(Ek,1));

Ek = Ek(2:length_of_Ek,:);
grad = diff(Ek)/k_step;
v = (q/h_bar).*grad;
Ek_cut = Ek(2:length_of_Ek-1,:);

y_resolution = 100;

% Recall v=(1/h_bar)dE/dk
% Chopping off the first element of Ek because it is 0.
Ek = Ek(2:length_of_Ek,:);
grad = diff(Ek)/k_step;
v = (q/h_bar).*grad;
% v(1,:) = zeros(1,6);
Ek_cut = Ek(2:length_of_Ek-1,:);
x = linspace(0,1e-2,size(Ek_cut,1));
num_of_x_indicies = length_of_Ek - 2;

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
integrand2 = (1/(4.*pi^2)).*fermi_across_y;

% MUST USE X TO DEFINE LIMITS OF INTEGRATION
fk_vk = trapz(k_x,integrand); % Integrate to infinity
fk = trapz(k_x,integrand2);
vth = nansum(fk_vk./fk);
disp(vth)
%Id=-q.*w.*fk_vk;

Vd_max = 4;
Vd_step = 0.1;
Id_Vd = zeros(int32(Vd_max/Vd_step), 1);
index = 1;
fk_vk_neg_list = zeros(int32(Vd_max/Vd_step), 1);
fermi_neg_list = zeros(int32(Vd_max/Vd_step), 6);
%disp(Vd_step)
for Vd = 0:Vd_step:Vd_max
    fermi_across_y_negative = zeros(num_of_x_indicies, 6);
    %disp(num_of_x_indicies)
    for x_index = 1:num_of_x_indicies
        Ek_y = zeros(y_resolution, 6);
        k_y = linspace(k_y_limit(x_index), -k_y_limit(x_index), y_resolution);
        for y_index=1:y_resolution
            Ek_y(y_index,:) = graphene_E_k(-k_x(x_index), k_y(y_index));
        end
        fermi = 1./(1+exp((Ek_y-(mu - Vd))./kbT)); % What is mu?
        %disp(index)
        %disp(x_index)
        %fermi_neg_list(index, x_index)=fermi;
        fermi_across_y_negative(x_index, :) = trapz(k_y, fermi);
    end

    integrand_neg = (1/(4*pi^2)).*fermi_across_y_negative.*v;
    integrand2_neg = (1/(4*pi^2)).*fermi_across_y_negative;
    fk_vk_neg = trapz(k_x, integrand_neg);
    fk_neg = trapz(k_x, integrand2_neg);  
    Id_Vd(index) = -q*w*(sum(fk_vk - fk_vk_neg));
    fk_vk_neg_list(index) = sum(fk_vk_neg);
    index = index + 1;
    
end

% % Back-scattering:


% v_inj = vth.*(1-exp(-Vd./kbT))./(1+exp(-Vd./kbT));
% Id=-q.*w.*v_inj.*sum(fk);

% plot(Vd,Id);
plot(Id_Vd);