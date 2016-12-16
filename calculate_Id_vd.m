%UNTITLED Summary of this function goes here
%   Detailed explanation goes hered
h_bar  = (6.626e-34)/(2*pi);
% h_bar        = 6.582e-16;
q            = 1.6e-19;
% q            = 1;
kbT          = 0.026; %eV
a_0          = 1.42e-10; %Graphene lattice constant

w            = 1e-6; % How wide is the transistor?
Vd           = linspace(0,4,30); %Volts
mu           = -0.3;
y_resolution = 100;
num_bands    = 6;

a            = 3/2;
b            = sqrt(3)/2;
kmax_x       = pi/(a); 
kmax_y       = 2*pi/(3*b);
kmin_y       = pi / (3*b);

% Ek = graphene_bandstructure();
length_of_Ek = size(Ek,1);

k_x = linspace(0, kmax_x, length_of_Ek);
k_x = k_x(3:length_of_Ek);
k_y_limit = linspace(kmax_y, kmin_y, length_of_Ek);
k_y_limit = k_y_limit(3:length_of_Ek);

k_step = pi/(a_0*size(Ek,1));
% Recall v=(1/h_bar)dE/dk
% Chopping off the first element of Ek because it is 0.
Ek = Ek(2:length_of_Ek,:);
grad = diff(Ek)/k_step;
v = (q/h_bar).*grad;


Ek_cut = Ek(2:length_of_Ek-1,:);
num_of_x_indicies = size(Ek_cut,1);

fermi_across_y       = zeros(num_of_x_indicies, num_bands);
fermi_across_y_holes = zeros(num_of_x_indicies, num_bands);

for x_index=1:num_of_x_indicies
    Ek_y = zeros(y_resolution, num_bands);
    k_y = linspace(k_y_limit(x_index), -k_y_limit(x_index), y_resolution);
    for y_index=1:y_resolution
        Ek_y(y_index,:) = graphene_E_k(k_x(x_index), k_y(y_index));
    end
    fermi = 1./(1+exp((Ek_y-mu)./kbT)); % What is mu?
    fermi_holes = exp((Ek_y+mu)./kbT)./(1+exp((Ek_y+mu)./kbT));
    
    fermi_across_y(x_index, :)       = trapz(k_y, fermi);
    fermi_across_y_holes(x_index, :) = trapz(k_y, fermi_holes);
end

integrand  = (1/(4.*pi^2)).*fermi_across_y.*v;
integrand2 = (1/(4.*pi^2)).*fermi_across_y;
integrand_holes = (1/(4.*pi^2)).*fermi_across_y_holes.*v;

% MUST USE X TO DEFINE LIMITS OF INTEGRATION
fk_vk = trapz(k_x,integrand); % Integrate to infinity
fk    = trapz(k_x,integrand2);
fk_vk_holes = trapz(k_x, integrand_holes);
vth = nansum(fk_vk./fk);
disp('Vth is: ');
disp(vth);
%Id=-q.*w.*fk_vk;

Id_Vd_electrons = zeros(length(Vd), 1);
Id_Vd_holes     = zeros(length(Vd), 1);
fk_vk_neg_list  = zeros(length(Vd), 1);
fermi_neg_list  = zeros(length(Vd), num_bands);

for index = 1:length(Vd)
    fermi_across_y_negative  = zeros(num_of_x_indicies, num_bands);
    fermi_across_y_neg_holes = zeros(num_of_x_indicies, num_bands);
    for x_index = 1:num_of_x_indicies
        Ek_y = zeros(y_resolution, num_bands);
        k_y = linspace(k_y_limit(x_index), -k_y_limit(x_index), y_resolution);
        for y_index=1:y_resolution
            Ek_y(y_index,:) = graphene_E_k(-k_x(x_index), k_y(y_index));
        end
        fermi = 1./(1+exp((Ek_y-(mu - Vd(index)))./kbT)); % What is mu?
        fermi_neg_holes =  exp((Ek_y+(mu - Vd(index)))./kbT)./(1+exp((Ek_y+(mu - Vd(index)))./kbT));
        %fermi_neg_list(index, x_index)=fermi;
        fermi_across_y_negative(x_index, :)  = trapz(k_y, fermi);
        fermi_across_y_neg_holes(x_index, :) = trapz(k_y, fermi_neg_holes);
    end

    integrand_neg       = (1/(4*pi^2)).*fermi_across_y_negative.*v;
    integrand2_neg      = (1/(4*pi^2)).*fermi_across_y_negative;
    integrand_neg_holes = (1/(4*pi^2)).*fermi_across_y_neg_holes.*v;
    
    fk_vk_neg       = trapz(k_x, integrand_neg);
    fk_vk_neg_holes = trapz(k_x, integrand_neg_holes);
    fk_neg          = trapz(k_x, integrand2_neg);  
    
    Id_Vd_electrons(index) = -q*w*(sum(fk_vk - fk_vk_neg));
    Id_Vd_holes(index)     = -q*w*(sum(fk_vk_holes - fk_vk_neg_holes));
    fk_vk_neg_list(index)  = sum(fk_vk_neg);
    
end

% % Back-scattering:


% v_inj = vth.*(1-exp(-Vd./kbT))./(1+exp(-Vd./kbT));
% Id=-q.*w.*v_inj.*sum(fk);

% plot(Vd,Id);
plot(Vd, Id_Vd_electrons);