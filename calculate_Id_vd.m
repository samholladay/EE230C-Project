%UNTITLED Summary of this function goes here
%   Detailed explanation goes hered
% h_bar  = (6.626e-34)/(2*pi); % SI
h_bar        = 6.582e-16; %eV
% q            = 1.6e-19;
q            = 1;
q_si         = 1.6e-19; 
kbT          = 0.026; %eV
a_0          = 1.42e-10; %Graphene lattice constant

w            = 1e-6; % How wide is the transistor?
Vd           = linspace(0,5.5,50); %Volts
mu           = 0.3;
y_resolution = 100;
num_bands    = 6;

a            = 3/2;
b            = sqrt(3)/2;
kmax_x       = pi/(a); 
kmax_y       = 2*pi/(3*b);
kmin_y       = pi / (3*b);

Ek = graphene_bandstructure();
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

E = zeros(num_of_x_indicies, y_resolution, num_bands);
E1 = zeros(num_of_x_indicies, y_resolution);
E2 = zeros(num_of_x_indicies, y_resolution);
E3 = zeros(num_of_x_indicies, y_resolution);
E4 = zeros(num_of_x_indicies, y_resolution);
E5 = zeros(num_of_x_indicies, y_resolution);
E6 = zeros(num_of_x_indicies, y_resolution);
for x_index = 1:num_of_x_indicies
    Ek_y = zeros(y_resolution, num_bands);
    k_y = linspace(k_y_limit(x_index), -k_y_limit(x_index), y_resolution);
    for y_index=1:y_resolution
        temp_E = graphene_E_k(-k_x(x_index), k_y(y_index));
        E(x_index, y_index, :) = temp_E;
        E1(x_index, y_index) = temp_E(1);
        E2(x_index, y_index) = temp_E(2);
        E3(x_index, y_index) = temp_E(3);
        E4(x_index, y_index) = temp_E(4);
        E5(x_index, y_index) = temp_E(5);
        E6(x_index, y_index) = temp_E(6);
    end
end

[Vx1, Vy] = gradient(E1);
[Vx2, Vy] = gradient(E2);
[Vx3, Vy] = gradient(E3);
[Vx4, Vy] = gradient(E4);
[Vx5, Vy] = gradient(E5);
[Vx6, Vy] = gradient(E6);
Vx1 = (q/h_bar) .* Vx1;
Vx2 = (q/h_bar) .* Vx2;
Vx3 = (q/h_bar) .* Vx3;
Vx4 = (q/h_bar) .* Vx4;
Vx5 = (q/h_bar) .* Vx5;
Vx6 = (q/h_bar) .* Vx6;

fk_vk_across_y       = zeros(num_of_x_indicies, num_bands);
fk_vk_across_y_holes = zeros(num_of_x_indicies, num_bands);
fermi_across_y       = zeros(num_of_x_indicies, num_bands);

for x_index=1:num_of_x_indicies
    Ek_y = E(x_index, :,:);
    fermi       = 1./(1+exp((Ek_y-mu)./kbT)); % What is mu?
    fermi_holes = 1./(1+exp((Ek_y+mu)./kbT));
    V_x = [Vx1(x_index,:)'; Vx2(x_index,:)'; Vx3(x_index,:)'; Vx4(x_index,:)'; Vx5(x_index,:)'; Vx6(x_index,:)'];
%     fermi_holes = 1 - fermi_holes;
    fermi_across_y(x_index, :)       = trapz(k_y, fermi);
    fk_vk_across_y(x_index, :)       = trapz(k_y, fermi.*V_x);
    fk_vk_across_y_holes(x_index, :) = trapz(k_y, fermi_holes.*V_x);
end

% integrand  = (1/(4.*pi^2)).*fk_vk_across_y.*v;
% integrand2 = (1/(4.*pi^2)).*fk_vk_across_y;
% integrand_holes = (1/(4.*pi^2)).*fk_vk_across_y_holes.*v;
integrand  = (1/(4.*pi^2)).*fk_vk_across_y;
integrand2 = (1/(4.*pi^2)).*fermi_across_y;
integrand_holes = (1/(4.*pi^2)).*fk_vk_across_y_holes;

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
        Ek_y = E(x_index, :,:);
        fermi = 1./(1+exp((Ek_y-(mu - Vd(index)))./kbT)); % What is mu?
        fermi_neg_holes = 1./(1+exp((Ek_y+(mu - Vd(index)))./kbT));
%         fermi_neg_holes = 1 - fermi_neg_holes;
        %fermi_neg_list(index, x_index)=fermi;
        fermi_across_y_negative(x_index, :)  = trapz(k_y, fermi);
        fermi_across_y_neg_holes(x_index, :) = trapz(k_y, fermi_neg_holes);
    end

    integrand_neg       = (1/(4*pi^2)).*fermi_across_y_negative.*v;
    integrand_neg_holes = (1/(4*pi^2)).*fermi_across_y_neg_holes.*v;
    integrantd_neg_new  = (1/(4*pi^2)).*fermi_across_y_negative.*v;
    
    fk_vk_neg       = trapz(k_x, integrand_neg);
    fk_vk_neg_holes = trapz(k_x, integrand_neg_holes);
    
    Id_Vd_electrons(index) = -q_si*w*(sum(fk_vk - fk_vk_neg));
    Id_Vd_holes(index)     =  q_si*w*(sum(fk_vk_holes - fk_vk_neg_holes));
    fk_vk_neg_list(index)  = sum(fk_vk_neg);
    
%     fk_vk_neg_new(index,:) = trapz(k_x, trapz(1/1+exp((E-(mu + Vd(index)))/kbT)));
    
end

% % Back-scattering:

figure();
plot(Vd, Id_Vd_electrons+Id_Vd_holes);
title(['Total Current']);
xlabel(['Drain Voltage']);
ylabel(['Drain Current (\mu m / \mu m)']);
figure();
plot(Vd, Id_Vd_electrons)
title(['Electrons']);
xlabel(['Drain Voltage']);
ylabel(['Drain Current (\mu m / \mu m)']);
figure();
plot(Vd, Id_Vd_holes)
title(['Holes'])
xlabel(['Drain Voltage']);
ylabel(['Drain Current (\mu m / \mu m)']);

% v_inj = vth.*(1-exp(-Vd./kbT))./(1+exp(-Vd./kbT));
% Id=-q.*w.*v_inj.*sum(fk);
%%
% subplot(3, 1, 1)
% plot(Vd, Id_Vd_electrons+Id_Vd_holes);
% title(['Total Current']);
% xlabel(['Drain Voltage']);
% ylabel(['Drain Current (\mu m / \mu m)']);
% subplot(3, 1, 2)
% plot(Vd, Id_Vd_electrons)
% title(['Electrons']);
% xlabel(['Drain Voltage']);
% ylabel(['Drain Current (\mu m / \mu m)']);
% subplot(3, 1, 3)
% plot(Vd, Id_Vd_holes)
% title(['Holes'])
% xlabel(['Drain Voltage']);
% ylabel(['Drain Current (\mu m / \mu m)']);

figure();
plot(Vd, Id_Vd_electrons+Id_Vd_holes);
title(['Total Current']);
xlabel(['Drain Voltage']);
ylabel(['Drain Current (\mu m / \mu m)']);
figure();
plot(Vd, Id_Vd_electrons)
title(['Electrons']);
xlabel(['Drain Voltage']);
ylabel(['Drain Current (\mu m / \mu m)']);
figure();
plot(Vd, Id_Vd_holes)
title(['Holes'])
xlabel(['Drain Voltage']);
ylabel(['Drain Current (\mu m / \mu m)']);