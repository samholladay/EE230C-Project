%UNTITLED Summary of this function goes here
%   Detailed explanation goes hered
clear all
% h_bar  = (6.626e-34)/(2*pi); % SI
h_bar        = 6.582e-16; %eV
% q            = 1.6e-19;
q            = 1;
q_si         = 1.6e-19; 
kbT          = 0.026; %eV
a_0          = 1.42e-10; %Graphene lattice constant

w            = 1e-6; % How wide is the transistor?
Vd           = linspace(0,1,40); %Volts
vgs           = 0.3;
y_resolution = 400;
num_bands    = 6;
delta = 0.01;

a            = 3/2;
b            = sqrt(3)/2;
kmax_x       = pi/(a); 
kmax_y       = 2*pi/(3*b);
kmin_y       = pi / (3*b);

% Ek = graphene_bandstructure();
x_resolution = 50;

k_x = linspace(0, kmax_x, x_resolution);
k_y_limit = linspace(kmax_y, kmin_y, x_resolution);

% Recall v=(1/h_bar)dE/dk

E = zeros(x_resolution, y_resolution, num_bands);
E1 = zeros(x_resolution, y_resolution);
E2 = zeros(x_resolution, y_resolution);
E3 = zeros(x_resolution, y_resolution);
E4 = zeros(x_resolution, y_resolution);
E5 = zeros(x_resolution, y_resolution);
E6 = zeros(x_resolution, y_resolution);

Vx1 = zeros(x_resolution, y_resolution);
Vx2 = zeros(x_resolution, y_resolution);
Vx3 = zeros(x_resolution, y_resolution);
Vx4 = zeros(x_resolution, y_resolution);
Vx5 = zeros(x_resolution, y_resolution);
Vx6 = zeros(x_resolution, y_resolution);

for x_index = 1:x_resolution
    Ek_y = zeros(y_resolution, num_bands);
    k_y = linspace(-k_y_limit(x_index), k_y_limit(x_index), y_resolution);
    for y_index=1:y_resolution
        temp_E = graphene_E_k(-k_x(x_index), k_y(y_index));
        temp_V = a_0 * (q/h_bar) .* k_to_v(-k_x(x_index), k_y(y_index), 0.1);
        E(x_index, y_index, :) = temp_E;
        E1(x_index, y_index) = temp_E(1);
        E2(x_index, y_index) = temp_E(2);
        E3(x_index, y_index) = temp_E(3);
        E4(x_index, y_index) = temp_E(4);
        E5(x_index, y_index) = temp_E(5);
        E6(x_index, y_index) = temp_E(6);
        Vx1(x_index, y_index) = temp_V(1);
        Vx2(x_index, y_index) = temp_V(2);
        Vx3(x_index, y_index) = temp_V(3);
        Vx4(x_index, y_index) = temp_V(4);
        Vx5(x_index, y_index) = temp_V(5);
        Vx6(x_index, y_index) = temp_V(6);
    end
end

% [Vx1, Vy] = gradient(E1);
% [Vx2, Vy] = gradient(E2);
% [Vx3, Vy] = gradient(E3);
% [Vx4, Vy] = gradient(E4);
% [Vx5, Vy] = gradient(E5);
% [Vx6, Vy] = gradient(E6);
% Vx1 = (q/(h_bar*k_step)) .* Vx1;
% Vx2 = (q/(h_bar*k_step)) .* Vx2;
% Vx3 = (q/(h_bar*k_step)) .* Vx3;
% Vx4 = (q/(h_bar*k_step)) .* Vx4;
% Vx5 = (q/(h_bar*k_step)) .* Vx5;
% Vx6 = (q/(h_bar*k_step)) .* Vx6;


fk_vk_across_y       = zeros(x_resolution, num_bands);
fk_vk_across_y_holes = zeros(x_resolution, num_bands);
fermi_across_y       = zeros(x_resolution, num_bands);

mu = channel_sc_potential(E, x_resolution, y_resolution, vgs, 0, -0.5*vgs);

for x_index=1:x_resolution
% for x_index=1:2
    Ek_y = [E1(x_index,:)' E2(x_index,:)' E3(x_index,:)' E4(x_index,:)' E5(x_index,:)' E6(x_index,:)'];
    k_y = linspace(-k_y_limit(x_index), k_y_limit(x_index), y_resolution);
    
    fermi       = 1./(1+exp((Ek_y-mu)./kbT)); % What is mu?
    fermi_holes = 1./(1+exp((Ek_y+mu)./kbT));
    V_x = [Vx1(x_index,:)' Vx2(x_index,:)' Vx3(x_index,:)' Vx4(x_index,:)' Vx5(x_index,:)' Vx6(x_index,:)'];
%     fermi_holes = 1 - fermi_holes;
    fermi_across_y(x_index, :)       = trapz(k_y/a_0, fermi);
    fk_vk_across_y(x_index, :)       = trapz(k_y/a_0, fermi.*V_x);
    fk_vk_across_y_holes(x_index, :) = trapz(k_y/a_0, fermi_holes.*V_x);
end

% integrand  = (1/(4.*pi^2)).*fk_vk_across_y.*v;
% integrand2 = (1/(4.*pi^2)).*fk_vk_across_y;
% integrand_holes = (1/(4.*pi^2)).*fk_vk_across_y_holes.*v;
integrand  = (1/(4.*pi^2)).*fk_vk_across_y;
integrand2 = (1/(4.*pi^2)).*fermi_across_y;
integrand_holes = (1/(4.*pi^2)).*fk_vk_across_y_holes;

% MUST USE X TO DEFINE LIMITS OF INTEGRATION
fk_vk = trapz(k_x/a_0,integrand); % Integrate to infinity
fk    = trapz(k_x/a_0,integrand2);
fk_vk_holes = trapz(k_x/a_0, integrand_holes);
vth = nansum(fk_vk./fk);
disp('Vth is: ');
disp(vth);
%Id=-q.*w.*fk_vk;

%%
Id_Vd_electrons = zeros(length(Vd), 1);
Id_Vd_holes     = zeros(length(Vd), 1);
fk_vk_neg_list  = zeros(length(Vd), 1);
fermi_neg_list  = zeros(length(Vd), num_bands);
vinj_Vd         = zeros(length(Vd), 1);

for index = 1:length(Vd)
    fk_vk_across_y       = zeros(x_resolution, num_bands);
    fk_vk_across_y_holes = zeros(x_resolution, num_bands);
    
    fk_vk_across_y_neg       = zeros(x_resolution, num_bands);
    fermi_neg_across_y = zeros(x_resolution, num_bands);
    fermi_forward_across_y = zeros(x_resolution, num_bands);
    fk_vk_across_y_neg_holes = zeros(x_resolution, num_bands);
    guess_mu = mu;
    for x_index = 1:x_resolution
        mu = channel_sc_potential(E, x_resolution, y_resolution, vgs, Vd(index), guess_mu);
        guess_mu = mu;
        k_y = linspace(-k_y_limit(x_index), k_y_limit(x_index), y_resolution);
        Ek_y = [E1(x_index,:)' E2(x_index,:)' E3(x_index,:)' E4(x_index,:)' E5(x_index,:)' E6(x_index,:)'];
        V_x = [Vx1(x_index,:)' Vx2(x_index,:)' Vx3(x_index,:)' Vx4(x_index,:)' Vx5(x_index,:)' Vx6(x_index,:)'];
        
        fermi_pos   = 1./(1+exp((Ek_y-mu)./kbT)); % What is mu?
        fermi_holes = 1./(1+exp((Ek_y+mu)./kbT));
    %     fermi_holes = 1 - fermi_holes;
        fk_vk_across_y(x_index, :)       = trapz(k_y/a_0, fermi_pos.*V_x);
        fermi_forward_across_y(x_index,:)= trapz(k_y/a_0, fermi_pos);
        fk_vk_across_y_holes(x_index, :) = trapz(k_y/a_0, fermi_holes.*V_x);
        
        fermi_neg       = 1./(1+exp((Ek_y-(mu - Vd(index)))./kbT)); % What is mu?
        fermi_neg_holes = 1./(1+exp((Ek_y+(mu - Vd(index)))./kbT));
%         fermi_neg_holes = 1 - fermi_neg_holes;
        %fermi_neg_list(index, x_index)=fermi;
        
        fk_vk_across_y_neg(x_index, :)  = trapz(k_y/a_0, fermi_neg.*V_x);
        fermi_neg_across_y(x_index, :)       = trapz(k_y/a_0, fermi_neg);
        fk_vk_across_y_neg_holes(x_index, :) = trapz(k_y/a_0, fermi_neg_holes.*V_x);
    end
    
    integrand  = (1/(4.*pi^2)).*fk_vk_across_y;
    integrand_holes = (1/(4.*pi^2)).*fk_vk_across_y_holes;
    integrand_fermi_forward = (1/(4.*pi^2)).*fermi_forward_across_y;

    % MUST USE X TO DEFINE LIMITS OF INTEGRATION
    fk_vk = trapz(k_x/a_0,integrand); % Integrate to infinity
    fk_forward = trapz(k_x/a_0,integrand_fermi_forward);
    fk_vk_holes = trapz(k_x/a_0, integrand_holes);

    integrand_neg       = (1/(4*pi^2)).*fk_vk_across_y_neg;
    integrand_fermi_neg = (1/(4.*pi^2)).*fermi_neg_across_y;
    integrand_neg_holes = (1/(4*pi^2)).*fk_vk_across_y_neg_holes;
    
    fk_vk_neg       = trapz(k_x/a_0, integrand_neg);
    fk_neg          = trapz(k_x/a_0, integrand_fermi_neg);
    fk_vk_neg_holes = trapz(k_x/a_0, integrand_neg_holes);
    
    Id_Vd_electrons(index) = -q_si*w*(sum(fk_vk - fk_vk_neg));
    Id_Vd_holes(index)     =  q_si*w*(sum(fk_vk_holes - fk_vk_neg_holes));
    fk_vk_neg_list(index)  = sum(fk_vk_neg);
    v_inj_fwd = nansum(fk_vk./fk_forward);
    v_inj_back = nansum(fk_vk_neg./fk_neg);
    vinj_Vd(index)         = v_inj_fwd - v_inj_back;
    
%     fk_vk_neg_new(index,:) = trapz(k_x, trapz(1/1+exp((E-(mu + Vd(index)))/kbT)));
    
end

figure();
plot(Vd, (Id_Vd_electrons+Id_Vd_holes));
title(['Total Current']);
xlabel(['Drain Voltage']);
ylabel(['Drain Current (A / \mu m)']);
figure();
plot(Vd, Id_Vd_electrons)
title(['Electrons']);
xlabel(['Drain Voltage']);
ylabel(['Drain Current (A / \mu m)']);
figure();
plot(Vd, Id_Vd_holes)
title(['Holes'])
xlabel(['Drain Voltage']);
ylabel(['Drain Current (A / \mu m)']);

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
xlabel(['Drain Voltage (V)']);
ylabel(['Drain Current (\mu m / \mu m)']);
figure();
plot(Vd, Id_Vd_electrons)
title(['Electrons']);
xlabel(['Drain Voltage (V)']);
ylabel(['Drain Current (\mu m / \mu m)']);
figure();
plot(Vd, Id_Vd_holes)
title(['Holes'])
xlabel(['Drain Voltage (V)']);
ylabel(['Drain Current (\mu m / \mu m)']);

figure();
plot(Vd, vinj_Vd)
title(['Injection Velocity'])
xlabel(['Drain Voltage (V)']);
ylabel(['Injection Velocity (m/s)']);