%Calculating ID_VD curves using effective mass approximation
clear all
vf = 1e6; %m/s Fermi Velocity of graphene
kbT_eV = 0.026; %in eV
a_0    = 1.42e-10; %Graphene lattice constant
h_bar_eV  = 6.582e-16; %eV*s
h_bar_J = 1.055e-34; %J*s
q      = 1.6e-19; %Coulomb
w      = 1e-6; % 1um wide transistor
l      = 100e-9; %100nm long transistor
Vd     = linspace(0,1,100); %Volts
mu_s     = 0; 

a            = 3/2;
b            = sqrt(3)/2;
kmax_x       = pi/(a); 
kmax_y       = 2*pi/(3*b);
kmin_y       = pi / (3*b);

% Ek = graphene_bandstructure();
x_resolution = 50;

k_x = linspace(0, kmax_x, x_resolution);
k_y_limit = linspace(kmax_y, kmin_y, x_resolution);
num_bands = 6;
% Recall v=(1/h_bar)dE/dk
y_resolution = 50;
E = zeros(x_resolution, y_resolution, num_bands);
E1 = zeros(x_resolution, y_resolution);
E2 = zeros(x_resolution, y_resolution);
E3 = zeros(x_resolution, y_resolution);
E4 = zeros(x_resolution, y_resolution);
E5 = zeros(x_resolution, y_resolution);
E6 = zeros(x_resolution, y_resolution);

for x_index = 1:x_resolution
    Ek_y = zeros(y_resolution, num_bands);
    k_y = linspace(-k_y_limit(x_index), k_y_limit(x_index), y_resolution);
    for y_index=1:y_resolution
        temp_E = graphene_E_k(-k_x(x_index), k_y(y_index));
%         temp_V = a_0 * (q/h_bar) .* k_to_v(-k_x(x_index), k_y(y_index), 0.1);
        E(x_index, y_index, :) = temp_E;
        E1(x_index, y_index) = temp_E(1);
        E2(x_index, y_index) = temp_E(2);
        E3(x_index, y_index) = temp_E(3);
        E4(x_index, y_index) = temp_E(4);
        E5(x_index, y_index) = temp_E(5);
        E6(x_index, y_index) = temp_E(6);
    end
end
% 

%Functions defining integrals for holes and electrons
fe = @(y, eta) y./(1 + exp(y - eta)); 
fh = @(y, eta) y.*(1./(1 + exp(y + eta)));

%Density of states in SI units
D = @(E) 2*abs(E)/(pi*(h_bar_eV).^2*vf.^2); 

%Simple fermi function
f = @(E) 1./(1 + exp((E)./kbT_eV)); 

%Initiat
Fplus_e = zeros(1, length(Vd));
Fminus_e = zeros(1, length(Vd));
Fplus_h = zeros(1, length(Vd));
Fminus_h = zeros(1, length(Vd));
Id = zeros(1, length(Vd));
Id_e = zeros(1, length(Vd));
Id_h = zeros(1, length(Vd));
Vg = 0;

%Idealized MOSFET capacitances to begin with, full gate control
Cg = 0.022*l*w; % Farads Hafnia gate, 10nm (thick), 25 dielectric
Cd = 0;
Cs = 0;
alpha_s = Cs ./ (Cg + Cd + Cs);
alpha_g = Cg ./ (Cg + Cd + Cs);
alpha_d = Cd ./ (Cg + Cd + Cs);

%Convergence Error Percentage
epsilon = 1e-5;


%%
%Id-Vd Plot
Vg = 0.0;
damping = 0.1;
vinj = zeros(size(Vd));
vthermal = zeros(size(Vd));
for i = 1:length(Vd)
    %Potential created by gate bias 
    UL = -(alpha_g*Vg + alpha_d*Vd(i));

    %Equillibrium number of electroncs
    N0_e = integral(@(y)(f(y).*D(y)), 0, Inf)*l*w; 
    %Equillibrium number of holes
    N0_h = integral(@(y)(1-f(-y)).*D(y), 0, Inf)*l*w; 

    N0 = N0_e - N0_h;

    Uc = q/(Cg + Cd + Cs);
    dN = 0;
    Uscf = UL + Uc*dN; 
    Ulast = Uscf*2;
    q = 1;
    q_si = 1.6e-19;
    U_0 = -0.15;
    new_U = 1;
    ind = 0;
%   while abs(new_U - U_0) > epsilon
%     while 1 < 0
%         ind = ind + 1;
%         %Electron Concentrations
%         N1_e = 0.5*integral(@(y)(f(y-(mu_s + U_0)).*D(y)), 0, Inf)*l*w; 
%         N2_e = 0.5*integral(@(y)(f(y-(mu_s + U_0 + Vd(i))).*D(y)), 0, Inf)*l*w;
% 
%         %Hole Concentrations
%         N1_h = 0.5*integral(@(y)((f(y+mu_s + U_0))).*D(y), 0, Inf)*l*w; 
%         N2_h = 0.5*integral(@(y)(f(y+mu_s + U_0 - Vd(i)).*D(y)), 0, Inf)*l*w;
% 
%         N1 = N1_e - N1_h;
%         N2 = N2_e - N2_h;
%         delta_N = (N1_h + N2_h - N1_e - N2_e)*l*w;
% 
%         new_U = -q*alpha_g*Vg - q*alpha_d*Vd(i) + q*q_si*delta_N/Cg;
%         U_0 = U_0 + damping * (new_U - U_0);
%     end
%     disp(ind);
    %U_0 = channel_sc_potential(E, x_resolution, y_resolution, Vg, Vd(i), 0);
    U_L = U_0;
%     disp(U_0);
    eta_s = (mu_s + U_L)/ kbT_eV;
    eta_d = (mu_s + Vd(i) + U_L) / kbT_eV; %pretty sure this is the correct conversion?
    %Electrons
    Fplus_e(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf) * integral(@(y)fe(y, eta_s), 0, Inf);
    Fminus_e(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf) * integral(@(y)fe(y, eta_d), 0, Inf);

    %vinj(i) = (Fplus_e(i) + Fplus_h(i))/(integral(@(E)(f(E - eta_s)),0, Inf));

    %Holes
    Fplus_h(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf) * integral(@(y)fh(y, eta_s), 0, Inf);
    Fminus_h(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf) * integral(@(y)fh(y, eta_d), 0, Inf);

    N_temp(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf^2) * integral(@(y)fh(y, eta_s), 0, Inf);
    N_temp_neg(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf^2) * integral(@(y)fh(y, eta_d), 0, Inf);
    N_temp_h(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf^2) * integral(@(y)fh(y, eta_s), 0, Inf);
    N_temp_neg_h(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf^2) * integral(@(y)fe(y, eta_d), 0, Inf);
    Id_e(i) = -2 * q_si * w *(Fplus_e(i)- Fminus_e(i));
    Id_h(i) = 2 * q_si * w *(Fplus_h(i)- Fminus_h(i));
    vinj(i) = (Id_e(i) + Id_h(i))./(2*q_si*w * (N_temp(i) + N_temp_neg_h(i) + N_temp(i) + N_temp_neg_h(i)));
%     vthermal(i) = Fplus_h(i)/N_temp(i);
end
Id = Id_e + Id_h;

%%
subplot(3, 1, 1)
plot(Vd, Id);
title(['Total Current']);
subplot(3, 1, 2)
plot(Vd, Id_e)
title(['Electrons']);
ylabel(['Drain Current (A / \mu m)']);
subplot(3, 1, 3)
plot(Vd, Id_h)
title(['Holes'])
xlabel(['Drain Voltage']);
figure();
plot(Vd, vinj);
title(['Injection Velocity'])
ylabel(['Injection Velocity (m/s)']);
xlabel(['Drain Voltage']);


%%

%Making Id-Vg Plot

Vd2 = 0.2;
Vg = linspace(-2, 2, 50);
Ig_e = zeros(size(Vg));
Ig_h = zeros(size(Vg));
Fgplus_e = zeros(size(Vg));
Fgminus_e = zeros(size(Vg));
U = [];
for i = 1:length(Vg)
    
    %Solving Self Consistent Potential Problem
   %Potential created by gate bias 
   UL = -(alpha_g*Vg(i) + alpha_d*Vd2);

   %Equillibrium number of electroncs
   N0_e = integral(@(y)(f(y).*D(y)), 0, Inf)*l*w; 
   %Equillibrium number of holes
   N0_h = integral(@(y)(f(y)).*D(y), 0, Inf)*l*w; 
   
   N0 = N0_e - N0_h;
   
    Uc = q/(Cg + Cd + Cs);
    dN = 0;
    Uscf = UL + Uc*dN; 
    Ulast = Uscf*2;
    q = 1;
    q_si = 1.6e-19;
    U_0 = -0.15;
    new_U = 1;
    ind = 0;
    while abs(new_U - U_0) > epsilon
       ind = ind + 1;
       %Electron Concentrations
       N1_e = 0.5*integral(@(y)(f(y-(mu_s + U_0)).*D(y)), 0, Inf)*l*w; 
       N2_e = 0.5*integral(@(y)(f(y-(mu_s + U_0 + Vd2)).*D(y)), 0, Inf)*l*w;
       
       %Hole Concentrations
       N1_h = 0.5*integral(@(y)((f(y+mu_s + U_0))).*D(y), 0, Inf)*l*w; 
       N2_h = 0.5*integral(@(y)(f(y+mu_s + U_0 - Vd2).*D(y)), 0, Inf)*l*w;
       
       N1 = N1_e - N1_h;
       N2 = N2_e - N2_h;
       delta_N = (N1_h + N2_h - N1_e - N2_e)*l*w;

       new_U = -q*alpha_g*Vg(i) - q*alpha_d*Vd2 + q*q_si*delta_N/Cg;
       U_0 = U_0 + damping * (new_U - U_0);
    end
   disp(ind);
   %U_0 = channel_sc_potential(E, x_resolution, y_resolution, Vg, Vd(i), 0);
   U_L = U_0;
   
   eta_s = (mu_s + U_0) / kbT_eV;
   eta_d = (mu_s + U_0 + Vd2) / kbT_eV; %pretty sure this is the correct conversion?
   %Electrons
   Fgplus_e(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf) * integral(@(y)fe(y, eta_s), 0, Inf);
   Fgminus_e(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf) * integral(@(y)fe(y, eta_d), 0, Inf);
   
   %Holes
   Fplus_h(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf) * integral(@(y)fh(y, eta_s), 0, Inf);
   Fminus_h(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf) * integral(@(y)fh(y, eta_d), 0, Inf);
   
   
   Ig_e(i) = -2 * q_si * w *(Fgplus_e(i)- Fgminus_e(i));
   Ig_h(i) = +2 * q_si * w *(Fplus_h(i)- Fminus_h(i));
end
Ig = Ig_e + Ig_h;


figure();


subplot(3, 1, 1)
plot(Vg, Ig);
title(['Total Current']);
subplot(3, 1, 2)
plot(Vg, Ig_e)
title(['Electrons']);
ylabel(['Drain Current (A / \mu m)']);
subplot(3, 1, 3)
plot(Vg, Ig_h)
xlabel(['Gate Voltage']);
title(['Holes'])