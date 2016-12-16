%Calculating ID_VD curves using effective mass approximation

vf = 1e6; %m/s Fermi Velocity of graphene
kbT_eV = 0.026; %in eV
a_0    = 1.42e-10; %Graphene lattice constant
h_bar_eV  = 6.582e-16; %eV*s
h_bar_J = 1.055e-34; %J*s
q      = 1.6e-19; %Coulomb
w      = 1e-6; % How wide is the transistor?
Vd     = linspace(0,-5.5,500); %Volts
mu_s     = 0; 
%Functions defining integrals for holes and electrons
fe = @(y, eta) y./(1 + exp(y - eta)); 
fh = @(y, eta) y.*(1./(1 + exp(y + eta)));
Fplus_e = zeros(1, length(Vd));
Fminus_e = zeros(1, length(Vd));
Fplus_h = zeros(1, length(Vd));
Fminus_h = zeros(1, length(Vd));
Id = zeros(1, length(Vd));
Id_e = zeros(1, length(Vd));
Id_h = zeros(1, length(Vd));
U = 0.3;
for i = 1:length(Vd)
   eta_s = (mu_s + U)/ kbT_eV;
   eta_d = (mu_s-Vd(i) + U) / kbT_eV; %pretty sure this is the correct conversion?
   %Electrons
   Fplus_e(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf) * integral(@(y)fe(y, eta_s), 0, Inf);
   Fminus_e(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf) * integral(@(y)fe(y, eta_d), 0, Inf);
   
   %Holes
%   Fplus_h(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar)).^2 * 1/(vf) * integral(@(y)fh(y, eta_s), 0, Inf);
%   Fminus_h(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar)).^2 * 1/(vf) * integral(@(y)fh(y, -eta_d), 0, Inf);
   
   
   Id_e(i) = 2 * q * w *(Fplus_e(i)- Fminus_e(i));
   Id_h(i) = 2 * q * w *(Fplus_h(i)- Fminus_h(i));
end
Id = Id_e + Id_h;


%Self Consistent Model generation

%Idealized MOSFET capacitances to begin with, full gate control
alpha_s = 0;
alpha_g = 1;
alpha_d = 0;


D = @(E) 2*abs(E) * q./(pi*(h_bar_J).^2*vf.^2); %Density of states in SI units
f = @(E, mu) 1./(1 + exp((E - mu)./kbT));
Vd2 = 0.1;
Vg = linspace(0, 1, 100);
Ig_e = zeros(size(Vg));
Ig_h = zeros(size(Vg));
Fgplus_e = zeros(size(Vg));
Fgminus_e = zeros(size(Vg));
for i = 1:length(Vg)
   U = (alpha_g*Vg(i) + alpha_d*Vd2);
   disp(U);
   eta_s = (mu_s + U) / kbT_eV;
   eta_d = (mu_s + U - Vd2) / kbT_eV; %pretty sure this is the correct conversion?
   %Electrons
   Fgplus_e(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf) * integral(@(y)fe(y, eta_s), 0, Inf);
   Fgminus_e(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf) * integral(@(y)fe(y, eta_d), 0, Inf);
   
   %Holes
%   Fplus_h(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar)).^2 * 1/(vf) * integral(@(y)fh(y, eta_s), 0, Inf);
%   Fminus_h(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar)).^2 * 1/(vf) * integral(@(y)fh(y, -eta_d), 0, Inf);
   
   
   Ig_e(i) = 2 * q * w *(Fgplus_e(i)- Fgminus_e(i));
   Ig_h(i) = 2 * q * w *(Fplus_h(i)- Fminus_h(i));
end
Ig = Ig_e + Ig_h;

subplot(3, 1, 1)
plot(Vd, Id);
title(['Total Current']);
xlabel(['Drain Voltage']);
ylabel(['Drain Current (\mu m / \mu m)']);
subplot(3, 1, 2)
plot(Vd, Id_e)
xlabel(['Drain Voltage']);
title(['Electrons']);
ylabel(['Drain Current (\mu m / \mu m)']);
subplot(3, 1, 3)
plot(Vd, Id_h)
title(['Holes'])
xlabel(['Drain Voltage']);
ylabel(['Drain Current (\mu m / \mu m)']);

figure();

subplot(3, 1, 1)
plot(Vg, Ig);
title(['Total Current']);
xlabel(['Gate Voltage']);
ylabel(['Drain Current (\mu m / \mu m)']);
subplot(3, 1, 2)
plot(Vg, Ig_e)
xlabel(['Gate Voltage']);
title(['Electrons']);
ylabel(['Drain Current (\mu m / \mu m)']);
subplot(3, 1, 3)
plot(Vg, Ig_h)
title(['Holes'])
xlabel(['Gate Voltage']);
ylabel(['Drain Current (\mu m / \mu m)']);