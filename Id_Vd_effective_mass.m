%Calculating ID_VD curves using effective mass approximation

vf = 1e6; %Fermi Velocity of graphene
kbT_eV = 0.026; %in eV
a_0    = 1.42e-10; %Graphene lattice constant
h_bar  = 6.582e-16; %eV*s
q      = 1.6e-19; %Coulomb
w      = 1e-6; % How wide is the transistor?
Vd     = linspace(0,0.5,500); %Volts
mu_s     = 0; 
f = @(y, eta) y./(1 + exp(y - eta)); 
Fplus = zeros(1, length(Vd));
Fminus = zeros(1, length(Vd));
I = zeros(1, length(Vd));

for i = 1:length(Vd)
   eta_s = mu_s / kbT_eV;
   eta_d = mu_s-Vd(i) / kbT_eV; %pretty sure this is the correct conversion?
   Fplus(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar)).^2 * 1/(vf) * integral(@(y)f(y, eta_s), 0, Inf);
   Fminus(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar)).^2 * 1/(vf) * integral(@(y)f(y, eta_d), 0, Inf);
   I(i) = -2 * q * w *(Fplus(i) - Fminus(i));
end

plot(Vd, I);