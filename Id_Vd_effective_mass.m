%Calculating ID_VD curves using effective mass approximation

vf = 1e6; %m/s Fermi Velocity of graphene
kbT_eV = 0.026; %in eV
a_0    = 1.42e-10; %Graphene lattice constant
h_bar_eV  = 6.582e-16; %eV*s
h_bar_J = 1.055e-34; %J*s
q      = 1.6e-19; %Coulomb
w      = 1e-6; % 1um wide transistor
l      = 100e-9; %100nm long transistor
Vd     = linspace(0,1,500); %Volts
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
Vg = 0;
for i = 1:length(Vd)
   %Potential created by gate bias 
   UL = -(alpha_g*Vg + alpha_d*Vd(i));

   %Equillibrium number of electroncs
   N0_e = integral(@(y)(f(y).*D(y)), 0, Inf)*l*w; 
   %Equillibrium number of holes
   N0_h = integral(@(y)((1-f(y)).*D(y)), 0, Inf)*l*w; 
   
   N0 = N0_e - N0_h;
   
   Uc = q/(Cg*l*w + Cd*l*w + Cs*l*w);
   dN = 0;
   Uscf = UL + Uc*dN; 
   Ulast = Uscf*2;
   
   ind = 0;
   while Ulast-Uscf > epsilon
       ind = ind + 1;
       %Electron Concentrations
       N1_e = 0.5*integral(@(y)(f(y+mu_s + Uscf).*D(y)), 0, Inf)*l*w; 
       N2_e = 0.5*integral(@(y)(f(y+mu_s + Uscf + Vd(i)).*D(y)), 0, Inf)*l*w;
       
       %Hole Concentrations
       N1_h = 0.5*integral(@(y)((1-f(-y+mu_s + Uscf)).*D(y)), 0, Inf)*l*w; 
       N2_h = 0.5*integral(@(y)((1-f(-y+mu_s + Uscf + Vd(i))).*D(y)), 0, Inf)*l*w;
       
       N1 = N1_e - N1_h;
       N2 = N2_e - N2_h;
       
       dN = N1 + N2 - N0;
       disp(dN);
       Ulast = Uscf;
       Uscf = Ulast + ((UL + Uc*dN) - Ulast)/2; %Setting to half of the difference?
       disp(['N1_e:' num2str(N1_e) ' N1_h:' num2str(N1_h) ' N2_e:' num2str(N2_e) ' N2_h:' num2str(N2_h)]);
       disp(['Iteration: ' num2str(ind) ' Uscf:' num2str(Uscf) ' UL:' num2str(UL)]); 
       disp(num2str(Ulast - Uscf));
   end
    
   eta_s = (mu_s + UL)/ kbT_eV;
   eta_d = (mu_s + Vd(i) + UL) / kbT_eV; %pretty sure this is the correct conversion?
   %Electrons
   Fplus_e(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf) * integral(@(y)fe(y, eta_s), 0, Inf);
   Fminus_e(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf) * integral(@(y)fe(y, eta_d), 0, Inf);
   
   %Holes
   Fplus_h(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf) * integral(@(y)fh(y, eta_s), 0, Inf);
   Fminus_h(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf) * integral(@(y)fh(y, eta_d), 0, Inf);
   
   
   Id_e(i) = 2 * q * w *(Fplus_e(i)- Fminus_e(i));
   Id_h(i) = 2 * q * w *(Fplus_h(i)- Fminus_h(i));
end
Id = Id_e + Id_h;


%Self Consistent Model generation

%Idealized MOSFET capacitances to begin with, full gate control
Cg = 0.022*l*w; % Farads Hafnia gate, 10nm (thick), 25 dielectric
Cd = 0;
Cs = 0;
alpha_s = Cs ./ (Cg + Cd + Cs);
alpha_g = Cg ./ (Cg + Cd + Cs);
alpha_d = Cd ./ (Cg + Cd + Cs);

%Convergence Error Percentage
epsilon = 1e-5;

D = @(E) 2*abs(E)/(pi*(h_bar_eV).^2*vf.^2); %Density of states in SI units
f = @(E) 1./(1 + exp((E)./kbT_eV)); %Simple fermi function
Vd2 = 0.4;
Vg = linspace(-2, 2, 10);
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
   N0_h = integral(@(y)((1-f(y)).*D(y)), 0, Inf)*l*w; 
   
   N0 = N0_e - N0_h;
   
   Uc = q/(Cg*l*w + Cd*l*w + Cs*l*w);
   dN = 0;
   Uscf = UL + Uc*dN; 
   Ulast = Uscf*2;
   
   ind = 0;
   while Ulast-Uscf > epsilon
       ind = ind + 1;
       %Electron Concentrations
       N1_e = 0.5*integral(@(y)(f(y-(mu_s + Uscf)).*D(y)), 0, Inf)*l*w; 
       N2_e = 0.5*integral(@(y)(f(y-(mu_s + Uscf + Vd2)).*D(y)), 0, Inf)*l*w;
       
       %Hole Concentrations
       N1_h = 0.5*integral(@(y)((1-f(-y-(mu_s + Uscf))).*D(y)), 0, Inf)*l*w; 
       N2_h = 0.5*integral(@(y)((1-f(-y-(mu_s + Uscf + Vd2))).*D(y)), 0, Inf)*l*w;
       
       N1 = N1_e - N1_h;
       N2 = N2_e - N2_h;
       
       dN = N1 + N2 - N0;
       disp(dN);
       Ulast = Uscf;
       Uscf = Ulast + ((UL + Uc*dN) - Ulast)/2; %Setting to half of the difference?
       disp(['N1_e:' num2str(N1_e) ' N1_h:' num2str(N1_h) ' N2_e:' num2str(N2_e) ' N2_h:' num2str(N2_h)]);
       disp(['Iteration: ' num2str(ind) ' Uscf:' num2str(Uscf) ' UL:' num2str(UL)]); 
       disp(num2str(Ulast - Uscf));
   end
   U = [U Uscf];
   
   eta_s = (mu_s + Uscf) / kbT_eV;
   eta_d = (mu_s + Uscf + Vd2) / kbT_eV; %pretty sure this is the correct conversion?
   %Electrons
   Fgplus_e(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf) * integral(@(y)fe(y, eta_s), 0, Inf);
   Fgminus_e(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf) * integral(@(y)fe(y, eta_d), 0, Inf);
   
   %Holes
   Fplus_h(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf) * integral(@(y)fh(y, eta_s), 0, Inf);
   Fminus_h(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf) * integral(@(y)fh(y, eta_d), 0, Inf);
   
   
   Ig_e(i) = 2 * q * w *(Fgplus_e(i)- Fgminus_e(i));
   Ig_h(i) = -2 * q * w *(Fplus_h(i)- Fminus_h(i));
end
Ig = Ig_e + Ig_h;

subplot(3, 1, 1)
plot(Vd, Id*1e3);
title(['Total Current']);
xlabel(['Drain Voltage']);
ylabel(['Drain Current (\mu m / \mu A)']);
subplot(3, 1, 2)
plot(Vd, Id_e*1e3)
xlabel(['Drain Voltage']);
title(['Electrons']);
ylabel(['Drain Current (\mu m / \mu A)']);
subplot(3, 1, 3)
plot(Vd, Id_h*1e3)
title(['Holes'])
xlabel(['Drain Voltage']);
ylabel(['Drain Current (\mu m /mA)']);

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