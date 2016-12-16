function[mu_s] = channel_sc_potential(E, x_resolution, y_resolution, vgs, vds, guess_U_0);
%Self Consistent Model generation
% a_0          = 1.42e-10; %Graphene lattice constant
% a            = 3/2;
% b            = sqrt(3)/2;
% kmax_x       = pi/(a); 
% kmax_y       = 2*pi/(3*b);
% kmin_y       = pi / (3*b);
% k_y_limit = linspace(kmax_y, kmin_y, x_resolution);
% k_x = linspace(0, kmax_x, x_resolution);
% num_bands = 6;
% kbT = 0.026;

w      = 1e-6; % 1um wide transistor
l      = 100e-9; %100nm long transistor

damping = 0.5;

%Idealized MOSFET capacitances to begin with, full gate control
Cg = 0.022*l*w; % Farads Hafnia gate, 10nm (thick), 25 dielectric
Cd = 0;
Cs = 0;
Ct = Cg + Cd + Cs;
alpha_s = Cs ./ (Ct);
alpha_g = Cg ./ (Ct);
alpha_d = Cd ./ (Ct);

%Convergence Voltage Difference
epsilon = 1e-4;

%Initial Guess
%Note: Everything is in eV.
%Everything follows from Datta, page 17 (32 in PDF).
q = 1;
q_si = 1.6e-19;
U_0 = guess_U_0;
new_U = 1;


while new_U - U_0 > epsilon
    [e, h] = find_concentrations(E, x_resolution, y_resolution, U_0);
    [e2, h2] = find_concentrations(E, x_resolution, y_resolution, U_0-vds);
    delta_N = (h + h2 - e - e2)*l*w;

    new_U = -q*alpha_g*vgs - q*alpha_d*vds + q*q_si*delta_N/Cg;
%     disp(new_U - U_0)
    U_0 = U_0 + damping * (new_U - U_0);
end
mu_s = U_0;
%%
% 
% D = @(E) 2*abs(E)/(pi*(h_bar_eV).^2*vf.^2); %Density of states in SI units
% f = @(E) 1./(1 + exp((E)./kbT_eV)); %Simple fermi function
% Vd2 = 0.4;
% Vg = linspace(-2, 2, 10);
% Ig_e = zeros(size(Vg));
% Ig_h = zeros(size(Vg));
% Fgplus_e = zeros(size(Vg));
% Fgminus_e = zeros(size(Vg));
% U = [];
% for i = 1:length(Vg)
%     
%     %Solving Self Consistent Potential Problem
%    
%    %Potential created by gate bias 
%    UL = -(alpha_g*Vg(i) + alpha_d*Vd2);
% 
%    %Equillibrium number of electroncs
%    N0_e = integral(@(y)(f(y).*D(y)), 0, Inf)*l*w; 
%    %Equillibrium number of holes
%    N0_h = integral(@(y)((1-f(y)).*D(y)), 0, Inf)*l*w; 
%    
%    N0 = N0_e - N0_h;
%    
%    Uc = q/(Cg*l*w + Cd*l*w + Cs*l*w);
%    dN = 0;
%    Uscf = UL + Uc*dN; 
%    Ulast = Uscf*2;
%    
%    ind = 0;
%    while Ulast-Uscf > epsilon
%        ind = ind + 1;
%        %Electron Concentrations
%        N1_e = 0.5*integral(@(y)(f(y+mu_s + Uscf).*D(y)), 0, Inf)*l*w; 
%        N2_e = 0.5*integral(@(y)(f(y+mu_s + Uscf + Vd2).*D(y)), 0, Inf)*l*w;
%        
%        %Hole Concentrations
%        N1_h = 0.5*integral(@(y)((1-f(-y+mu_s + Uscf)).*D(y)), 0, Inf)*l*w; 
%        N2_h = 0.5*integral(@(y)((1-f(-y+mu_s + Uscf + Vd2)).*D(y)), 0, Inf)*l*w;
%        
%        N1 = N1_e - N1_h;
%        N2 = N2_e - N2_h;
%        
%        dN = N1 + N2 - N0;
%        disp(dN);
%        Ulast = Uscf;
%        Uscf = Ulast + ((UL + Uc*dN) - Ulast)/2; %Setting to half of the difference?
%        disp(['N1_e:' num2str(N1_e) ' N1_h:' num2str(N1_h) ' N2_e:' num2str(N2_e) ' N2_h:' num2str(N2_h)]);
%        disp(['Iteration: ' num2str(ind) ' Uscf:' num2str(Uscf) ' UL:' num2str(UL)]); 
%        disp(num2str(Ulast - Uscf));
%    end
%    U = [U Uscf];
%    
%    eta_s = (mu_s + Uscf) / kbT_eV;
%    eta_d = (mu_s + Uscf - Vd2) / kbT_eV; %pretty sure this is the correct conversion?
%    %Electrons
%    Fgplus_e(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf) * integral(@(y)fe(y, eta_s), 0, Inf);
%    Fgminus_e(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf) * integral(@(y)fe(y, eta_d), 0, Inf);
%    
%    %Holes
%    Fplus_h(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf) * integral(@(y)fh(y, eta_s), 0, Inf);
%    Fminus_h(i) = 1 / (2*pi^2) * (kbT_eV / (h_bar_eV)).^2 * 1/(vf) * integral(@(y)fh(y, eta_d), 0, Inf);
%    
%    
%    Ig_e(i) = 2 * q * w *(Fgplus_e(i)- Fgminus_e(i));
%    Ig_h(i) = -2 * q * w *(Fplus_h(i)- Fminus_h(i));
% end
% Ig = Ig_e + Ig_h;
