h_bar_SI  = (6.626e-34)/(2*pi); % SI
h_bar        = 6.582e-16; %eV
% q            = 1.6e-19;
q            = 1;
q_si         = 1.6e-19; 
kbT          = 0.026; %eV
a_0          = 1.42e-10; %Graphene lattice constant

w            = 1e-6; % How wide is the transistor?
Vd           = linspace(0,1.5,20); %Volts
mu           = 0.0;
y_resolution = 100;
num_bands    = 6;

a            = 3/2;
b            = sqrt(3)/2;
kmax_x       = pi/(a); 
kmax_y       = 2*pi/(3*b);
kmin_y       = pi / (3*b);

%Ek = graphene_bandstructure(); CHANGE THIS
length_of_Ek = size(Ek,1);

k_x = linspace(0, kmax_x, length_of_Ek);
k_x = k_x(3:length_of_Ek);
k_y_limit = linspace(kmax_y, kmin_y, length_of_Ek);
k_y_limit = k_y_limit(3:length_of_Ek);

k_step = pi/(a_0*size(Ek,1));
% Recall v=(1/h_bar)dE/dk
% Chopping off the first element of Ek because it is 0.
Ek = Ek(2:length_of_Ek, 2:length_of_Ek, :); CHANGE THIS LINE BECAUSE Ek HAS DIMENSION 3
grad_x = diff(Ek,kx_step,1);
grad_y = diff(Ek, ky_step,2);
v = (1/h_bar).*(grad_x+grad_y);

Ek_cut = Ek(2:length_of_Ek-1,2:length_of_Ek-1,:);
num_of_x_indicies = size(Ek_cut,1);
for band=1:num_bands
    fermi = 1./(1+exp((Ek(:,:,band)-mu)./kbT)); % What is mu?
    fermi_holes = 1./(1+exp((Ek(:,:,band)+mu)./kbT));
    integrand_elec  = (1/(4.*pi^2)).*fermi.*v;
    integrand_holes = (1/(4.*pi^2)).*fermi_holes.*v;
    integrand_fk_elec = (1/(4.*pi^2)).*fermi;
    integrand_fk_holes = (1/(4.*pi^2)).*fermi_holes;
    fk_vk_elec = kx_step.*ky_step.*trapz(kx,trapz(ky, integrand_elec, 2));
    fk_vk_holes = kx_step.*ky_step.*trapz(kx,trapz(ky, integrand_holes, 2));
    fk_elec = kx_step.*ky_step.*trapz(kx,trapz(ky, integrand_fk_elec, 2));
    fk_holes = kx_step.*ky_step.*trapz(kx,trapz(ky, integrand_fk_holes, 2));
    vth = (fk_vk_elec./fk_elec) + (fk_vk_holes./fk_holes);
end

for Vd=0:Vd_step:Vd_max
    for band2=1:num_bands
        fermi_back = 1./(1+exp((Ek(:,:,band2)-(mu-Vd))./kbT)); % What is mu?
        fermi_holes_back = 1./(1+exp((Ek(:,:,band2)+(mu-Vd))./kbT));
        integrand_elec  = (1/(4.*pi^2)).*fermi_back.*v;
        integrand_holes = (1/(4.*pi^2)).*fermi_holes.*v;
        fk_vk_elec = kx_step.*ky_step.*trapz(kx,trapz(ky, integrand_elec, 2));
        fk_vk_holes = kx_step.*ky_step.*trapz(kx,trapz(ky, integrand_holes, 2));
    end
end

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