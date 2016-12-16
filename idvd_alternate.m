h_bar_SI  = (6.626e-34)/(2*pi); % SI
h_bar        = 6.582e-16; %eV
% q            = 1.6e-19;
q            = 1;
q_si         = 1.6e-19; 
kbT          = 0.026; %eV
a_0          = 1.42e-10; %Graphene lattice constant

w            = 1e-6; % How wide is the transistor?
Vd           = linspace(0,2,20); %Volts
mu           = 0.0;
num_bands    = 6;

a            = 3/2;
b            = sqrt(3)/2;
kmax_x       = pi/(a); 
kmax_y       = 2*pi/(3*b);
kmin_y       = pi / (3*b);
num_of_x_indicies = 200;
y_resolution = 200;


Ek = E_k_relationship(num_of_x_indicies,y_resolution);
disp(size(Ek));
length_of_Ek = size(Ek,1);

k_x = linspace(0, kmax_x, length_of_Ek);
%k_x = k_x(3:length_of_Ek);
k_y_limit = linspace(kmax_y, kmin_y, length_of_Ek);
k_y = linspace(k_y_limit(x_index), -k_y_limit(x_index), y_resolution);
%k_y_limit = k_y_limit(3:length_of_Ek);

k_step = pi/(a_0*size(Ek,1));
ky_step = (2*k_y_limit/a_0);
% Recall v=(1/h_bar)dE/dk
% Chopping off the first element of Ek because it is 0.
%Ek = Ek(2:length_of_Ek, 2:length_of_Ek, :); % CHANGE THIS LINE BECAUSE Ek HAS DIMENSION 3
grad_x = diff(Ek,1,1)/k_step;
grad_y = diff(Ek,1,2)/k_step;
disp(grad_x)
%disp(grad_y)
grad_x = grad_x(:,2:y_resolution,:);
grad_y = grad_y(2:y_resolution,:,:);
v = (1/h_bar).*(grad_x+grad_y);
%disp(v)

Ek = Ek(2:length_of_Ek,2:length_of_Ek,:);
num_of_x_indicies = size(Ek,1);
total_fk_vk_elec_forward = zeros(num_of_x_indicies,6);
total_fk_vk_holes_forward = zeros(num_of_x_indicies,6);
total_vth = zeros(num_of_x_indicies,6);
k_x = k_x(2:y_resolution);
for band=1:num_bands
    fermi = 1./(1+exp((Ek(:,:,band)-mu)./kbT)); % What is mu?
    fermi_holes = 1./(1+exp((Ek(:,:,band)+mu)./kbT));
    integrand_elec  = (1/(4.*pi^2)).*fermi.*v(:,:,band);
    integrand_holes = (1/(4.*pi^2)).*fermi_holes.*v(:,:,band);
    integrand_fk_elec = (1/(4.*pi^2)).*fermi;
    integrand_fk_holes = (1/(4.*pi^2)).*fermi_holes;
    disp(size(k_x))
    disp(size(integrand_elec))
    fk_vk_elec = k_step.*k_step.*trapz(k_x,trapz(integrand_elec, 2));
    fk_vk_holes = k_step.*k_step.*trapz(k_x,trapz(integrand_holes, 2));
    fk_elec = k_step.*k_step.*trapz(k_x,trapz(integrand_fk_elec, 2));
    fk_holes = k_step.*k_step.*trapz(k_x,trapz(integrand_fk_holes, 2));
    vth = (fk_vk_elec./fk_elec) + (fk_vk_holes./fk_holes);
    total_fk_vk_elec_forward(:,band)=fk_vk_elec;
    total_fk_vk_holes_forward(:,band)=fk_vk_holes;
    total_vth(:,band)=vth;
end
disp(total_fk_vk_elec_forward)
elec_forward_current = sum(total_fk_vk_elec_forward,2);
hole_forward_current = sum(total_fk_vk_holes_forward,2);
vth = sum(total_vth,2);
disp(vth)

Id_Vd_electrons_back = zeros(Vd_max/Vd_step,1);
Id_Vd_holes_back= zeros(Vd_max/Vd_step,1);
for index = 1:length(Vd)
    total_fk_vk_elec_back = zeros(num_of_x_indicies,6);
    total_fk_vk_holes_back = zeros(num_of_x_indicies,6);
    for band2=1:num_bands
        fermi_back = 1./(1+exp((Ek(:,:,band2)-(mu-Vd(index)))./kbT)); % What is mu?
        fermi_holes_back = 1./(1+exp((Ek(:,:,band2)+(mu-Vd(index)))./kbT));
        integrand_elec  = (1/(4.*pi^2)).*fermi_back.*v(:,:,band);
        integrand_holes = (1/(4.*pi^2)).*fermi_holes.*v(:,:,band);
        fk_vk_elec = k_step.*k_step.*trapz(k_x,trapz(integrand_elec, 2));
        fk_vk_holes = k_step.*k_step.*trapz(k_x,trapz(integrand_holes, 2));
        total_fk_vk_elec_back(:,band)=fk_vk_elec;
        total_fk_vk_holes_back(:,band)=fk_vk_holes;
    end
    elec_back_current = sum(total_fk_vk_elec_back,2);
    hole_back_current = sum(total_fk_vk_holes_back,2);
    Id_Vd_electrons_back(index)=sum(elec_back_current);
    Id_Vd_holes_back(index)=sum(hole_back_current);
end
Id_Vd_electrons = -q_si.*w.*(Id_Vd_electrons_back+sum(elec_forward_current));
Id_Vd_holes = -q_si.*w.*(Id_Vd_holes_back+sum(hole_forward_current));

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