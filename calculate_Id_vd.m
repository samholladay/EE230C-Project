function [Id, Vd] = calculate_Id_vg(Ek)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
kbT = 0.026;
a = 1.42e-10; %Graphene lattice constant
h_bar=(6.626e-34)/(2*pi);
q=1.6e-19;
w=1e-6; % How wide is the transistor?
Vd=linspace(0,0.5,500);
mu=0;
k_step = pi/(a*size(Ek,1));

% Recall v=(1/h_bar)dE/dk
grad = diff(Ek)/k_step;
v = (1/h_bar).*grad;
length=size(Ek,1);
Ek_cut=Ek(2:length,:);
fermi=1./(1+exp((Ek_cut-mu)./kbT)); % What is mu?
integrand = (1/(2.*pi)).*fermi.*v;
x=linspace(0,1e-10,size(Ek_cut,1));
% MUST USE X TO DEFINE LIMITS OF INTEGRATION
fk_vk =trapz(x,integrand); % Integrate to infinity
fk = trapz(x,fermi);
vth = nansum(fk_vk./fk);
disp(vth)
%Id=-q.*w.*fk_vk;

v_inj=vth.*(1-exp(-Vd./kbT))./(1+exp(-Vd./kbT));
Id=-q.*w.*v_inj;

plot(Vd,Id);
%xlabel('Vd');
end

