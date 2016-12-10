function [Id, Vd] = calculate_Id_vg(Ek)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
kbT = 0.026;
h_bar=(6.626e-34)/(2*pi);
q=1.6e-19;
w=10e-9; % How wide is the transistor?
ss=0.01; % Step size
Vd=0:ss:8;
mu=0;
k_step = 1;

% Recall v=(1/h_bar)dE/dk
grad = diff(Ek)/k_step;
v = (1/h_bar).*grad;
length=size(Ek,1);
Ek_cut=Ek(2:length,:);
fermi=1./(1+exp((Ek_cut-mu)./kbT)); % What is mu?
integrand = (1/(2.*pi)).*fermi.*v;
x=linspace(0,1e6,size(Ek_cut,1));
% MUST USE X TO DEFINE LIMITS OF INTEGRATION
fk_vk =trapz(x,integrand); % Integrate to infinity
fk = trapz(x,fermi);
vth = nansum(fk_vk./fk);
%Id=-q.*w.*fk_vk;
Vd=linspace(0,2,500);
v_inj=vth.*(1-exp(-Vd./kbT))./(1+exp(-Vd./kbT));
Id=-q.*w.*v_inj;

plot(Vd,Id);
%xlabel('Vd');
clear all;
end

