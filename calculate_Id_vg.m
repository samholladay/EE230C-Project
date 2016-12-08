function [Id, Vd] = calculate_Id_vg(k, Ek)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
clear all;
kbT = 0.026;
h_bar=(6.626e-34)/(2*pi);
q=1.6e-19;
w=10e-9; % How wide is the transistor?
ss=0.01; % Step size
Vd=0:ss:8;

% Recall v=(1/h_bar)dE/dk
grad = diff(Ek)/ss;
v = (1/h_bar).*grad;
fermi=1./(1+exp((Ek-mu)./kbT)); % What is mu?
integrand = (1/(2*pi))*fermi*v;
fk_vk =integral(integrand,0,1e15); % Integrate to infinity
Id=-q.*w.*fk_vk;
%v_inj=vth.*(1-exp(-Vd./kbT))./(1-exp(-Vd./kbT));

plot(Vd,Id);
xlabel('Vd');
end

