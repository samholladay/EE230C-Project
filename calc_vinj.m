function [ output_args ] = calc_vinj(k, Ek)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
kbT = 0.026;
h_bar=(6.626e-34)/(2*pi);
q=1.6e-19;

vth_m=sqrt((2*kbT)/(pi*m_eff));
v_inj_m=vth_m.*(1-exp(-Vd./kbT))./(1-exp(-Vd./kbT));

end

