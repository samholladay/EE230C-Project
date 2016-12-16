function[e_concentration, h_concentration] = find_concentrations(E, x_res, y_res, mu);
a_0          = 1.42e-10; %Graphene lattice constant
a            = 3/2;
b            = sqrt(3)/2;
kmax_x       = pi/(a); 
kmax_y       = 2*pi/(3*b);
kmin_y       = pi / (3*b);
k_y_limit = linspace(kmax_y, kmin_y, x_res);
k_x = linspace(0, kmax_x, x_res);
num_bands = 6;
kbT = 0.026;

E1 = reshape(E(:,:,1),x_res, y_res);
E2 = reshape(E(:,:,2),x_res, y_res);
E3 = reshape(E(:,:,3),x_res, y_res);
E4 = reshape(E(:,:,4),x_res, y_res);
E5 = reshape(E(:,:,5),x_res, y_res);
E6 = reshape(E(:,:,6),x_res, y_res);

fermi_across_y       = zeros(x_res, num_bands);
fermi_across_y_holes = zeros(x_res, num_bands);
for x_index=1:x_res
% for x_index=1:2
    Ek_y = [E1(x_index,:)' E2(x_index,:)' E3(x_index,:)' E4(x_index,:)' E5(x_index,:)' E6(x_index,:)'];
    k_y = linspace(-k_y_limit(x_index), k_y_limit(x_index), y_res);
    
    fermi       = 1./(1+exp((Ek_y-mu)./kbT)); % What is mu?
    fermi_holes = 1./(1+exp((Ek_y+mu)./kbT));
%     fermi_holes = 1 - fermi_holes;
    fermi_across_y(x_index, :)       = trapz(k_y/a_0, fermi);
    fermi_across_y_holes(x_index, :) = trapz(k_y/a_0, fermi_holes);
end

integrand = (1/(4.*pi^2)).*fermi_across_y;
integrand_holes = (1/(4.*pi^2)).*fermi_across_y_holes;
e_concentration    = sum(trapz(k_x/a_0,integrand));
h_concentration    = sum(trapz(k_x/a_0,integrand_holes));
