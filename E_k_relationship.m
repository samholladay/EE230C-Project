function [E_k] = E_k_relationship(num_of_x_indicies, y_resolution);
a            = 3/2;
b            = sqrt(3)/2;
kmax_x       = pi/(a); 
kmax_y       = 2*pi/(3*b);
kmin_y       = pi / (3*b);

length_of_Ek = num_of_x_indicies;

k_x = linspace(0, kmax_x, length_of_Ek);
k_y_limit = linspace(kmax_y, kmin_y, length_of_Ek);
num_bands = 6;
% Recall v=(1/h_bar)dE/dk
% Chopping off the first element of Ek because it is 0.

E_k = zeros(num_of_x_indicies, y_resolution, num_bands);
for x_index = 1:num_of_x_indicies
    k_y = linspace(k_y_limit(x_index), -k_y_limit(x_index), y_resolution);
    for y_index=1:y_resolution
        temp_E = graphene_E_k(-k_x(x_index), k_y(y_index));
        E_k(x_index, y_index, :) = temp_E;
    end
end