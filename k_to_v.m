function[vx] = k_to_v(k_x, k_y, delta);
%delta is normalized to lattice constant
E_plus_x = graphene_E_k(k_x+delta, k_y);
E_minus_x = graphene_E_k(k_x-delta, k_y);
vx = (E_plus_x - E_minus_x)/delta;