function [Id_cr, Vd_cr] = contact_resistance(Id_holes, Id_electrons, Vd, channel_length)

% channel-length in nanometers
% contact res originally in Ohm-um
contact_res_elec = (-0.13878*(channel_length/1000))+293.266;
contact_res_holes = (-0.11943*(channel_length/1000))+209.557;
damping = 0.2;

%Convergence Voltage Difference
epsilon = 1e-4;
for v_index=1:length(Vd)
    V0 = Vd(v_index);
    Vd_loss_0_elec = Id_electrons(v_index)*contact_res_elec;
    Vd_loss_0_holes = Id_holes(v_index)*contact_res_holes;
    new_V = V0-(Vd_loss_0_elec+Vd_loss_0_holes);
    while new_V - V0 > epsilon
        [i_hole,i_electrons] = find_id_for_vd(new_V);
        Vd_loss_elec = i_electrons*contact_res_elec;
        Vd_loss_holes = i_hole*contact_res_holes;
        new_V = (new_V*damping) + (new_V-(Vd_loss_elec+Vd_loss_holes));
    end
end