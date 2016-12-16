function [Id_cr, Vd_cr] = contact_res(Id_holes, Id_electrons, Vd, channel_length)

% channel-length in nanometers
% contact res originally in Ohm-um
contact_res_elec = (-0.13878*(channel_length/1000))+293.266;
contact_res_holes = (-0.11943*(channel_length/1000))+209.557;
for v_index=1:length(Vd)
    Ih = Id_holes(v_index);
    Ie = Id_electrons(v_index);
    V = Vd(v_index);
    current_loss = -(V*contact_res_elec);
end