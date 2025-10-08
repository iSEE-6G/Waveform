initializations;
[x_ru, y_ru, x, y, d, d_from_origin, d_from_destination] = location_of_channel_nodes(bounds_ru, bounds, num_ue, num_ru, num_scatterers, distance_thres);
[s, P_all, faseis, fi_rand] = create_complex_gains(lambda, d, d_from_origin, d_from_destination, subpaths_per_cluster);

sensing_requirements;
scheduler_config;

[pilot_positions_subcarriers, pilot_positions_time, total_pilots_per_ru] = allocate_pilots(num_ru, num_subcarriers, num_ofdm_symbols_per_slot, pilot_step_freq, pilot_step_time);
[fd_iT, fd_iR, direct_freq_doppler] = doppler_creation(lambda, num_scatterers, u_MaxTR, u_Max, num_ru, num_ue);

[time_vector, delay_axis] = create_time_and_delay_axes(BW, ofdm_params, d ,d_from_destination, d_from_origin);

s_new = apply_doppler(time_vector, direct_freq_doppler, fd_iT, fd_iR, s);

create_ula_array(lambda, num_elements_tx, num_elements_rx, x_ru, y_ru, x, y, x_scat, y_scat)

ofdm_params.num_pilots = total_pilots_per_ru;
%% Assume Downlink -- for now:
data = create_data_for_users('full-buffer', num_ru, num_ue, ofdm_params, scheduler);

% if nargin<3

% end
