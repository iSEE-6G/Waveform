function [time_vector, delay_axis, num_delay_bins, delay_bin_size] = create_time_and_delay_axes(BW, ofdm_params, d, d_from_destination, d_from_origin)

c = 3e8; %speed of light
cplen = ofdm_params.cplen;
total_ofdm_symbols = ofdm_params.num_slots*ofdm_params.num_ofdm_symbols_per_slot;
num_subcarriers = ofdm_params.num_subcarriers;
num_scatterers = size(d_from_destination, 1);
num_ue = size(d_from_destination, 2);
num_ru = size(d_from_origin, 2);
delay_bin_size = 1/BW; 

% Time Axis:
timestep = 1/BW*(num_subcarriers+cplen); % equal to the time length of an ofdm symbol
time_vector=0:timestep:(total_ofdm_symbols-1)*timestep; % ORIZW ENA DIANYSMA XRONOY GIA POLY MIKRO DIASTIMA PX 10 EIS TIN MEION 6H SECOND 'H 10MICRO SEC

% Delays:
% delay_axis = 0: 1/BW: 10/BW;
num_delay_bins = 0;
delay_los = d / c;
for ii = 1 : num_ue
    d_from_dest_tmp = d_from_destination(:, ii);
    for ll = 1 : num_ru
        d_from_origin_tmp = d_from_origin(:, ll);
        delay = (d_from_dest_tmp + d_from_origin_tmp)/ c;

        tmp = ceil(max(delay) / delay_bin_size);
        num_delay_bins = max([num_delay_bins, tmp]);
    end
end

% max_delay_length = max([delay_los; delay]);
% delay1 = delay-delay_los;
delay_axis = 0: delay_bin_size: (num_delay_bins-1)/BW;
