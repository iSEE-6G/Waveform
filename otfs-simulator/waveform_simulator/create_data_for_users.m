function data = create_data_for_users(type, num_ru, num_ue, ofdm_params, scheduler)

%% Create random data:
bits_total = ofdm_params.num_subcarriers*ofdm_params.num_ofdm_symbols_per_slot*ofdm_params.num_slots*scheduler.bits_per_subcarrier;

data = cell(num_ru, 1);
if strcmp(type, 'full-buffer')

    for ii = 1 : num_ru

        data{ii} = randi([0, 1], 1, bits_total);
    end
end
