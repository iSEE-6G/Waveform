function [pilot_positions_subcarriers, pilot_positions_time, total_pilots_per_ru] = allocate_pilots(num_rus, num_subcarriers, num_symbols, pilot_step_freq, pilot_step_time)

if num_rus>pilot_step_time
    error('This allocation cannot be done if rus are more than pilot_step_time');
end

% pilot positions:
pilot_positions_subcarriers = 0:pilot_step_freq:num_subcarriers-1;
pilot_positions_subcarriers = repmat(pilot_positions_subcarriers, num_rus, 1);

pilot_positions_time = [];
for ii = 1 : num_rus
    tmp = (1+ii):pilot_step_time:num_symbols;
    pilot_positions_time = [pilot_positions_time; tmp];
end

total_pilots_per_ru = size(pilot_positions_time,2)*size(pilot_positions_subcarriers,2);
