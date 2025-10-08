%% Initializations:
freq = 28e9;
slot_dur = 0.125e-3;
num_scatterers = 10;
num_subcarriers = 3360;
num_slots = 240;
BW = 120e3 * num_subcarriers;
num_elements_tx = 1;
num_elements_rx = 4;

% mikos kymatos 
c = 3 * (10^8);
lambda = c/freq;

% Distance threshold between Tx and Rx
distance_thres = 200;

% Location of Rx - Tx is assumed at (0,0)
while 1
  x = 10 + rand(1,1) * 100;
  y = 10 + rand(1,1) * 100;
  d = sqrt(x^2 + y^2);
  if (d <= distance_thres)
      break;
  end
end



SNRdB = 23;
SNRmargin = 7;
SNR = 10^(SNRdB/10);





% Scatterer locations, distance from origin (Tx) and distance from Rx:
x_scat = zeros(num_scatterers,1);
y_scat = zeros(num_scatterers,1);
d_from_origin = zeros(num_scatterers,1);
d_from_destination = zeros(num_scatterers,1);
d_total = zeros(num_scatterers,1);

for kk = 1:num_scatterers
    while 1
        x_scat(kk) = 1 + rand(1,1) * (d-1)-1;
        y_scat(kk) = 1 + rand(1,1) * (d-1)-1;

        d_from_origin(kk) = sqrt(x_scat(kk)^2 + y_scat(kk)^2);
        d_from_destination(kk) = sqrt((abs(x_scat(kk)-x)^2 + abs(y_scat(kk)-y)^2));

        d_sum = d_from_destination(kk) + d_from_origin(kk);

        if (d_from_origin(kk) < d) && ~ismember(d_sum, d_total(1:kk-1)) && ~any(abs(d_total(1:kk-1)-d_sum)<=2)
            d_total(kk) = d_sum;
            break;
        end
    end
end

% Losses for the LoS:
P = (lambda / (4 * pi * d)) ^ 2;

% Losses for the other paths:
P_all = 0.4*(lambda ./ (4 * pi * (d_from_origin + d_from_destination) )).^ 2;

% Phase due to distance for LoS:
fash = 2*pi/lambda*d; 

% Phase due to distance for multipath:
fash_all = 2*pi/lambda*(d_from_origin + d_from_destination); 

% From path power to complex gain:
% First take the square root for gain:
a = sqrt(P); %LoS:
% Multipath:
a_all = sqrt(P_all);

% Then add the phase shift due to distance:
s = a * exp(1j*fash); %LoS
s_all = a_all .* exp(1j*fash_all); %multipath

% Create the random phase shift:
fi = 0 + (2 * pi) * rand(num_scatterers, 10);
fi = sum(fi, 2);

% Then add the random phase shift:
s_all = s_all .* exp(1j*fi); 

%% Doppler and movement:
% velocity for Tx and Rx (max 50 m/s)
u_MaxTR = 0;
u_Max = 70;
u_T = rand(1,1) * u_MaxTR; 
u_R = rand(1,1) * u_MaxTR; 

% Direction of movement:
f_T = 0 + (2 * pi) * rand(1,1);  
f_R = 0 + (2 * pi) * rand(1,1); 
% Relative direction of Rx-Tx movement
f_relative_TR = f_R - f_T;

% Velocity in each axis
u_Tx = u_T * cos(f_T); 
u_Ty = u_T * sin(f_T);
u_Rx = u_R * cos(f_R);
u_Ry = u_R * sin(f_R);

relative_velocity_x = u_Rx - u_Tx;
relative_velocity_y = u_Ry - u_Ty;
relative_velocity = sqrt((relative_velocity_x ^ 2) + (relative_velocity_y ^ 2)); % // sxetikh taxythta

%% H dieuthinsi tis sxetikhs taxutitas einai: f_relative_TR opote:
direct_freq_doppler = (relative_velocity / lambda) * cos(f_relative_TR);

% Movement of the scatterers:
u_scatterers = -u_Max + rand(num_scatterers,1) * 2*u_Max;
u_scatterers = -45:10:50;
f_scat = 2*pi*rand(num_scatterers,1);
% TWRA BRISKW TIS SXETIKES GWNIES KINISIS APO POMPO SE SKEDASTI:
f_relative_Tscat = f_scat - f_T;
% TWRA BRISKW TIS SXETIKES GWNIES KINISIS APO SKEDASTI SE POMPO:
f_relative_scatR = f_R - f_scat;

u_scatterers_x = u_scatterers(kk).*cos(f_scat);
u_scatterers_y = u_scatterers(kk).*sin(f_scat);
relative_scatterers_Tx = u_scatterers_x - u_Tx; % // sxetikh taxythta apo pompo se skedasth
relative_scatterers_Ty = u_scatterers_y - u_Ty;
relative_scatterers_Rx = u_Rx - u_scatterers_x; % // sxetikh taxythta apo skedasth se dekth
relative_scatterers_Ry = u_Ry - u_scatterers_y;

relative_velocity_scatterers_tx = sqrt((relative_scatterers_Tx.^ 2) + (relative_scatterers_Ty.^ 2)); 
relative_velocity_scatterers_rx = sqrt((relative_scatterers_Rx.^ 2) + (relative_scatterers_Ry.^ 2));

% KAI ANTISTOIXA DYO DOPPLER SHIFTS!: PROSOXI BAZOYME THN SXETIKI GWNIA (STO ALLAKSA):
fd_iT = (relative_velocity_scatterers_tx / lambda).* cos(f_relative_Tscat);
fd_iR = (relative_velocity_scatterers_rx / lambda).* cos(f_relative_scatR);

%% OFDM Parameters
FFT = 4096; 
cplen = 240;
Ts = num_subcarriers * 1 / BW; 
Tcp = cplen * 1 / BW; 
OFDM_symbol_time = Ts + Tcp; 
Tb = slot_dur * num_slots;
num_of_sensing_signals = 2 * 240;
num_ofdm_symbols_per_slot = 14;
mod_rank = 16; % {4, 16, 32, 64, ...}
bit_per_symbol = log2(mod_rank);
num_ofdm_symbols = 10;
total_ofdm_symbols = num_ofdm_symbols_per_slot * num_slots;

% pilot positions:
pilot_steps = 7;
pilot_positions_subcarriers = 0:pilot_steps:num_subcarriers-1;
pilot_positions_time = [2, 9];
total_pilots = length(pilot_positions_time)*length(pilot_positions_subcarriers);

bits_total = num_subcarriers * bit_per_symbol * total_ofdm_symbols;
data = randi([0, 1], 1, bits_total);

% create the table for OFDM symbols
symbol_matrix = zeros(num_subcarriers, num_ofdm_symbols_per_slot);
% create QPSK symbols for the pilots
for slot = pilot_positions_time
    symbol_matrix(pilot_positions_subcarriers+1, slot+1) = qammod(randi([0, 3], length(pilot_positions_subcarriers), 1), 4, 'InputType', 'integer');
end

symbol_matrix = repmat(symbol_matrix, 1, num_slots).';

% Time Axis:
timestep = 1/BW*(num_subcarriers+cplen); % equal to the time length of an ofdm symbol
time_vector=0:timestep:(total_ofdm_symbols-1)*timestep; % ORIZW ENA DIANYSMA XRONOY GIA POLY MIKRO DIASTIMA PX 10 EIS TIN MEION 6H SECOND 'H 10MICRO SEC

% TO DIRECT:
s = s * exp(2j*pi*direct_freq_doppler*time_vector);

% TA SCATTERERS:
s_all = s_all .* exp(2j*pi*fd_iT*time_vector) .* exp(2j*pi*fd_iR*time_vector);

% Delays:
delay_axis = 0: 1/BW: 10/BW;
delay_los = d / c;
delay = (d_from_destination + d_from_origin)/ c;
max_delay_length = max([delay_los; delay]);
delay1 = delay-delay_los;
delay_bin_size = 1/BW; 
num_delay_bins = ceil(max(delay) / delay_bin_size);
delay_axis = 0: delay_bin_size: (num_delay_bins-1)/BW;

gwnia_T = atan2(y, x);
gwnia_R = pi + gwnia_T;
gwnia_scat = atan2(y_scat, x_scat);
gwnia_scat_R = atan2(y_scat-y, x_scat-x);

d_antenna_tx = (0:num_elements_tx-1)*lambda/2;

steervec_los_tx = exp(-1j * (2 * pi * d_antenna_tx / lambda) * cos(gwnia_T));
steervec_scat_tx = zeros(num_scatterers, num_elements_tx);
for ii = 1 : num_scatterers
    steervec_scat_tx(ii,:) = exp(-1j * (2 * pi * d_antenna_tx / lambda) * cos(gwnia_scat(ii)));
end
steervec_tx = [steervec_los_tx; steervec_scat_tx];

d_antenna_rx = (0:num_elements_rx-1)*lambda/2;
steervec_los_rx = exp(-1j * (2 * pi * d_antenna_rx / lambda) * cos(gwnia_R));
steervec_scat_rx = zeros(num_scatterers, num_elements_rx);
for ii = 1 : num_scatterers
    steervec_scat_rx(ii,:) = exp(1j * (2 * pi * d_antenna_rx / lambda) * cos(gwnia_scat(ii)- gwnia_R));
end
steervec_rx = [steervec_los_rx; steervec_scat_rx];

% Initialize H with the additional dimension for delay bins
H = zeros(num_elements_rx, num_elements_tx, length(time_vector), num_delay_bins);

steervec_rx_H = steervec_rx';
delay = [delay_los; delay];

% Compute the channel matrix H for each time step and delay bin
for tt = 1:size(time_vector,2)
    s_all_t = s_all(:, tt);
    s_combined = [s(tt); s_all_t]; % LoS gain and scatterers gain
    
    low = 0;
    high = delay_bin_size;
    for delaybin = 1:num_delay_bins
        ind = (delay>=low) & (delay<high);
        if sum(ind)>0
            Sigma = diag(s_combined(ind));
            H(:, :, tt, delaybin) = steervec_rx_H(:,ind) * Sigma * steervec_tx(ind,:);
        end
        low = high;
        high = high + delay_bin_size;
    end
end

pilot_grid = zeros(num_subcarriers, num_slots*num_ofdm_symbols_per_slot);
for ii = 1 : num_slots
    pilot_grid(pilot_positions_subcarriers+1, (ii-1)*num_ofdm_symbols_per_slot + pilot_positions_time+1) = 1;
end
pilot_grid = pilot_grid';

%% Transmitter
qam = reshape(data, bit_per_symbol, bits_total/bit_per_symbol);
mean_qam = mean(abs(qammod(0:mod_rank-1, mod_rank)).^2);
qam = qammod(qam, mod_rank, 'InputType', 'bit');

serial_parallel = reshape(qam, num_subcarriers, []).';
serial_parallel(symbol_matrix~=0) = sqrt(mean_qam/2)*symbol_matrix(symbol_matrix~=0);
serial_parallel = serial_parallel/sqrt(mean_qam);


serial_parallel = [zeros(size(serial_parallel,1), (FFT-num_subcarriers)/2), serial_parallel, zeros(size(serial_parallel,1), (FFT-num_subcarriers)/2)]; % Add zero padding
signal = FFT/sqrt(num_subcarriers)*ifft(serial_parallel, FFT, 2);

add_cp = [signal(:, end-cplen+1:end) signal];
add_cp = add_cp.';

parallel_serial = add_cp(:);
x_tx = repmat(parallel_serial, 1, size(H,2));


%% CHANNEL:
% H (Rx antennas, Tx antennas, time, delay) 

sum_y_rx = zeros(length(parallel_serial)+ num_delay_bins-1, num_elements_rx);
sum_Y_rx = zeros(FFT, num_elements_rx, length(time_vector));
Hh = zeros(FFT, num_elements_rx, num_elements_tx, length(time_vector));
sum_Hh = zeros(FFT, num_elements_rx, length(time_vector));
snr_per_subcarrier = zeros(FFT, num_elements_rx, length(time_vector));

for tt = 1:length(time_vector)
    for rx_ant = 1 : size(H,1)
        for tx_ant = 1 : size(H,2)
            h = squeeze(H(rx_ant,tx_ant,tt,:));
            X = serial_parallel(tt,:);X=X.';
            Htmp = fft(h, FFT);
            Hh(:, rx_ant, tx_ant, tt) = Htmp; 
            % y_rx = conv(h, x_tx(:, tx_ant));
            Y = Htmp.*X;
            % sum_y_rx(:, rx_ant) = sum_y_rx(:, rx_ant) + y_rx;
            sum_Y_rx(:, rx_ant, tt) = Y;
            sum_Hh(:, rx_ant, tt) = sum_Hh(:, rx_ant, tt) + Hh(:, rx_ant, tx_ant, tt);
        end
        snr_per_subcarrier(:, rx_ant, tt) = abs(sum_Hh(:, rx_ant, tt)).^2;
    end
    display(tt)
end


% SUM the results y_rx1 y_rx2... y_rx4 (for rx ant 1) ... rx 

% sum_y_rx = sum(sum_y_rx,2);
% sum_Hh = sum(sum_Hh, 2);
% Exeis tesseris eksodo apo ton pompo -- analoga me tis keraies

%% Add noise
if num_elements_rx==1  
    signal_power = mean(squeeze(mean(abs(sum_Y_rx).^2, 1)));
else
    signal_power = mean(mean(squeeze(mean(abs(sum_Y_rx).^2, 1)),1));
end

sig_pwr_per_subc = mean(abs(sum_Y_rx).^2, 3);


noise_power = signal_power / SNR;
noise = sqrt(noise_power/2) * (randn(size(sum_Y_rx)) + 1i * randn(size(sum_Y_rx)));

snr_per_subcarrier = snr_per_subcarrier./noise_power;
snr_per_subcarrier = 10*log10(snr_per_subcarrier);

Rx_signal = sum_Y_rx + noise;
num_ofdm_symbols = size(serial_parallel,1);

%% Receiver
% reverse_parallel= zeros(FFT+cplen, size(serial_parallel,1), num_elements_rx);
% for rx_ant = 1:num_elements_rx`
%     tmp = reshape(rx_signal(1:end-num_delay_bins+1, rx_ant), FFT+cplen, []);
%     reverse_parallel(:,:,rx_ant) = tmp;
% end
% remove_cp = reverse_parallel(cplen + 1:end, :, :);

received_signal= zeros(FFT, size(serial_parallel,1), num_elements_rx);
for rx_ant = 1:num_elements_rx
    % received_signal(:,:,rx_ant) = sqrt(num_subcarriers)/FFT*fft(squeeze(remove_cp(:,:,rx_ant)), FFT, 1);
    received_signal(:,:,rx_ant) =  Rx_signal(:,rx_ant, :);
end

%% DRx:
DRx_big = squeeze(received_signal((FFT-num_subcarriers)/2+1:end-(FFT-num_subcarriers)/2, :, :));
DTx_big = serial_parallel(:, (FFT-num_subcarriers)/2+1:end-(FFT-num_subcarriers)/2);
DTx = zeros(num_of_sensing_signals);
DRx = zeros(num_of_sensing_signals,num_of_sensing_signals,  num_elements_rx);
counter = 1;
for ii = 3 : 7: num_ofdm_symbols
    DTx(counter, :) = DTx_big(ii, logical(pilot_grid(ii,:)));
    DRx(:, counter, :) = DRx_big(logical(pilot_grid(ii,:)), ii, :);
    counter = counter + 1;
end

DRx = permute(DRx, [2 1 3]);

%% HERE EQUALIZATION!!

x = zeros(size(received_signal));
for rx_ant = 1:num_elements_rx
    for sym = 1 : num_ofdm_symbols
        x(:, sym, rx_ant) = squeeze(received_signal(:, sym , rx_ant)) ./ squeeze(sum_Hh(:, rx_ant, sym));
    end
end

% range resolution, max range, velocity resolution, max velocity
range_resolution = c / (2 * BW); 
max_range = (c * FFT) / (2 * 120e3 * num_subcarriers); 
velocity_resolution = c / (2 * freq * Tb); 
max_velocity = (c * 120e3 * num_of_sensing_signals) / (2 * freq * total_ofdm_symbols); 

%%
D = zeros(size(DRx));
for ii = 1 : num_elements_rx
    D(:,:,ii) = DRx(:,:,ii)./DTx;
end

D_idft = ((ifft(D, [], 2)));% ifftshift(ifft(D, [], 2));
D_dft = fftshift(fft(D_idft, [], 1), 1); %fftshift(fft(D_idft, [], 1));

% D_idft = (ifft(D, [], 2));% ifftshift(ifft(D, [], 2));
% D_dft = (fft(D_idft, [], 1)); %fftshift(fft(D_idft, [], 1));

snr_thres = SNRdB - SNRmargin;
for rx_ant = 1 : num_elements_rx
    remove_first_subcarriers = squeeze(x((FFT-num_subcarriers)/2+1:end - (FFT-num_subcarriers)/2,:,rx_ant));
    tx_symbols = squeeze(serial_parallel(:, (FFT-num_subcarriers)/2+1:end - (FFT-num_subcarriers)/2)).';
    snr_per_sub = squeeze(snr_per_subcarrier((FFT-num_subcarriers)/2+1:end - (FFT-num_subcarriers)/2, rx_ant, :));
    ind = (snr_per_sub>snr_thres);
    ind2 = ind & pilot_grid;
    
    qam_demodulation = qamdemod(sqrt(mean_qam)*remove_first_subcarriers(ind2), mod_rank, 'OutputType', 'bit');
    qam_tx = qamdemod(sqrt(mean_qam)*tx_symbols(ind2), mod_rank, 'OutputType', 'bit');

    errors = sum(qam_demodulation~=qam_tx);
    error_rate = errors / sum(sum(ind2));
end

D = (D_dft(:, :, 1));
A = 10*log10(abs(D).^2);
A(A<max(max(A))-SNRdB*1.5) = -300;

% Apply Gaussian smoothing
A_smooth = imgaussfilt(A, 3); % Sigma = 1 (adjust as needed)

% Find local maxima
maximaMask = imregionalmax(A_smooth);
[row, col] = find(maximaMask);

% Display results
imagesc(A_smooth); hold on;
plot(col, row, 'r*', 'MarkerSize', 10); % Mark local maxima
colorbar;
title('Local Maxima After Gaussian Smoothing');

sorted_delay_real = sort(delay);
estimated = sort(FFT/num_subcarriers*col/BW);
f_est = sort(-1/slot_dur + (row-1)/slot_dur/num_slots);
sorted_fd = sort([fd_iT+fd_iR; 0]);

err = 0;
err1 = 0;
for ii = 1 : length(estimated)
    val = min(abs(sorted_delay_real - estimated(ii)));
    val1 = min(abs(f_est(ii) - sorted_fd));
    err = err + val;
    err1 = err1 + val1;
end
err = err/length(estimated);
err1 = err1/length(estimated);
dist_err = err*3e8;
freq_err = err1;
keyboard

real_num_of_paths = 1 + num_scatterers;
estimated_num_of_scatterers = length(estimated);

if real_num_of_paths >  estimated_num_of_scatterers
    prob_false_alarm = 0;
    prob_missed_detection = (real_num_of_paths - estimated_num_of_scatterers)/estimated_num_of_scatterers;
elseif real_num_of_paths <  estimated_num_of_scatterers
    prob_false_alarm = (estimated_num_of_scatterers-real_num_of_paths)/estimated_num_of_scatterers;
    prob_missed_detection = 0;
else
    prob_missed_detection = 0;
    prob_false_alarm = 0;
end