tic
%% OTFS parameters%%%%%%%%%%
% N: number of symbols in time
N = 64;
% M: number of subcarriers in frequency
M = 64;
% M_mod: size of QAM constellation
M_mod = 4;
M_bits = log2(M_mod);
% average energy per data symbol
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));
% OTFS variant: (RZP / RCP / CP / ZP)
variant='ZP';

%% Preamble generation:
preamble = pskmod(randi([0, M_mod-1], 1, M*N/2), M_mod,pi/M_mod);
preamble = (ifft(ifftshift(preamble), M*N/2)*sqrt(M*N/2));
preamble = repmat(preamble, 1, 2);

% average energy per data symbol
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));
% number of symbols per frame
N_syms_perfram = N*M;
% number of bits per frame
N_bits_perfram = N*M*M_bits;

SNR_dB = 10;
SNR = 10.^(SNR_dB/10);
noise_var_sqrt = sqrt(1./SNR);
sigma_2 = abs(eng_sqrt*noise_var_sqrt).^2;
%%
rng(1)
N_fram = 10^2;
err_ber = zeros(length(SNR_dB),1);

        %% random input bits generation%%%%%
        data_info_bit = randi([0,1],N_bits_perfram,1);
        data_temp = bi2de(reshape(data_info_bit,N_syms_perfram,M_bits));
        save('data_info_bit.mat')
        x = qammod(data_temp,M_mod,'gray');
        x = reshape(x,N,M);
        
        %% OTFS modulation%%%%
        s = OTFS_modulation(N,M,x);

        s = [ preamble, zeros(1,20), s' zeros(1,1000)];
        s = s/max([max(real(s), max(imag(s)))]);
        %% SDR
% Transmitter
htx = comm.SDRuTransmitter;
htx.Platform="B200";
htx.SerialNum = '315C79D';

htx.CenterFrequency = 2.4e9;
htx.Gain = 60;
htx.MasterClockRate = 30e6;
htx.InterpolationFactor =4;


while true 
    htx(s.');
end