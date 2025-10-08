rng('shuffle')
%% OTFS parameters%%%%%%%%%%
% N: number of symbols in time
N = 32;
% M: number of subcarriers in frequency
M = 32;
% M_mod: size of QAM constellation
M_mod = 4;
M_bits = log2(M_mod);

length_ZP = M/16; % ZP length (required only for CP-OTFS)
length_CP = 0;

M_data=M-length_ZP;
% data positions of OTFS delay-Doppler domain data symbols  in the 2-D grid
data_grid=zeros(M,N);
data_grid(1:M_data,1:N)=1;
% number of symbols per frame
N_syms_perfram = sum(sum(data_grid));
% number of bits per frame
N_bits_perfram = N_syms_perfram*M_bits;

% Time and frequency resources
car_fre = 4*10^9;% Carrier frequency
delta_f = 15*10^3; % subcarrier spacing: 15 KHz
T = 1/delta_f; %one time symbol duration in OTFS frame

SNR_dB=20;

%% Normalized DFT matrix
Fn=dftmtx(N);  % Generate the DFT matrix
Fn=Fn./norm(Fn);  % normalize the DFT matrix
current_frame_number=zeros(1,length(SNR_dB));

% Generate 32x32 Pilot grid
X = Generate_2D_pilot_grid(N,M);

Xtf=ISFFT(N,M,X);
save("xtf.mat","Xtf");
%% Preamble generation:
preamble = pskmod(randi([0, M_mod-1], 1, M*N/2), M_mod,pi/M_mod);
preamble = (ifft(ifftshift(preamble), M*N/2)*sqrt(M*N/2));
preamble = repmat(preamble, 1, 2);

%% OTFS modulation%%%%
X_tilda=X*Fn';                     %equation (2) in [R1]
s = reshape(X_tilda,N*M,1);        %equation (4) in [R1]

s = [ preamble, zeros(1,20), s' zeros(1,1000)];
s = s/max([max(real(s), max(imag(s)))]);

% Transmitter
htx = comm.SDRuTransmitter;
htx.Platform="B200";
htx.SerialNum = '315C79D';
    
htx.CenterFrequency = 850e6;
htx.Gain = 60;
htx.MasterClockRate = 30e6;
htx.InterpolationFactor =4;


while true 
    htx(s.');
end