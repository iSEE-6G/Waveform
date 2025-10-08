%clear all
% Receiver
hrx= comm.SDRuReceiver;
hrx.Platform="B210";
hrx.SerialNum = '3150302';
hrx.CenterFrequency = 850e6;
hrx.MasterClockRate = 30e6;
hrx.DecimationFactor = 4;
hrx.Gain = 5;

hrx.OutputDataType = 'double';
hrx.SamplesPerFrame = 2*4096;

rxLog = dsp.SignalSink;
for i=1:1000
    Rx = hrx(); 
    rxLog(Rx);    
end
Rx = rxLog.Buffer;
release(rxLog);
release(hrx);

rng('shuffle')
gap=20;
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
%% Preamble generation:
preamble = pskmod(randi([0, M_mod-1], 1, M*N/2), M_mod,pi/M_mod);
preamble = (ifft(ifftshift(preamble), M*N/2)*sqrt(M*N/2));
preamble = repmat(preamble, 1, 2);
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

trans_info_bit = randi([0,1],N_syms_perfram*M_bits,1);

%%2D QAM symbols generation %%%%%%%%        
data=qammod(reshape(trans_info_bit,M_bits,N_syms_perfram), M_mod,'gray','InputType','bit');
X = Generate_2D_data_grid(N,M,data,data_grid);

%% Preamble processing
%% Preamble processing
P = zeros(2*4096,1);
for n = 1:2*4096
    P(n) = sum ( Rx(n:n+N*M/2-1)...
        .*conj( Rx(n+N*M/2:n+N*M-1) ) );
end

[val, ind] = max(P);
sync_point = ind;
fprintf('Synchronization point: %d \n', sync_point);
% time sync
Rx = Rx(sync_point:sync_point+5*4096);
%% Frequency offset estimation & Compensation
freq_error = -angle(val)/pi;

% compensation
Rx_corrected = Rx.*exp(-2*1i*pi*freq_error*(sync_point:sync_point+length(Rx)-1)/N*M);

Rx = Rx_corrected(length(preamble)+ gap+1:length(preamble)+ gap+M*N);

H_tf = reshape(SFFT(N,M,zeros_rx), N,M);

%% OTFS demodulation%%%%
Y_tilda=reshape(Rx,M,N);     
Y = Y_tilda*Fn;             

 %% MRC delay-time detection in [R1,R2,R3]
        
[est_info_bits_1tap,data_1tap] = TF_single_tap_equalizer(N,M,M_mod,20,data_grid,Y,H_tf);

%% errors count%%%%%
load("t_data.mat")
errors = sum(xor(est_info_bits_1tap,trans_info_bit));
[s_er,~]=size(est_info_bits_1tap);
ber=errors/s_er