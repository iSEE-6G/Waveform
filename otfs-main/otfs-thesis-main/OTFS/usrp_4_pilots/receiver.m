%clear all
% Receiver
hrx= comm.SDRuReceiver;
hrx.Platform="B210";
hrx.SerialNum = '3150302';
hrx.CenterFrequency = 850e6;
hrx.MasterClockRate = 30e6;
hrx.DecimationFactor = 4;
hrx.Gain = 30;

hrx.OutputDataType = 'double';
hrx.SamplesPerFrame = 1024;

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
N = 32;
% M: number of subcarriers in frequency
M = 32;
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
N = 32;
% M: number of subcarriers in frequency
M = 32;
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


%% Normalized DFT matrix
Fn=dftmtx(N);  % Generate the DFT matrix
Fn=Fn./norm(Fn);  % normalize the DFT matrix

%% Preamble processing
P = zeros(2*hrx.SamplesPerFrame,1);
for n = 1:2*hrx.SamplesPerFrame
    P(n) = sum ( Rx(n:n+N*M/2-1)...
        .*conj( Rx(n+N*M/2:n+N*M-1) ) );
end

[val, ind] = max(P);
sync_point = ind;
fprintf('Synchronization point: %d \n', sync_point);
% time sync
Rx = Rx(sync_point:sync_point+5*hrx.SamplesPerFrame);
%% Frequency offset estimation & Compensation
freq_error = -angle(val)/pi;

% compensation
Rx_corrected = Rx.*exp(-2*1i*pi*freq_error*(sync_point:sync_point+length(Rx)-1)/N*M);
Rx = Rx_corrected(length(preamble)+ gap+1:length(preamble)+ gap+M*N);


%% OTFS demodulation%%%%
Y_tilda=reshape(Rx,M,N);     
Y = Y_tilda*Fn;             
load("xtf.mat");
Ytf=ISFFT(N,M,Y);

        %est = SFFT(N,M,reshape(channel,N,M));
        
%         Ytf_corrected=Ytf./est;
%         % Delay Doppler
%         Y_new=SFFT(N,M,Ytf_corrected);

        figure(1)
bar3(abs(SFFT(N,M,est)))
title("Estimated channel")
        figure(2)
        bar3(abs(Y))
title("Received pilot")

