%clear all
% Receiver
hrx= comm.SDRuReceiver;
hrx.Platform="B210";
hrx.SerialNum = '3150302';
hrx.CenterFrequency = 850e6;
hrx.MasterClockRate = 30e6;
hrx.DecimationFactor = 4;
hrx.Gain = 40;

hrx.OutputDataType = 'double';
hrx.SamplesPerFrame = 2*1024;

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
% number of bits per frame
N_bits_perfram = N_syms_perfram*M_bits;


%% Preamble generation:
preamble = pskmod(randi([0, M_mod-1], 1, M*N/2), M_mod,pi/M_mod);
preamble = (ifft(ifftshift(preamble), M*N/2)*sqrt(M*N/2));
preamble = repmat(preamble, 1, 2);

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

est=reshape(Rx_corrected(length(preamble)+ gap+1:length(preamble)+ gap+M*N),N,M);
est_tf=ISFFT(N,M,est);

H_tf=rot90(est_tf);
Rx = Rx_corrected(length(preamble)+ gap+1+M*N:length(preamble)+ gap+2*M*N);


%% OTFS demodulation%%%%
Y_tilda=reshape(Rx,M,N);     
Y = Y_tilda*Fn;             

% go to time frequency
%  Ytf=ISFFT(N,M,Y);
%  Ytf_corrected=Ytf./est;
%  
%  Y_new=SFFT(N,M,Ytf_corrected);
%  
  %data_array=reshape(data_grid,1,N*M);
  %[~,data_pos]=find(data_array==1);

%rec_data=reshape(qamdemod(Y(data_pos),M_mod,'gray','OutputType','bit'),[],1);

 %% MRC delay-time detection in [R1,R2,R3]
%% Different detection methods
        
        n_ite_MRC=15; % maximum number of MRC detector iterations  (Algorithm 2 in [R1])
        n_ite_algo3=15; % maximum number of matched filtered Guass Seidel (Algorithm 3 in [R1])
        n_ite_MPA=15; % maximum number of MPA detector iterations
        %damping parameter - reducing omega improves error performance at the cost of increased detector iterations
        omega=1;
        if(M_mod>=64)
            omega=0.25;     % set omega to a smaller value (for example: 0.05) for modulation orders greater than 64-QAM
        end
        decision=1; %1-hard decision, 0-soft decision
        init_estimate=1; %1-use the TF single tap estimate as the initial estimate for MRC detection, 0-initialize the symbol estimates to 0 at the start of MRC iteration
        %(Note: it is recommended to set init_estimate to 0 for higher order modulation schemes like 64-QAM as the single tap equalizer estimate may be less accurate)
[est_info_bits_1tap,data_1tap]=TF_single_tap_equalizer(N,M,M_mod,0,data_grid,Y,H_tf);

%% errors count%%%%%
load("t_data.mat")
errors = sum(xor(est_info_bits_1tap,trans_info_bit));
[s_er,~]=size(trans_info_bit);
ber=errors/s_er