%clear all
% Receiver
hrx= comm.SDRuReceiver;
hrx.Platform="B210";
hrx.SerialNum = '3150302';
hrx.CenterFrequency = 2.4e9;
hrx.MasterClockRate = 30e6;
hrx.DecimationFactor = 4;
hrx.Gain = 60;

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


% number of symbol
N = 32;
% number of subcarriers
M = 32;
% size of constellation
M_mod = 4;
M_bits = log2(M_mod);
% average energy per data symbol
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));
% number of symbols per frame
N_syms_perfram = N*M;
% number of bits per frame
N_bits_perfram = N*M*M_bits;



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
Rx = Rx(sync_point:sync_point-1+2*M*N+20);
%% Frequency offset estimation & Compensation
freq_error = -angle(val)/pi;

% compensation
Rx_corrected = Rx.*exp(-2*1i*pi*freq_error*(sync_point:sync_point+length(Rx)-1)/N*M);
Rx_corrected = Rx_corrected(N*M+1:end);

Rx=Rx_corrected(1,1:hrx.SamplesPerFrame);

%% OTFS demodulation%%%%
y = OTFS_demodulation(N,M,Rx');
        
%% output bits and errors count%%%%%
data_demapping = qamdemod(y,M_mod,'gray');
data_info_est = reshape(de2bi(data_demapping,M_bits),N_bits_perfram,1);
load('data_info_bit.mat')
errors = sum(xor(data_info_est,data_info_bit));

[size_data,~] = size(data_info_bit);
BER=errors/size_data
  