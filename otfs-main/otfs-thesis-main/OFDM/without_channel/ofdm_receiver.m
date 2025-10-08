%function [BER,SNR] = ofdm_receiver() 
%clear all
% Receiver
hrx= comm.SDRuReceiver;
hrx.Platform="B210";
hrx.SerialNum = '3150302';
hrx.CenterFrequency = 2.5e9;
hrx.MasterClockRate = 30e6;
hrx.DecimationFactor = 4;
hrx.Gain = 40;

hrx.OutputDataType = 'double';
hrx.SamplesPerFrame = 8936;

rxLog = dsp.SignalSink;
for i=1:1000
    Rx = hrx(); 
    rxLog(Rx);    
end

%% Parameter declaration 
Ncp = 128; %number of cyclic prefix symbols
%SNRdb = 15;
theta = 256;
Nfft = 1024;
Nsym = 6; % number of ofdm symbols

nPilots=9; % number of pilots per ofdm symbol
nSide = 32;

pilotIdx=nSide+1:(Nfft-nSide)/(nPilots-1):Nfft-nSide;
pilotIdx(nPilots) = Nfft-nSide;
pilot_symbols = exp(1i*pi/nPilots*(-nPilots:2:nPilots-1)); 
M=4;

SNRdb=20;

%% Preamble generation:
preamble = pskmod(randi([0, M-1], 1, Nfft/2), M,pi/M);
preamble = (ifft(ifftshift(preamble), Nfft/2)*sqrt(Nfft/2));
preamble = repmat(preamble, 1, 2);

%% OFDM symbol generation


data_idx = setdiff(1:Nfft, pilotIdx); % index of data symbols
data_idx = setdiff(data_idx, 1:nSide);
data_idx = setdiff(data_idx, Nfft-nSide+1:Nfft);
data_idx = setdiff(data_idx, Nfft/2+1);


Rx = rxLog.Buffer;
release(rxLog);
release(hrx);

sym = 1;
%% Preamble processing
P = zeros(2*hrx.SamplesPerFrame,1);
for n = 1:2*hrx.SamplesPerFrame
    P(n) = sum ( Rx(n:n+Nfft/2-1)...
        .*conj( Rx(n+Nfft/2:n+Nfft-1) ) );
end

[val, ind] = max(P);
sync_point = ind;
fprintf('Synchronization point: %d \n', sync_point);
% time sync
Rx = Rx(sync_point:sync_point-1+8936);
 
%% Frequency offset estimation & Compensation
freq_error = -angle(val)/pi;

% compensation
Rx_corrected = Rx.*exp(-2*1i*pi*freq_error*(sync_point:sync_point+length(Rx)-1)/Nfft);
Rx_corrected = Rx_corrected(Nfft+1:end);
OFDMsym_r = Rx_corrected(Ncp+(Nfft)*(sym-1)+1:(Ncp+Nfft)*sym); % total ofdm-> (pilots+nfft+ncp)
for sym = 1:Nsym
    %total OFDM symbol
    OFDMsym_r = Rx_corrected((sym*Ncp+(Nfft)*(sym-1)+1):(sym*Ncp+(Nfft)*sym));
    OFDMsym_r_freq = fftshift(fft(OFDMsym_r,Nfft)/sqrt(Nfft));
    symbols_r = OFDMsym_r_freq(pilotIdx); % get pilots from received ofdm symbol
    channelEst(:,sym) = symbols_r./pilot_symbols; % do channel estimation
    interp_channel(:, sym) = interp1(pilotIdx, channelEst(:,sym), 1:Nfft, "pchip");
    OFDMsym_r_data = OFDMsym_r_freq./interp_channel(:, sym).'; % do channel equalization
    
    %keyboard
    Rx_n((Nfft-nPilots-2*nSide-1)*(sym-1)+1:(Nfft-nPilots-2*nSide-1)*sym) = OFDMsym_r_data(data_idx);% add data to total received data
end

qpsk_demod=pskdemod(Rx_n,M,pi/M);
load('qpsk_mod.mat')
%qpsk_demod=qamdemod(Rx_n,M);
received_data=qpsk_demod;
load('data.mat');
total_error = sum(xor(received_data,data));
[~, ndata] = size(data);

BER = total_error/ndata

evm = lteEVM(Rx_n,qpsk_mod);
SNR = 10*log10(1/evm.RMS)
%end