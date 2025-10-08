%clear all
%% Parameter declaration 

Ncp = 128; %number of cyclic prefix symbols
FreqOffset = 0.1;
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
data = randi([0 M-1],1,Nsym*(Nfft-nPilots-2*nSide-1)); %generate data
save('data')
qpsk_mod = pskmod(data,M, pi/M); %modulate data
save('qpsk_mod')

data_idx = setdiff(1:Nfft, pilotIdx); % index of data symbols
data_idx = setdiff(data_idx, 1:nSide);
data_idx = setdiff(data_idx, Nfft-nSide+1:Nfft);
data_idx = setdiff(data_idx, Nfft/2+1);
 

Tx = zeros(1,Nsym*(Nfft+Ncp)); 
OFDMsym = zeros(1,Nfft);
for sym = 1:Nsym
    OFDMsym(data_idx) = qpsk_mod((Nfft-nPilots-2*nSide-1)*(sym-1)+1:((Nfft-nPilots-2*nSide-1)*sym)); % add data to ofdm symbol
    OFDMsym(pilotIdx) = pilot_symbols; % add pilots to ofdm symbol
    OFDMsym = (ifft(ifftshift(OFDMsym))*sqrt(Nfft));
    % add ofdm symbols to a total array for transmit
    Tx((Nfft+Ncp)*(sym-1)+1:(Nfft+Ncp)*sym) = [OFDMsym(Nfft-Ncp+1:Nfft) OFDMsym]; 
end
Tx = [ preamble, Tx zeros(1,1000)];
Tx = Tx/max(real(abs(Tx)));

%% SDR
% Transmitter
htx = comm.SDRuTransmitter;
htx.Platform="B200";
htx.SerialNum = '315C79D';

htx.CenterFrequency = 2.5e9;
htx.Gain = 60;
htx.MasterClockRate = 30e6;
htx.InterpolationFactor =4;


while true 
    htx(Tx.');
end