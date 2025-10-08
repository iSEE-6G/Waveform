function [tot] = ofdm(SNRdb,M, nPilots, Nfft, ...
    Nsym)

%% Parameter declaration 
Ncp = 128; %number of cyclic prefix symbols
FreqOffset = 0.15;
%SNRdb = 15;
theta = 256;

SNR = 10.^(-SNRdb/10);



pilotIdx=[150 400 900 1000];
pilot_symbols = exp(1i*pi/4*[-3 -1 1 3]); 

%% OFDM symbol generation
data = randi([0 M-1],1,Nsym*(Nfft)); %generate data
qpsk_mod = pskmod(data,M,pi/M); %modulate data

data_idx = setdiff(1:Nfft+nPilots, pilotIdx); % index of data symbols

Tx = zeros(1,Nsym*(Nfft+Ncp)); 

for sym = 1:Nsym
    OFDMsym(data_idx) = ifft(qpsk_mod(Nfft*(sym-1)+1:(Nfft*sym)),Nfft)*sqrt(Nfft); % add data to ofdm symbol
    OFDMsym(pilotIdx) = pilot_symbols(1:4); % add pilots to ofdm symbol

    % add ofdm symbols to a total array for transmit
    Tx((Nfft+Ncp+nPilots)*(sym-1)+1:(Nfft+Ncp+nPilots)*sym) = [OFDMsym(Nfft-Ncp+1:Nfft) OFDMsym]; 
end

%% AWGN channel
noise = sqrt(SNR/2)*(randn(1,Nsym*(Nfft+Ncp+nPilots))+1i*randn(1,Nsym*(Nfft+Ncp+nPilots)));
Tx_n=Tx+noise;
ricianChan = comm.RicianChannel('KFactor',60);

Rx=ricianChan(Tx_n');
Rx=Rx';
%Rx = exp(1i*2*pi*FreqOffset*(0:length(Tx)-1)./Nfft).*Tx+noise;

%% ML estimation of timing and frequency offset

PHI_sum = zeros(1,Nsym*(Nfft+Ncp)-Nfft);
GM_sum = zeros(1,Nsym*(Nfft+Ncp)-Nfft);

for n = theta:Nsym*(Nfft+Ncp)-(Nfft+Ncp)
    PHI=0;GM=0;
    for m = n:n+Ncp-1    
        PHI = PHI+ (Rx(m)*conj(Rx(m)) + Rx(m+Nfft)*conj(Rx(m+Nfft)));
        GM = GM+ Rx(m)*conj(Rx(m+Nfft));    
    end
    PHI_sum(n) = abs(GM)- (SNR/(SNR+1))*PHI;
    GM_sum(n) = -angle(GM)/(2*pi);
end

%% Estimation results display

%subplot(2,1,1);plot(PHI_sum);title('Estimation of timing offset');grid on;
%subplot(2,1,2);plot(GM_sum);title('Estimation of frequency offset');grid on;

[value,thetaEst]=findpeaks(PHI_sum,'minpeakdistance',Nfft);
%disp(thetaEst);
%disp(GM_sum(thetaEst));

%% Receiver
Rx_corrected = exp(-1i*2*pi*GM_sum(thetaEst(end))*(0:length(Rx)-1)./Nfft).*Rx; 
Rx_n = zeros(1,Nsym*(Nfft));
sym=1;
OFDMsym_r = Rx_corrected(Ncp+(Nfft+nPilots)*(sym-1)+1:(Ncp+Nfft+nPilots*sym)); % total ofdm-> (pilots+nfft+ncp)

symbols_r = OFDMsym_r(pilotIdx); %get pilots for ofdm symbol 1
channelEst(:,1) = symbols_r./pilot_symbols; % calculate channel estimation

% channel equalization for OFDM symbol 1
OFDMsym_r_data = (fft(OFDMsym_r(data_idx),Nfft)*sqrt(Nfft))./mean(channelEst(:,sym));
Rx_n((Nfft)*(sym-1)+1:(Nfft)*sym) = OFDMsym_r_data; % add data to total received data
for sym = 2:Nsym
    %total OFDM symbol
    OFDMsym_r = Rx_corrected((sym*Ncp+(Nfft+nPilots)*(sym-1)+1):(sym*Ncp+(Nfft+nPilots)*sym));

    symbols_r = OFDMsym_r(pilotIdx); % get pilots from received ofdm symbol
    channelEst(:,sym) = symbols_r./pilot_symbols; % do channel estimation
    OFDMsym_r_data = (fft(OFDMsym_r(data_idx),Nfft)*sqrt(Nfft))./mean(channelEst(:,sym)); % do channel equalization
    
    
    Rx_n((Nfft)*(sym-1)+1:(Nfft)*sym) = OFDMsym_r_data;% add data to total received data
end

qpsk_demod=pskdemod(Rx_n,M,pi/M);
received_data=qpsk_demod;

tot = sum(xor(received_data,data));
%tot=biterr(received_data, data)/(Nsym*(Nfft));
%fprintf("Totan bit Error rate: %.1f%%\n", tot*100);
end

