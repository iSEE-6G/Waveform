clear all
%% Parameter declaration 
Ncp = 128; % number of cyclic prefix symbols
FreqOffset = 0.1;
theta = 256;
Nfft = 1024;
Nsym = 6; % number of ofdm symbols

% pilots
nPilots=9; % number of pilots per ofdm symbol
pilotIdx=1:Nfft/(nPilots-1):Nfft;pilotIdx(nPilots) = Nfft;
pilot_symbols = exp(1i*pi/nPilots*(-nPilots:2:nPilots-1)); 

M=4; 

SNRdb=20;

%% Preamble generation:
preamble = pskmod(randi([0, M-1], 1, Nfft/2), M,pi/M);
preamble = ifft(preamble, Nfft/2)*sqrt(Nfft/2);
preamble = repmat(preamble, 1, 2);

%% OFDM symbol generation
data = randi([0 M-1],1,Nsym*(Nfft-nPilots)); %generate data
qpsk_mod = pskmod(data,M, pi/M); %modulate data

data_idx = setdiff(1:Nfft, pilotIdx); % index of data symbols

Tx = zeros(1,Nsym*(Nfft+Ncp)); 

for sym = 1:Nsym
    OFDMsym(data_idx) = qpsk_mod((Nfft-nPilots)*(sym-1)+1:((Nfft-nPilots)*sym)); % add data to ofdm symbol
    OFDMsym(pilotIdx) = pilot_symbols; % add pilots to ofdm symbol
    OFDMsym = ifft(OFDMsym)*sqrt(Nfft);
    % add ofdm symbols to a total array for transmit
    Tx((Nfft+Ncp)*(sym-1)+1:(Nfft+Ncp)*sym) = [OFDMsym(Nfft-Ncp+1:Nfft) OFDMsym]; 
end
Tx = [ preamble, Tx];
Tx = Tx/max(real(abs(Tx)));




%% AWGN channel
snr = 10^(-SNRdb/10);
noise = sqrt(snr/2)*(randn(1,length(Tx))+1i*randn(1,length(Tx)));
Rx = exp(1i*2*pi*FreqOffset*(0:length(Tx)-1)./Nfft).*Tx+noise;
% Rx=Tx;
%% ML estimation of timing and frequency offset
PHI_sum = zeros(1,Nsym*(Nfft+Ncp)-Nfft);
GM_sum = zeros(1,Nsym*(Nfft+Ncp)-Nfft);

for n = theta:Nsym*(Nfft+Ncp)-(Nfft+Ncp)
    PHI=0;GM=0;
    for m = n:n+Ncp-1    
        PHI = PHI+ (Rx(m)*conj(Rx(m)) + Rx(m+Nfft)*conj(Rx(m+Nfft)));
        GM = GM+ Rx(m)*conj(Rx(m+Nfft));    
    end
    PHI_sum(n) = abs(GM)- (snr/(snr+1))*PHI;
    GM_sum(n) = -angle(GM)/(2*pi);
end

%% Estimation results display

%subplot(2,1,1);plot(PHI_sum);title('Estimation of timing offset');grid on;
%subplot(2,1,2);plot(GM_sum);title('Estimation of frequency offset');grid on;

[value,thetaEst]=findpeaks(PHI_sum,'minpeakdistance',Nfft);
%disp(thetaEst);
%disp(GM_sum(thetaEst));

%% Receiver
%Rx_n = zeros(1,Nsym*(Nfft));
sym=1;
%% Preamble processing
P = zeros(length(Tx)-Nfft,1);
for n = 1:length(P)
    P(n) = sum ( Rx(n:n+Nfft/2-1)...
        .*conj( Rx(n+Nfft/2:n+Nfft-1) ) );
end
 
[val, ind] = max(P);
sync_point = ind;
fprintf('Synchronization point: %d \n', sync_point);
% time sync
Rx = Rx(sync_point:end);
 
%% Frequency offset estimation & Compensation
freq_error = -angle(val)/pi+0.05;
% compensation
Rx_corrected = Rx.*exp(-2*1i*pi*freq_error*(sync_point:sync_point+length(Rx)-1)/Nfft);
Rx_corrected = Rx_corrected(Nfft+1:end);
OFDMsym_r = Rx_corrected(Ncp+(Nfft)*(sym-1)+1:(Ncp+Nfft)*sym); % total ofdm-> (pilots+nfft+ncp)
for sym = 1:Nsym
    %total OFDM symbol
    OFDMsym_r = Rx_corrected((sym*Ncp+(Nfft)*(sym-1)+1):(sym*Ncp+(Nfft)*sym));
    OFDMsym_r_freq = (fft(OFDMsym_r,Nfft)/sqrt(Nfft));
    symbols_r = OFDMsym_r_freq(pilotIdx); % get pilots from received ofdm symbol
    channelEst(:,sym) = symbols_r./pilot_symbols; % do channel estimation
    interp_channel(:, sym) = interp1(pilotIdx, channelEst(:,sym), 1:Nfft, "pchip");
    OFDMsym_r_data = OFDMsym_r_freq./interp_channel(:, sym).'; % do channel equalization
   
    Rx_n((Nfft-nPilots)*(sym-1)+1:(Nfft-nPilots)*sym) = OFDMsym_r_data(data_idx);% add data to total received data
end

qpsk_demod=pskdemod(Rx_n,M,pi/M);
received_data=qpsk_demod;

total_error = sum(xor(received_data,data))