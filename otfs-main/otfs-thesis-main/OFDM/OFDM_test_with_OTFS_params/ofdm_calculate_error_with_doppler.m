function [tot] = ofdm_calculate_error_with_doppler(SNRdb,M)

%% Parameter declaration 
Ncp = 128; %number of cyclic prefix symbols
FreqOffset = 0;
%SNRdb = 15;
theta = 256;
Nfft = 1024;
Nsym = 6; % number of ofdm symbols

% average energy per data symbol
eng_sqrt = (M==2)+(M~=2)*sqrt((M-1)/6*(2^2));

nPilots=9; % number of pilots per ofdm symbol

pilotIdx=1:Nfft/(nPilots-1):Nfft;pilotIdx(nPilots) = Nfft;
pilot_symbols = exp(1i*pi/nPilots*(-nPilots:2:nPilots-1)); 
%M=4;

% Time and frequency resources
car_fre = 4*10^9;% Carrier frequency
delta_f = 15*10^3; % subcarrier spacing: 15 KHz
T = 1/delta_f; %one time symbol duration in OTFS frame

% SNR and variance of the noise
% SNR = P/\sigma^2; P: avg. power of albhabet transmitted
snr = 10.^(SNRdb/10);
sigma_2 = (abs(eng_sqrt)^2)./snr;


%% OFDM symbol generation
data = randi([0 M-1],1,Nsym*(Nfft-nPilots)); %generate data
%qpsk_mod = pskmod(data,M, pi/M); %modulate data
qpsk_mod = qammod(data,M); %modulate data

data_idx = setdiff(1:Nfft, pilotIdx); % index of data symbols

Tx = zeros(1,Nsym*(Nfft+Ncp+nPilots)); 

for sym = 1:Nsym
    OFDMsym(data_idx) = qpsk_mod((Nfft-nPilots)*(sym-1)+1:((Nfft-nPilots)*sym)); % add data to ofdm symbol
    OFDMsym(pilotIdx) = pilot_symbols; % add pilots to ofdm symbol
    OFDMsym = ifft(OFDMsym)*sqrt(Nfft);
    % add ofdm symbols to a total array for transmit
    Tx((Nfft+Ncp)*(sym-1)+1:(Nfft+Ncp)*sym) = [OFDMsym(Nfft-Ncp+1:Nfft) OFDMsym]; 
end

%% AWGN channel
% snr = 10^(-SNRdb/10);
% noise = sqrt(snr/2)*(randn(1,Nsym*(Nfft+Ncp+nPilots))+1i*randn(1,Nsym*(Nfft+Ncp+nPilots)));
% %Rx = exp(1i*2*pi*FreqOffset*(0:length(Tx)-1)./Nfft).*Tx+noise;
% Rx=Tx+noise;
max_speed=50;  % km/hr
        [chan_coef,delay_taps,Doppler_taps,taps]=Generate_delay_Doppler_channel_parameters(Nsym,Nfft,car_fre,delta_f,T,max_speed);
        L_set=unique(delay_taps);       
        gs=Gen_discrete_time_channel(Nsym,Nfft,taps,delay_taps,Doppler_taps,chan_coef);  %equation (14) in [R1]
        % Generate discrete-time baseband channel in TDL form (Eq. (2.22))

        %% channel output%%%%%             
        r=zeros(Nsym*Nfft,1);        
        l_max=max(delay_taps);
        for q=1:Nsym*Nfft
            for l=(L_set+1)
                if(q>=l)
                    r(q)=r(q)+gs(l,q)*Tx(q-l+1);  %equation (18) in [R1]
                end
            end
        end
        noise= sqrt(sigma_2/2)*(randn(size(Tx)) + 1i*randn(size(Tx)));
        Rx=r+noise;
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
%Rx_corrected = exp(-1i*2*pi*GM_sum(thetaEst(end))*(0:length(Rx)-1)./Nfft).*Rx; 
%Rx_n = zeros(1,Nsym*(Nfft));
Rx_corrected=Rx;
sym=1;
OFDMsym_r = Rx_corrected(Ncp+(Nfft+nPilots)*(sym-1)+1:(Ncp+Nfft+nPilots*sym)); % total ofdm-> (pilots+nfft+ncp)

symbols_r = OFDMsym_r(pilotIdx); %get pilots for ofdm symbol 1
channelEst(:,1) = symbols_r./pilot_symbols; % calculate channel estimation

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

%qpsk_demod=pskdemod(Rx_n,M,pi/M);
qpsk_demod=qamdemod(Rx_n,M);    

received_data=qpsk_demod;

tot = sum(xor(received_data,data));
