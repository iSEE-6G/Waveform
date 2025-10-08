%% OTFS parameters%%%%%%%%%%
% N: number of symbols in time
N = 16;
% M: number of subcarriers in frequency
M = 64;
% M_mod: size of QAM constellation
M_mod = 4;
M_bits = log2(M_mod);
% average energy per data symbol
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));

%% delay-Doppler grid symbol placement
% ZP length  should be set to greater than or equal to maximum value
% of delay_taps
length_ZP = M/16;
% data positions of OTFS delay-Doppler domain data symbols  in the 2-D grid
M_data = M-length_ZP;
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

% SNR and variance of the noise
% SNR = P/\sigma^2; P: avg. power of albhabet transmitted
SNR_dB = 5:5:30;
SNR = 10.^(SNR_dB/10);
sigma_2 = (abs(eng_sqrt)^2)./SNR;
%% Initializing simulation error count variables

N_fram = 1000;

%% Normalized DFT matrix
Fn=dftmtx(N);  % Generate the DFT matrix
Fn=Fn./norm(Fn);  % normalize the DFT matrix
current_frame_number=zeros(1,length(SNR_dB));


trans_info_bit = randi([0,1],N_syms_perfram*M_bits,1);
        %%2D QAM symbols generation %%%%%%%%
        data=qammod(reshape(trans_info_bit,M_bits,N_syms_perfram), M_mod,'gray','InputType','bit');
        [X,data_pos] = Generate_2D_data_grid(N,M,data,data_grid);
        
        
        %% OTFS modulation%%%%
        X_tilda=X*Fn';
        s = reshape(X_tilda,N*M,1);

         %% OTFS demodulation%%%%
        Y_tilda=reshape(s,M,N);
        Y = Y_tilda*Fn;
        y=reshape(Y.',N*M,1);

        end_data=qamdemod(y(data_pos),M_mod,'gray','OutputType','bit');
        end_data-trans_info_bit
            