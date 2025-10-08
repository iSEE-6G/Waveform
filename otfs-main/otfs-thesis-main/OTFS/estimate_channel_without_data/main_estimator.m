
clear all
%% OTFS parameters%%%%%%%%%%
% N: number of symbols in time
N = 64;
% M: number of subcarriers in frequency
M = 64;
% M_mod: size of QAM constellation
M_mod = 64;
M_bits = log2(M_mod);
% average energy per data symbol
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));
variant='ZP';

length_ZP = M-1; % ZP length (required only for ZP-OTFS)
    length_CP = 0;

    M_data=M-length_ZP;
% data positions of OTFS delay-Doppler domain data symbols  in the 2-D grid
data_grid=zeros(M,N);
data_grid(12,12)=1;
% number of symbols per frame
N_syms_perfram = sum(sum(data_grid));
% number of bits per frame
N_bits_perfram = N_syms_perfram*M_bits;


% Time and frequency resources
car_fre = 4*10^9;% Carrier frequency
delta_f = 10*10^3; % subcarrier spacing: 15 KHz
T = 1/delta_f; %one time symbol duration in OTFS frame

% SNR and variance of the noise
% SNR = P/\sigma^2; P: avg. power of albhabet transmitted
SNR_dB = 45;
SNR = 10.^(SNR_dB/10);
sigma_2 = (abs(eng_sqrt)^2)./SNR;

%% Normalized DFT matrix
Fn=dftmtx(N);  % Generate the DFT matrix
Fn=Fn./norm(Fn);  % normalize the DFT matrix
iFm = dftmtx(M)';
iFm = iFm/norm(iFm); 

trans_info_bit = randi([0,1],N_syms_perfram*M_bits,1);
%%2D QAM symbols generation %%%%%%%%
X = Generate_2D_data_grid(N,M);

%% OTFS modulation%%%%
        X_tilda=X*Fn';
        s = reshape(X_tilda,N*M,1);
        
        %% OTFS channel generation%%%% The user can use either the synthetic channel model or the 3GPP channel by uncommenting the corresonding piece of code.        
        %% synthetic channel model with equal power paths with delays [0,l_max] and Dopplers [-k_max,k_max]
%                 taps=4;
%                 l_max=delay_spread;
%                 k_max=4;
%                 chan_coef=1/sqrt(2)*(randn(1,taps)+1i.*randn(1,taps));
%                 delay_taps=randi(taps,[1,l_max+1]);  
%                 delay_taps=sort(delay_taps-min(delay_taps));  %% random delay shifts in the range [0,l_max] 
%                 Doppler_taps=k_max-2*k_max*rand(1,taps);   %% uniform Doppler profile [-k_max,k_max]
        
        % channel model following 3GPP standard
        max_speed=50;  % km/hr
        [chan_coef,delay_taps,Doppler_taps,taps]=Generate_delay_Doppler_channel_parameters(N,M,car_fre,delta_f,T,max_speed);
        
        chan_coef = [1 -0.5i 0.25];
        delay_taps = [0 1 2];
        Doppler_taps = [1 3 -3];
        
        L_set=unique(delay_taps);
        l_max=max(L_set);
        
        %% channel output%%%%%                      
        r=zeros(N*M,1);
        sigma_2 = (mean(abs(s).^2))./SNR;
        noise = sqrt(sigma_2/2)*(randn(size(s)) + 1i*randn(size(s))); 
        gs = Gen_time_domain_channel_OTFSvariants(N,M,delay_taps,Doppler_taps,chan_coef,length_CP,variant);
        for q=0:N*M-1  
            for l=L_set
                if(q>=l)
                    r(q+1)=r(q+1)+gs(l+1,q+1)*s(q-l+1);    % see Section 4.5.4, Chapter 4, [R1]
                end
            end            
        end
        r=r+noise;

%% OTFS demodulation%%%%
        Y_tilda=reshape(r,M,N);
        Y = Y_tilda*Fn;
        Y2 = (iFm*Y_tilda)*Fn;
%         Y = Y_tilda; 
        y=reshape(Y.',N*M,1);
        
         subplot(1,2,1);
         bar3(abs(X))
         subplot(1,2,2);
        bar3(abs(Y))