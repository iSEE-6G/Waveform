
%%

%clc
clear all
%close all
tic
%% OTFS parameters%%%%%%%%%%
% number of symbol
N = 8;
% number of subcarriers
M = 8;
% size of constellation
M_mod = 4;
M_bits = log2(M_mod);
% average energy per data symbol
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));
% number of symbols per frame
N_syms_perfram = N*M;
% number of bits per frame
N_bits_perfram = N*M*M_bits;

SNR_dB = 30;
SNR = 10.^(SNR_dB/10);
noise_var_sqrt = sqrt(1./SNR);
sigma_2 = abs(eng_sqrt*noise_var_sqrt).^2;
%%
rng(1)
N_fram = 10^2;
err_ber = zeros(length(SNR_dB),1);

        %% random input bits generation%%%%%
        data_info_bit = randi([0,1],N_bits_perfram,1);
        data_temp = bi2de(reshape(data_info_bit,N_syms_perfram,M_bits));
        x = qammod(data_temp,M_mod,'gray');
        x = reshape(x,N,M);
        
        %% OTFS modulation%%%%
        s = OTFS_modulation(N,M,x);

        %% OTFS channel output%%%%%
        r = OTFS_channel_output(N,M,sigma_2,s);
        
        %% OTFS demodulation%%%%
        y = OTFS_demodulation(N,M,r);
        
        %% output bits and errors count%%%%%
        data_demapping = qamdemod(y,M_mod,'gray');
        data_info_est = reshape(de2bi(data_demapping,M_bits),N_bits_perfram,1);
        errors = sum(xor(data_info_est,data_info_bit));
        BER=errors/N_bits_perfram
        
