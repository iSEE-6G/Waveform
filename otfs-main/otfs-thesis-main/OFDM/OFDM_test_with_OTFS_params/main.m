clear
%close all
%% BER - Doppler Delay Diagram
SNR=5:5:25;
N_frame=10^2; % Number of frames
nPilots=9; % number of pilots per ofdm symbol
Nfft = 1024-nPilots;
Nsym = 6; % number of ofdm symbols
Ncp = 128;
err_ber= zeros(size(SNR));

%----4-----
k=1;
for i=SNR
    for jframe=1:N_frame
        err=ofdm_calculate_error_with_doppler(i,4);
        err_ber(k) = err+err_ber(k)
        jframe
    end
    
    k=k+1;
end
err_ber_fram = err_ber/((Nfft+nPilots+Ncp)*Nsym)./N_frame;
semilogy(SNR,err_ber_fram,'-*','LineWidth',2)
grid

hold on

%{

%----16-----
k=1;
for i=SNR
    for jframe=1:N_frame
        err=ofdm_calculate_error_with_doppler(i,16, nPilots, Nfft,Nsym);
        err_ber(k) = err+err_ber(k);
        jframe
    end
    k=k+1;
end
err_ber_fram = err_ber/(Nfft*Nsym)./N_frame;
semilogy(SNR,err_ber_fram,'-*','LineWidth',2)
grid
hold on

%----64-----
k=1;
for i=SNR
    for jframe=1:N_frame
        err=ofdm_calculate_error_with_doppler(i,64, nPilots, Nfft,Nsym);
        err_ber(k) = err+err_ber(k);
        jframe
    end
    k=k+1;
end
err_ber_fram = err_ber/(Nfft*Nsym)./N_frame;
semilogy(SNR,err_ber_fram,'-*','LineWidth',2)
grid

legend({'OFDM 4QAM','OFDM 16QAM','OFDM 64QAM'},'Location','southwest')


title("OFDM BER-SNR Diagram")
xlabel("SNR")
ylabel("BER")
%}