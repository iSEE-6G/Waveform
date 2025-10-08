clear
%close all
%% BER - Doppler Delay Diagram
SNR=0:2:20;
N_frame=10^2; % Number of frames
nPilots=4; % number of pilots per ofdm symbol
Nfft = 1024-nPilots;
Nsym = 6; % number of ofdm symbols

err_ber= zeros(size(SNR));

M=4;


%----4-----
k=1;
for i=SNR
    for jframe=1:N_frame
        err=ofdm(i,M, nPilots, Nfft,Nsym);
        err_ber(k) = err+err_ber(k)
        jframe
    end
    k=k+1;
end
err_ber_fram = err_ber/(Nfft*Nsym)./N_frame
semilogy(SNR,err_ber_fram,'-*','LineWidth',2)
grid

hold on

%}
%----16-----
err_ber= zeros(size(SNR));
k=1;
for i=SNR
    for jframe=1:N_frame
        err=ofdm(i,16, nPilots, Nfft,Nsym);
        err_ber(k) = err+err_ber(k)
        jframe
    end
    k=k+1;
end
err_ber_fram = err_ber/(Nfft*Nsym)./N_frame;
semilogy(SNR,err_ber_fram,'-*','LineWidth',2)
grid

hold on

%----32-----
err_ber= zeros(size(SNR));
k=1;
for i=SNR
    for jframe=1:N_frame
        err=ofdm(i,32, nPilots, Nfft,Nsym);
        err_ber(k) = err+err_ber(k)
        jframe
    end
    k=k+1;
end
err_ber_fram = err_ber/(Nfft*Nsym)./N_frame;
semilogy(SNR,err_ber_fram,'-*','LineWidth',2)
grid

hold on


%----64-----
err_ber= zeros(size(SNR));
k=1;
for i=SNR
    for jframe=1:N_frame
        err=ofdm(i,64, nPilots, Nfft,Nsym);
        err_ber(k) = err+err_ber(k)
        jframe
    end
    k=k+1;
end
err_ber_fram = err_ber/(Nfft*Nsym)./N_frame;
semilogy(SNR,err_ber_fram,'-*','LineWidth',2)
grid
%{
legend({'OFDM QPSK','OFDM 8-PSK','OFDM 16-PSK','OFDM 64-PSK'},'Location','southwest')


title("OFDM BER-SNR Diagram")
xlabel("SNR")
ylabel("BER")
%}