clear
%close all
%% BER - Doppler Delay Diagram
SNR=0:2:20;
N_frame=10^2; % Number of frames
nPilots=4; % number of pilots per ofdm symbol
Nfft = 1024-nPilots;
Nsym = 6; % number of ofdm symbols

err_ber= zeros(size(SNR));

% create rayleigh channel
rayleighchan = comm.RayleighChannel( ...
    'SampleRate',10e3, ...
    'PathDelays',[0 1.5e-4], ...
    'AveragePathGains',[2 3], ...
    'NormalizePathGains',true, ...
    'MaximumDopplerShift',30, ...
    'DopplerSpectrum',{doppler('Gaussian',0.6),doppler('Flat')}, ...
    'RandomStream','mt19937ar with seed', ...
    'Seed',22, ...
    'PathGainsOutputPort',true);
rayleighchan.RandomStream = 'Global stream';
rng(22)

%----4-----
k=1;
for i=SNR
    for jframe=1:N_frame
        err=ofdm_calculate_error_with_doppler(i,4, nPilots, Nfft,Nsym, rayleighchan);
        err_ber(k) = err+err_ber(k)
        jframe
    end
    k=k+1;
end
err_ber_fram = err_ber/(Nfft*Nsym)./N_frame
semilogy(SNR,err_ber_fram,'-*','LineWidth',2)
grid

hold on



%----16-----
k=1;
for i=SNR
    for jframe=1:N_frame
        err=ofdm_calculate_error_with_doppler(i,16, nPilots, Nfft,Nsym);
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
k=1;
for i=SNR
    for jframe=1:N_frame
        err=ofdm_calculate_error_with_doppler(i,64, nPilots, Nfft,Nsym);
        err_ber(k) = err+err_ber(k)
        jframe
    end
    k=k+1;
end
err_ber_fram = err_ber/(Nfft*Nsym)./N_frame;
semilogy(SNR,err_ber_fram,'-*','LineWidth',2)
grid
%{
legend({'OFDM 4QAM','OFDM 16QAM','OFDM 64QAM'},'Location','southwest')


title("OFDM BER-SNR Diagram")
xlabel("SNR")
ylabel("BER")
%}