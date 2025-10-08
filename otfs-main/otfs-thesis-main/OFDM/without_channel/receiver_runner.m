clear 
clc
gain=0:5:50;

j=1;
for i=gain
    i
    for nFrame=1:3
        nFrame
        error_frames(nFrame) = ofdm_receiver(i);
    end
    error(j) = mean(error_frames);
    j=j+1;
end

semilogy(gain,error,'-or','linewidth',1.5); hold on 
xlabel('SNR [dB]')
ylabel('BER')
%title('OFDM-4QAM')
axis([0 65 10^(-4) 1])

