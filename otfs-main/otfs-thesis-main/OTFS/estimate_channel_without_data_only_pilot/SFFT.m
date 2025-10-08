
function Y = SFFT(N,M,Ytf)

Y = ifft(fft(Ytf).').'/sqrt(N/M); % SFFT
end