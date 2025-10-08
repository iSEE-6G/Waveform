
function Xtf = ISFFT(N,M,x)
%% OTFS Modulation: 1. ISFFT, 2. Heisenberg transform
Xtf = fft(ifft(x).').'/sqrt(M/N); %%%ISFFT

end