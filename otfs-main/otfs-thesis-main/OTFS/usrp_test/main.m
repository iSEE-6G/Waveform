rng('shuffle')
%% OTFS parameters%%%%%%%%%%
% N: number of symbols in time
N = 64;
% M: number of subcarriers in frequency
M = 64;
% M_mod: size of QAM constellation
M_mod = 4;
M_bits = log2(M_mod);
% average energy per data symbol
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));
% OTFS variant: (RZP / RCP / CP / ZP)
variant='ZP';
       
length_ZP = M/16; % ZP length (required only for CP-OTFS)
length_CP = 0;

M_data=M-length_ZP;
% data positions of OTFS delay-Doppler domain data symbols  in the 2-D grid
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

%% Normalized DFT matrix
Fn=dftmtx(N);  % Generate the DFT matrix
Fn=Fn./norm(Fn);  % normalize the DFT matrix

trans_info_bit = randi([0,1],N_syms_perfram*M_bits,1);

%% Preamble generation:
preamble = pskmod(randi([0, M_mod-1], 1, M*N/2), M_mod,pi/M_mod);
preamble = (ifft(ifftshift(preamble), M*N/2)*sqrt(M*N/2));
preamble = repmat(preamble, 1, 2);

%%2D QAM symbols generation %%%%%%%%        
data=qammod(reshape(trans_info_bit,M_bits,N_syms_perfram), M_mod,'gray','InputType','bit');
X = Generate_2D_data_grid(N,M,data,data_grid);
        

%% OTFS modulation%%%%
X_tilda=X*Fn';                     %equation (2) in [R1]
s = reshape(X_tilda,N*M,1);        %equation (4) in [R1]
ss = s;
gap = 20;
s = [ preamble, zeros(1,gap), zeros(1,4096) s.' zeros(1,1000)];
%s = s/max([max(real(s), max(imag(s)))]);

%% chan
Rx=[zeros(1,4096) s zeros(1,2*4096)];
%% Preamble processing
P = zeros(2*4096,1);
for n = 1:2*4096
    P(n) = sum ( Rx(n:n+N*M/2-1)...
        .*conj( Rx(n+N*M/2:n+N*M-1) ) );
end

[val, ind] = max(P);
sync_point = ind;
fprintf('Synchronization point: %d \n', sync_point);
% time sync
Rx = Rx(sync_point:end);
%% Frequency offset estimation & Compensation
freq_error = -angle(val)/pi;

% compensation
Rx_corrected = Rx.*exp(-2*1i*pi*freq_error*(sync_point:sync_point+length(Rx)-1)/N*M);

zeros_rx = Rx_corrected(length(preamble)+ gap+1:length(preamble)+ gap+M*N);
Rx = Rx_corrected(length(preamble)+ gap+1+M*N:length(preamble)+ gap+2*M*N);

% Rx=Rx_corrected(1,1:4096);

%% OTFS demodulation%%%%
Y_tilda=reshape(Rx,M,N);
Y = Y_tilda*Fn;

data_r = reshape(Y,[],1);
data_array=reshape(data_grid,1,N*M);
[~,data_pos]=find(data_array>0);

data_r=data_r(data_pos);


received_bits=reshape(qamdemod(data_r,M_mod,'gray','OutputType','bit'),[],1);

%% errors count%%%%%
errors = sum(xor(received_bits,trans_info_bit));
[s_er,~]=size(received_bits);
ber=errors/s_er
