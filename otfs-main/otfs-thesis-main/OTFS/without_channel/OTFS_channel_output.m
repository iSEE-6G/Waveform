
%%
function r = OTFS_channel_output(N,M,sigma_2,s)
%% wireless channel and noise 

s = [s(N*M+1:N*M);s];%add one cp
noise = sqrt(sigma_2/2)*(randn(size(s)) + 1i*randn(size(s)));
r = s + noise;
r = r(1:(N*M));%discard cp
end