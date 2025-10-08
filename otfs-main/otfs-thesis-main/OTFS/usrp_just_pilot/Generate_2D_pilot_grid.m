
function X = Generate_2D_pilot_grid(N,M)
x_vec=zeros(N*M,1);
X=reshape(x_vec,M,N);
X(16,16)=qammod(1,16);
end