
function X = Generate_2D_pilot_grid(N,M)
x_vec=zeros(N*M,1);
X=reshape(x_vec,M,N);
X(32,32)=qammod(1,256);
end