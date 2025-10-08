
function X = Generate_2D_data_grid(N,M)
x_vec=zeros(N*M,1);
X=reshape(x_vec,M,N);
X(8,8)=qammod(1,16);
X(8,24)=qammod(1,16);
X(24,8)=qammod(1,16);
X(24,24)=qammod(1,16);
end