
function X = Generate_2D_data_grid(N,M)
x_vec=zeros(N*M,1);
X=reshape(x_vec,M,N);
X(12,12)=qammod(1,256);
end