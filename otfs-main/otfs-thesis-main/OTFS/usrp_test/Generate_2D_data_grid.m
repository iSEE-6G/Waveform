
function [X,pilot_pos] = Generate_2D_data_grid(N,M,x_data,data_grid)
    x_vec=zeros(N*M,1);
    data_array=reshape(data_grid,1,N*M);
    [~,data_pos]=find(data_array==1);
    x_vec(data_pos)=x_data;

    [~,pilot_pos]=find(data_array==2);
    x_vec(pilot_pos)=sqrt(N)*qammod(1,4);
    X=reshape(x_vec,M,N);


end