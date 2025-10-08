
function [H,H_tilda,P]= Gen_DD_and_DT_channel_matrices(N,M,G,Fn)
P=zeros(N*M,N*M);                  
for j=1:N
    for i=1:M
        E=zeros(M,N);
        E(i,j)=1;
        P((j-1)*M+1:j*M,(i-1)*N+1:i*N)=E;  % row-column interleaver matrix in equation (35) in [R1]
    end
end
H_tilda=(P'*G*P);  
H=kron(eye(M),Fn)*(P'*G*P)*kron(eye(M),Fn'); %using equation (33) and (40) in [R1]
end
