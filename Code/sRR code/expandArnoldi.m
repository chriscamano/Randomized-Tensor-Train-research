function [B,d] = expandArnoldi(A,B,m)
%EXPANDARNOLDI Summary of this function goes here
%   A : current krylov subspace
%   m : amount to expand Krylov subspace by  

%global nbit

d=size(B,2); 
H=zeros(d+m,d);%compute size of current Krylov space

for i=d:(d+m)-1 
    z = A*B(:,i); 
    %nbit = nbit + 1;
    delta = norm(z); % 
    %% orthogonalize
    h1 = B(:,1:i)'*z;
    z = z - B(:,1:i)*h1;
    %% reorthogonalize
    h2 = 0;
    if norm(z) < 0.5*delta    
        h2 = B(:,1:i)'*z;
        z = z - B(:,1:i)*h2;
    end
    %% update H
    H(1:i,i) = h1 + h2;
    H(i+1,i) = norm(z);
    %% expand subspace
    B(:,i+1) = z/H(i+1,i);
end;
d=d+m;                              %update size in main program

end

