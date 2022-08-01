
function [B,d] = rBlockKrylov(A,B,del_p,bsize)
% This function is an incomplete blocked version of the k truncated arnoldi
% process that functions by using a blocked version of the krylov vector
% construction process. Issues currently are unidentified and most likley
% related to vectorization syntax errors. 
% 


%   A     : current krylov subspace
%   B     : Krylov subspace
%   del_p : amount to expand Krylov subspace by
%   bsize : size of krylov blocks
%global nbit

p=size(B,2)/bsize; %krylov depth
for i=p:(p+del_p)-1
    q=i*bsize;
    z = A*B(:,(i-1)*bsize+1:q);                                        %multiply A by ith block of B
    %delta = norm(z); %
    
    %% orthogonalize
    h1 = B(:,1:q)'*z; %h1
    z = z - B(:,1:q)*h1;
    
    %% reorthogonalize
    % if norm(z) < 0.5*delta
    h2 = B(:,1:q)'*z;
    z = z - B(:,1:q)*h2;
    %  end
    
    %% update H
    %H(1:q,(i-1)*bsize+1:q) = h1 + h2;
%     H(1:q,(i-1)*bsize+1:q) = h1;
%     H(i+bsize,1:q) = vecnorm(z);
    
    %% expand subspace
    B(:,q+1:q+bsize) = orth(z./vecnorm(z));
    
    %     end
end
d=p+del_p;                              %update size in main program

%% BAD!!!
B = orth(B);

end

