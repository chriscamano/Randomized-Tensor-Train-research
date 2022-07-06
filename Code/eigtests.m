

%generate a random 1000x1000 symmetric adjaceny matrix
N=1000;
A = randn(N,N) ;
%copy upper triangle to lower triangle 
A = A - tril(A,-1) + triu(A,1)';

B=tt_matrix(A);


%Basic Tests
%____________________________________________________________________


%QR dmrg_eig (base case)
%____________________________________________________________________
[x,theta,out]=dmrg_eig_qr(B,1e-6)
%Ensure that compute dmrg results are consistent by checking norm for
%Ax=lambdaX
norm(A*full(x))-norm(theta*full(x))

%SVD dmrg_eig (base case)
%____________________________________________________________________

%rSVD dmrg_eig (base case)
%____________________________________________________________________

