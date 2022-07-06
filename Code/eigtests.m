

%generate a random 1000x1000 symmetric adjaceny matrix
N=1000;
A = randn(N,N) ;
%copy upper triangle to lower triangle 
A = A - tril(A,-1) + triu(A,1)';

B=tt_matrix(A);
B.d

%Basic Tests
%____________________________________________________________________


%QR dmrg_eig (base case)
%____________________________________________________________________
disp('╔═════════════════════════╗')
disp( "        QR based DMRG eigen test")
disp('╚═════════════════════════╝')
[x,theta,out]=dmrg_eig_qr(B,1e-6);
%Ensure that compute dmrg results are consistent by checking norm for
%Ax=lambdaX
fprintf(' \n━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
disp   ('        norm test')
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n')
norm(A*full(x))-norm(theta*full(x))

%SVD dmrg_eig (base case)
%____________________________________________________________________
disp('╔═════════════════════════╗')
disp( "        SVD based DMRG eigen test")
disp('╚═════════════════════════╝')
[x,theta,out]=dmrg_eig_svd(B,1e-6);
%Ensure that compute dmrg results are consistent by checking norm for
%Ax=lambdaX
fprintf(' \n━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
disp   ('        norm test')
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n')
norm(A*full(x))-norm(theta*full(x))


%rSVD dmrg_eig (base case)
%____________________________________________________________________
disp('╔═════════════════════════╗')
disp( "        rSVD based DMRG eigen test")
disp('╚═════════════════════════╝')
[x,theta,out]=dmrg_eig_rsvd(B,1e-6);
%Ensure that compute dmrg results are consistent by checking norm for
%Ax=lambdaX
fprintf(' \n━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
disp   ('        norm test')
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n')
norm(A*full(x))-norm(theta*full(x))
