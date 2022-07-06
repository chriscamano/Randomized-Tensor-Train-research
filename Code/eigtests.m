
% 
% %generate a random 1000x1000 symmetric adjaceny matrix
% N=1000;
% A = randn(N,N) ;
% %copy upper triangle to lower triangle
% B=tt_matrix(A - tril(A,-1) + triu(A,1)');
% 

n = 5;                                           % number of spins
H = HamHeis(n);                                  % Heisenberg Hamiltonian
tt_H = tt_matrix(H,eps,2*ones(1,n),2*ones(1,n))  % convert to TT format
norm(H - full(tt_H))



%Basic Tests
%____________________________________________________________________


%QR dmrg_eig (base case)
%____________________________________________________________________
disp('╔═════════════════════════╗')
disp( "        QR based DMRG eigen test")
disp('╚═════════════════════════╝')
[x,theta,out]=dmrg_eig_qr(tt_H,1e-6);

%Ensure that compute dmrg results are consistent by checking norm for
%Ax=lambdaX
fprintf(' \n━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
disp   ('        norm test')
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n')
norm(H*full(x))-norm(theta*full(x))

%SVD dmrg_eig (base case)
%____________________________________________________________________
disp('╔═════════════════════════╗')
disp( "        SVD based DMRG eigen test")
disp('╚═════════════════════════╝')
[x,theta,out]=dmrg_eig_svd(tt_H,1e-6);
%Ensure that compute dmrg results are consistent by checking norm for
%Ax=lambdaX
fprintf(' \n━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
disp   ('        norm test')
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n')
norm(H*full(x))-norm(theta*full(x))


%rSVD dmrg_eig (base case)
%____________________________________________________________________
disp('╔═════════════════════════╗')
disp( "        rSVD based DMRG eigen test")
disp('╚═════════════════════════╝')
[x,theta,out]=dmrg_eig_rsvd(tt_H,1e-6);
%Ensure that compute dmrg results are consistent by checking norm for
%Ax=lambdaX
fprintf(' \n━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
disp   ('        norm test')
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n')
norm(H*full(x))-norm(theta*full(x))
