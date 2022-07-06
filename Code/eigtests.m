
% 
% %generate a random 1000x1000 symmetric adjaceny matrix
% N=1000;
% A = randn(N,N) ;
% %copy upper triangle to lower triangle
% B=tt_matrix(A - tril(A,-1) + triu(A,1)');
% 

n = 10;                                           % number of spins
H = HamHeis(n);                                  % Heisenberg Hamiltonian
tt_H = tt_matrix(H,eps,2*ones(1,n),2*ones(1,n))  % convert to TT format
norm(H - full(tt_H));


%compute average time
qr_time=0;
svd_time=0;
rsvd_time=0;

% for q=1:200
% disp(q);
% n = 10;                                           % number of spins
% H = HamHeis(n);                                  % Heisenberg Hamiltonian
% tt_H = tt_matrix(H,eps,2*ones(1,n),2*ones(1,n));  % convert to TT format
% norm(H - full(tt_H));
% 
% tic;
% [x,theta,out]=dmrg_eig_qr(tt_H,1e-6);
% qr_time=qr_time+toc;  
% tic;
% [x,theta,out]=dmrg_eig_svd(tt_H,1e-6);
% svd_time=svd_time+toc;
% 
% tic;
% [x,theta,out]=dmrg_eig_rsvd(tt_H,1e-6);
% rsvd_time=rsvd_time+toc;
% end
% 
% qr_time=qr_time/1000
% svd_time=svd_time/1000
% rsvd_time=rsvd_time/1000
% 
% plot([qr_time svd_time rsvd_time],'ro')


%Basic Tests
%____________________________________________________________________
%[x,theta,out]=dmrg_eig(tt_H,1e-6,'numblocks',1);

%QR dmrg_eig (base case)
%____________________________________________________________________
disp('╔═════════════════════════╗')
disp( "        QR based DMRG eigen test")
disp('╚═════════════════════════╝')

tic;
[x,theta,out]=dmrg_eig_qr(tt_H,1e-6);
toc;

%Ensure that compute dmrg results are consistent by checking norm for
%Ax=lambdaX
fprintf(' \n━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
disp   ('        norm test')
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n')
norm(H*full(x)-theta*full(x))

% %SVD dmrg_eig (base case)
% %____________________________________________________________________
% disp('╔═════════════════════════╗')
% disp( "        SVD based DMRG eigen test")
% disp('╚═════════════════════════╝')
% 
% tic;
% [x,theta,out]=dmrg_eig_svd(tt_H,1e-6);
% toc;
% 
% %Ensure that compute dmrg results are consistent by checking norm for
% %Ax=lambdaX
% fprintf(' \n━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
% disp   ('        norm test')
% fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n')
% norm(H*full(x)-theta*full(x))
% 
% 
% 
% %rSVD dmrg_eig (base case)
% %____________________________________________________________________
% disp('╔═════════════════════════╗')
% disp( "        rSVD based DMRG eigen test")
% disp('╚═════════════════════════╝')
% 
% tic;
% [x,theta,out]=dmrg_eig_rsvd(tt_H,1e-6);
% toc;
% 
% %Ensure that compute dmrg results are consistent by checking norm for
% %Ax=lambdaX
% fprintf(' \n━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
% disp   ('        norm test')
% fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n')
% norm(H*full(x)-theta*full(x))
