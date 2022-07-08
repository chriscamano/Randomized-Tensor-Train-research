
%% input matrix construction
n = 10;                                          % number of spins
H = tt_HamHeis(n);                                  % Heisenberg Hamiltonian
L=HamHeis(n);
norm(full(H))-norm(L)
%tt_H = tt_matrix(H,eps,2*ones(1,n),2*ones(1,n))  % convert to TT format
theta=eigs(full(H),1,'smallestreal')
class(H)
%consider dynamic rank adjustment for svd versus rsvd during the actual
%algorithm adjusment .

%% Control dmrg test
tic;
[x,theta,out]=dmrg_eig(H,1e-10,'numblocks',1);
res_qr=norm(full(H)*full(x)-theta*full(x))
toc;
%% QR based DMRG test
disp('╔═════════════════════════╗')
disp( "        QR based DMRG eigen test")
disp('╚═════════════════════════╝')

tic;
[x_qr,theta_qr,out_qr]=dmrg_eig_qr(H,1e-10,'numblocks',1);
toc;

%Ensure that compute dmrg results are consistent by checking norm for Ax=lambdaX
fprintf(' \n━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
disp   ('        norm test')
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n')
res_qr=norm(full(H)*full(x_qr)-theta_qr*full(x_qr))

%% SVD based DMRG test
disp('╔═════════════════════════╗')
disp( "        SVD based DMRG eigen test")
disp('╚═════════════════════════╝')

tic;
[x_svd,theta_svd,out_svd]=dmrg_eig_svd(H,1e-10);
toc;

%Ensure that compute dmrg results are consistent by checking norm for Ax=lambdaX
fprintf(' \n━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
disp   ('        norm test')
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n')
res_svd=norm(full(H)*full(x_svd)-theta_svd*full(x_svd))

%% Randomized SVD Test
disp('╔═════════════════════════╗')
disp( "        rSVD based DMRG eigen test")
disp('╚═════════════════════════╝')

tic;
[x_rsvd,theta_rsvd,out_rsvd]=dmrg_eig_rsvd(H,1e-10);
toc;

%Ensure that compute dmrg results are consistent by checking norm for Ax=lambdaX
fprintf(' \n━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
disp   ('        norm test')
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n')
res_rsvd=norm(full(H)*full(x_rsvd)-theta_rsvd*full(x_rsvd))









%compute average time
% qr_time=0;
% svd_time=0;
% rsvd_time=0;

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
