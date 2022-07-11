
% input matrix construction
n = 100;                                          % number of spins
H = tt_HamHeis(n);

%% Heisenberg Hamiltonian
%L=HamHeis(n);
%norm(full(H))-norm(L)
%tt_H = tt_matrix(H,eps,2*ones(1,n),2*ones(1,n))  % convert to TT format
%theta=eigs(L,1,'smallestreal')

% consider dynamic rank adjustment for svd versus rsvd during the actual
% algorithm adjusment .

% % Control dmrg test
% tic;
% [x,theta,out]=dmrg_eig(H,1e-10,'numblocks',1);
% toc;
% fprintf(' \n━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
% disp   ('        norm test')
% fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n')
% 
% res=norm(H*x-theta*x)
% %% QR based DMRG test
% disp('╔═════════════════════════╗')
% disp( "        QR based DMRG eigen test")
% disp('╚═════════════════════════╝')
% 
% tic;
% [x_qr,theta_qr,out_qr]=dmrg_eig_qr(H,1e-10,'numblocks',1);
% toc;
% 
% %Ensure that compute dmrg results are consistent by checking norm for Ax=lambdaX
% fprintf(' \n━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
% disp   ('        norm test')
% fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n')
% res_qr=norm(H*x_qr-theta_qr*x_qr)
% 
% %% SVD based DMRG test
% disp('╔═════════════════════════╗')
% disp( "        SVD based DMRG eigen test")
% disp('╚═════════════════════════╝')
% 
% tic;
% [x_svd,theta_svd,out_svd]=dmrg_eig_svd(H,1e-10);
% toc;
% 
% %Ensure that compute dmrg results are consistent by checking norm for Ax=lambdaX
% fprintf(' \n━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
% disp   ('        norm test')
% fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n')
% res_svd=norm(H*x_svd-theta_svd*x_svd)

% x= tt_rand(H.n, H.d, 5 , -1);
% H=tt_matrix(H*x,[2,2])

%% Randomized SVD Test
disp('╔═════════════════════════╗')
disp( "        rSVD based DMRG eigen test")
disp('╚═════════════════════════╝')

tic;
[x_rsvd,theta_rsvd,out_rsvd]=dmrg_eig_rsvd(H,1e-10,'r',1,'numblocks',1);
toc;


%Ensure that compute dmrg results are consistent by checking norm for Ax=lambdaX
fprintf(' \n━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')
disp   ('        norm test')
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n')

res_rsvd=norm(H*x_rsvd-theta_rsvd*x_rsvd)





%%


% 
% % compute average time
% qr_time=0;
% svd_time=0;
% rsvd_time=0;
% 
% for q=1:100
% disp(q);                                          % number of spins
% H = tt_HamHeis(50);                               % Heisenberg Hamiltonian
% 
% tic;
% [x,theta,out]=dmrg_eig_qr(H,1e-10,'verb',0);
% qr_time=qr_time+toc;  
% tic;
% [x,theta,out]=dmrg_eig_svd(H,1e-10,'verb',0);
% svd_time=svd_time+toc;
% tic;
% [x,theta,out]=dmrg_eig_rsvd(H,1e-10,'verb',0);
% rsvd_time=rsvd_time+toc;
% end
% 
% qr_time=qr_time/100
% svd_time=svd_time/100
% rsvd_time=rsvd_time/100
% 
% disp(qr_time)
% disp(svd_time)
% disp(rsvd_time)
% plot([qr_time svd_time rsvd_time],'ro')
