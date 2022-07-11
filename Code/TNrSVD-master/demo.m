clear all
load Erdos972
A=Problem.A;
P=symrcm(A); % Symmetric reverse Cuthill-McKee permutation to force entries closer to the main diagonal
Ap=A(P,P);
% MPO cores have dimensions 343x343, 2x2, 2x2, 2x2, 2x2 
n=[7^3,7^3;2*ones(4,2)]; 
d=size(n,1);

% convert the sparse matrix into its exact MPO-representation using
% Algorithm 4.1 of "Computing low-rank approximations of large-scale
% matrices with the tensor network randomized SVD"
v=tic;
ATN=matrix2mpo(A,n);
toc(v)
% check MPO-ranks 
ATN.n
ATN=roundTN(ATN,1e-15);
% ATN=parroundTN(ATN); % round the MPO by parallel column removal
% check MPO-ranks after rounding
ATN.n

% reference singular values
K=10;                  % we want a rank-10 approximation
[Us,Ss,Vs]=svds(A,K);
sref=diag(Ss);
v=tic;
[UTN,S,VTN,err]=qTNrSVD(ATN,2*K,1e-10,1e-2);
toc(v)
shat=diag(S(1:K,1:K));
% check relative error of computed singular values
 norm(sref-shat)/norm(sref)