function [UTN,S,VTN]=TNrSVD(ATN,k,q,tol)
% [UTN,S,VTN]=TNrSVD(ATN,k,q,tol)
% -------------------------------
% Computes a k-rank approximation of a given matrix in MPO-form ATN 
% using the tensor network randomized SVD. Both the orthogonal U and V
% matrices are also returned in MPO-form and S is a diagonal matrix of
% size kxk.
%
% UTN       =   Tensor Network, MPO that represents the orthogonal U matrix in the SVD,
%
% S         =   matrix, diagonal matrix that contains approximations to singular values of A,
%
% VTN       =   Tensor Network, MPO that represents the orthogonal V matrix in the SVD,
%
% ATN       =   Tensor Network, MPO that represents the original matrix A,
%
% k         =   scalar, size of the low-rank approximation (oversampling included),
%
% q 		=	scalar, exponent in subspace iterations (A^T*A)^q,
%
% tol 		=	scalar, tolerance used in the SVD-based MPO-rounding.
%
% Reference
% ---------
%
% Computing low-rank approximations of large-scale matrices with the tensor network randomized SVD.
%
% 06-07-2016, Kim Batselier, Wenjian Yu, Luca Daniel, Ngai Wong

d=size(ATN.n,1); % number of cores in TN
n=[ATN.n(:,3) [k;ones(d-1,1)] ];
OTN=randnTN(n);
YTN=AwithO(ATN,OTN);
[QTN,~]=qrTN(YTN,tol); % qrTN does the rounding
for i=1:q
    QTN=contractab(ATN,QTN,[2 2]); 
    [QTN,~]=qrTN(QTN,tol);
    QTN=contractab(ATN,QTN,[3 2]);
    [QTN,~]=qrTN(QTN,tol);    
end
BTN=QtwithA(QTN,ATN); 
[U,S,VTN]=svdTN(BTN,tol); %svdTN does the rounding
UTN=cmodeprod(QTN,U',3,1); % scaling of first core
end
