function [UTN,S,VTN,err]=qTNrSVD(ATN,k,roundtol,relerrortol)
% [UTN,S,VTN]=TNrSVD(ATN,k,roundtol,relerrortol)
% ----------------------------------------------
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
% roundtol  =	scalar, tolerance used in the SVD-based MPO-rounding.
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
[QTN,~]=qrTN(YTN,roundtol); % qrTN does the rounding
BTN=QtwithA(QTN,ATN); 
[U,S,VTN]=svdTN(BTN,roundtol); %svdTN does the rounding
sprev=zeros(k,1);
s=diag(S);
q=0;
err(q+1,1)=max( abs(s(1:k/2).^2-sprev(1:k/2).^2));
err(q+1,2)=max( abs(s(1:k/2).^2-sprev(1:k/2).^2)/s(k/2)^2);
err(q+1,3)=max( abs(s(1:k/2).^2-sprev(1:k/2).^2)/s(1)^2);
% UTN=cmodeprod(QTN,U',3,1); % scaling of first core
% test=contractab(ATN,UTN,[3,2]);AtU=contract(test);
% test=cmodeprod(VTN,S,2,1);VS=contract(test)';
% err=norm(AtU(:)-VS(:))/norm(S(:));
while err(q+1,2) > relerrortol
    q=q+1;
    QTN=contractab(ATN,QTN,[2,2]); 
    [QTN,~]=qrTN(QTN,roundtol);
    QTN=contractab(ATN,QTN,[3,2 ]);
    [QTN,~]=qrTN(QTN,roundtol);   
    BTN=QtwithA(QTN,ATN); 
    [U,S,VTN]=svdTN(BTN,roundtol); %svdTN does the rounding
    sprev=s;
    s=diag(S);
    err(q+1,1)=max( abs(s(1:k/2).^2-sprev(1:k/2).^2));
    err(q+1,2)=max( abs(s(1:k/2).^2-sprev(1:k/2).^2)/s(k/2)^2);
    err(q+1,3)=max( abs(s(1:k/2).^2-sprev(1:k/2).^2)/s(1)^2);
end
UTN=cmodeprod(QTN,U',3,1); % scaling of first core
end
