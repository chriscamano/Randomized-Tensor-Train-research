function [U,S,VTN]=svdTN(ATN,tol)
% [U,S,VTN]=svdTN(ATN,tol)
% ------------------------
% Computes the SVD of a k x n^d matrix A in MPO-form ATN, where orthogonal
% U matrix is k x k, S is k x k and orthogonal V is n^d x k and given in
% MPO-form VTN.
%
% U 	=	matrix, orthogonal matrix containing left singular vectors,
%
% S 	=	matrix, diagonal matrix containing singular values,
%
% VTN 	=	Tensor Network, MPO that represents the right singular vectors,
%
% ATN 	=	Tensor Network, MPO that represents original kxn^d matrix A,
%
% tol 	=	scalar, tolerance used in SVD-basd rounding of the MPO.
%
% Reference
% ---------
%
% Computing low-rank approximations of large-scale matrices with the tensor network randomized SVD.
%
% 06-07-2016, Kim Batselier, Wenjian Yu, Luca Daniel, Ngai Wong

[dTN,d]=size(ATN.n);
VTN.core=cell(1,dTN);
VTN.n=ATN.n;
delta=tol/sqrt(d-1);

% if ATN.n(1,end) is too large, then better do a right-to-left sweep to cut
% down the rank of the first core
if ATN.n(1,end) >= 1000
    for i=dTN:-1:2
        [U,S,V]=svd(reshape(ATN.core{i},[ATN.n(i,1),prod(ATN.n(i,2:end))]),'econ');
        s=diag(S);
        r=0;
        while norm(s(end:-1:end-r)) < delta 
            r=r+1;
        end
        r=length(s)-r;
        ATN.core{i}=reshape(V(:,1:r)',[r,ATN.n(i,2:end)]);
        ATN.n(i,:)=[r,ATN.n(i,2:end)];
        ATN.core{i-1}=reshape(reshape(ATN.core{i-1},[prod(ATN.n(i-1,1:end-1)),ATN.n(i-1,end)])*U(:,1:r)*S(1:r,1:r),[ATN.n(i-1,1:end-1),r]);    
        ATN.n(i-1,:)=[ATN.n(i-1,1:end-1),r];
    end        
end


% left to right compression via truncated SVD
for i=1:dTN-1
    [U,S,V]=svd(reshape(ATN.core{i},[prod(ATN.n(i,1:end-1)),ATN.n(i,end)]),'econ');
    s=diag(S);
    r=0;
    while norm(s(end:-1:end-r)) < delta 
        r=r+1;
    end
    r=length(s)-r;
    ATN.core{i}=reshape(U(:,1:r),[ATN.n(i,1),prod(ATN.n(i,2:end-1)),r]);
    ATN.n(i,end)=r;
    ATN.core{i+1}=reshape((S(1:r,1:r)*V(:,1:r)')*reshape(ATN.core{i+1},[ATN.n(i+1,1),prod(ATN.n(i+1,2:end))]),[r,prod(ATN.n(i+1,2:end-1)),ATN.n(i+1,end)]);   
    ATN.n(i+1,1)=r;
end

for i=dTN:-1:2
    [U,S,V]=svd(reshape(ATN.core{i},[ATN.n(i,1),prod(ATN.n(i,2:end))]),'econ');
    s=diag(S);
    r=0;
    while norm(s(end:-1:end-r)) < delta 
        r=r+1;
    end
    r=length(s)-r;
    VTN.core{i}=reshape(V(:,1:r)',[r,VTN.n(i,2:end)]);
    VTN.n(i,:)=[r,VTN.n(i,2:end)];
    ATN.core{i-1}=reshape(reshape(ATN.core{i-1},[prod(ATN.n(i-1,1:end-1)),ATN.n(i-1,end)])*U(:,1:r)*S(1:r,1:r),[ATN.n(i-1,1:end-1),r]);    
    ATN.n(i,:)=[r,VTN.n(i,2:end)];
    ATN.n(i-1,:)=[ATN.n(i-1,1:end-1),r];
    VTN.n(i-1,:)=[VTN.n(i-1,1:end-1),r];
end
[U,S,V]=svd(reshape(ATN.core{1},[ATN.n(1,2),prod(ATN.n(1,3:end))]),'econ');
VTN.core{1}=reshape(V',VTN.n(1,:));
end
