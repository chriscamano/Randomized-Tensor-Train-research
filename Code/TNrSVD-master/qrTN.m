function [QTN,R]=qrTN(ATN,tol)
% [QTN,R]=qrTN(ATN,tol)
% -----------------
% Computes a QR decomposition of a matrix A in the TN format, returning the
% orthogonal Q also in the TN format. The assumption is that A is a n^d x k
% matrix where k << n^d. The second mode of the A matrix needs to be split
% over the TN-cores in such a way that the first core has column dimension k
% and the other cores have dimension 1.
%
% QTN       =	Tensor Network for orthogonal matrix Q,
%
% R	 		=	matrix, upper triangular k x k matrix,
%
% ATN 		=	Tensor Network of a n^d x k matrix.
%
% Reference
% ---------
%
% 11/2016, Kim Batselier, Ngai Wong

[dTN,d]=size(ATN.n);
QTN.core=cell(1,dTN);
QTN.n=ATN.n;
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
    temp=reshape(ATN.core{i},[ATN.n(i,1),prod(ATN.n(i,2:end))]);    
    [Q,S,V]=svd(temp','econ');
    s=diag(S);
    r=0;
    while norm(s(end:-1:end-r)) < delta 
        r=r+1;
    end
    r=length(s)-r;    
    QTN.core{i}=reshape(Q(:,1:r)',[r,QTN.n(i,2:end)]);
    QTN.n(i,:)=[r,QTN.n(i,2:end)];
    ATN.core{i-1}=reshape(reshape(ATN.core{i-1},[prod(ATN.n(i-1,1:end-1)),ATN.n(i-1,end)])*V(:,1:r)*S(1:r,1:r),[ATN.n(i-1,1:end-1),r]);    
    ATN.n(i,:)=[r,QTN.n(i,2:end)];
    ATN.n(i-1,:)=[ATN.n(i-1,1:end-1),r];
    QTN.n(i-1,:)=[QTN.n(i-1,1:end-1),r];
end
temp=permute(reshape(ATN.core{1},ATN.n(1,:)),[d,1:d-1]);
temp=reshape(temp,[ATN.n(1,end)*prod(ATN.n(1,1:end-2)),ATN.n(1,end-1)]);
[Q,R]=qr(temp,0);
QTN.core{1}=permute(reshape(Q,[ATN.n(1,end),ATN.n(1,1:end-2),ATN.n(1,end-1)]),[2:d,1]);


end