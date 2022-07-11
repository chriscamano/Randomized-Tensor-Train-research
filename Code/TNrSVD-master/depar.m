function [B,T]=depar(A,varargin)
% Finds matrices B and T such that A=B*T and B has at most as many columns 
% as A and no two columns of B are parallel to each other.

[p,q]=size(A);
if ~isempty(varargin)
    tol=varargin{1};
else
    tol=p*q*eps;
end

% remove zero columns from A
norms=sqrt(sum(A.^2,1));
qnzI=find(norms > tol)';             % contains the indices of nonzero columns in original A
norms=norms(qnzI);
% A=A(:,qnzI)*diag(1./norms);   % remove zero columns and normalize remaining columns
A=bsxfun(@rdivide,A(:,qnzI),norms);   % remove zero columns and normalize remaining columns

cosines=tril(A'*A,-1); % compute cosines of angles between vectors
[Ir,Ic]=find(abs(cosines-1)<tol);
toremove=unique(Ir);
tokeep=1:size(A,2);
tokeep(toremove)=[];
% B=A(:,tokeep)*diag(norms(tokeep));
B=bsxfun(@times,A(:,tokeep),norms(tokeep));

T=zeros(size(B,2)*q,1);
[LIA,LOCB]=ismember(Ic,tokeep);
lindices=sub2ind([size(B,2),q],LOCB(LIA),qnzI(Ir(LIA)));
T(lindices)=norms(Ir(LIA))./norms(Ic(LIA));
T=reshape(T,[size(B,2),q]);
T(:,qnzI(tokeep))=eye(size(B,2));

end