function [B,T]=deparallel(A,varargin)
% Finds matrices B and T such that A=B*T and B has at most as many columns 
% as A and no two columns of B are parallel to each other.

[p,q]=size(A);
if ~isempty(varargin)
    tol=varargin{1};
else
    tol=p*q*eps;
end

% remove zero columns from A
T=zeros(1,q);
nzI=find(sum(abs(A),1) > tol);

% initialize
B(:,1)=A(:,nzI(1));
T(1,nzI(1))=1;
for j=2:length(nzI)
    % check whether column A(:,j) is parallel to any of the columns in B
    notparallel=zeros(1,size(B,2));
    for i=1:size(B,2)
        % find nonzero entries in both columns
        Ia=abs(A(:,nzI(j)))>tol;
        Ib=abs(B(:,i))>tol;
%         Ia=find(abs(A(:,nzI(j)))>tol);
%         Ib=find(abs(B(:,i))>tol);
        if sum(Ia)~=sum(Ib) || sum(abs(Ia-Ib))~=0 
%         if length(Ia)~=length(Ib) || sum(Ia-Ib)~=0 
            notparallel(i)=1;
        else
            ratio=A(Ia,nzI(j))./B(Ib,i);
            scalefactor=ratio(1);
            ratio=ratio/scalefactor; % normalize with respect to first element
%             if abs(sum(ratio)-length(Ia)) > tol
            if abs(sum(ratio)-sum(Ia)) > tol                
                % not parallel
                notparallel(i)=1;
            else
                % parallel
                T(i,nzI(j))=scalefactor;
            end
        end
    end 
    if sum(notparallel)==size(B,2)
        % column A(:,j) was not parallel to any of the columns in B
        % so we add it to B
        B=[B A(:,nzI(j))];
        T(size(B,2),nzI(j))=1;        
    end
end

end