function R=maxRank(A,n)
% R=maxRank(A,n)
% --------------
% Computes the MPO-rank for given MPO-tensor dimensions n obtained from 
% matrix2mpo.m without computing the result. This function is useful to
% determine feasible block sizes.
%
% R 	=	scalar, MPO-rank,
%
% A 	=	matrix, matrix that needs to be converted into an MPO,
%
% n 	=	matrix, n(i,:) contains the dimensions of the ith MPO-tensor.
%
% Reference
% ---------
%
% Computing low-rank approximations of large-scale matrices with the tensor network randomized SVD.
%
% 06-07-2016, Kim Batselier, Wenjian Yu, Luca Daniel, Ngai Wong

tileX=n(1,1);
tileY=n(1,2);
lrows=zeros(1,prod(n(2:end,2)));

% first keep track how many nonzero block matrices there are and where they
% are
R=0; % the total maximal MPO-rank
rows=cell(1,prod(n(2:end,2)));
for i=1:prod(n(2:end,2)) % i counts the columns
   % determine which blocks on this "column" are nonzero
   I=find(sum(abs(A(:,(i-1)*tileY+1:i*tileY)),2));
   if ~isempty(I)
       I2=find(rem(I,tileX)==0);
       I=fix(I/tileX);
       I=I+1;
       I(I2)=I(I2)-1;
       rows{i}=[I(1);I(find(diff(I))+1)]';
       R=R+length(rows{i});
   end
end
end
