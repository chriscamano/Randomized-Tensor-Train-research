function OTN=randnTN(n)
% OTN=randnTN(n)
% --------------
% Constructs a unit-rank MPO that represents a particular random matrix.
% A unit MPO-rank MPO can always be constructed to represent a particular
% n^d x K random matrix via Lemma 5.1 of "Computing low-rank approximations
% of large-scale matrices with the tensor network randomized SVD".
%
% OTN	=	Tensor Network, MPO that represents a particular random matrix,
%
% n 	=	matrix, n(i,:) contains the dimensions of the ith MPO-tensor.
%
% Reference
% ---------
%
% Computing low-rank approximations of large-scale matrices with the tensor network randomized SVD.
%
% 06-07-2016, Kim Batselier, Wenjian Yu, Luca Daniel, Ngai Wong

r=1;
d=size(n,1);
OTN.n=[ones(d,1) n ones(d,1)];
for i=1:d-1
    OTN.n(i,end)=r;
    OTN.n(i+1,1)=r;
end
for i=1:d
   OTN.core{i}=reshape(randn(OTN.n(i,:)),OTN.n(i,:)); 
end

end
