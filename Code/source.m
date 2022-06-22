% Implemented algorithm described on page 2301 of original paper 
% by oseledetss Algorithm 1. TT-SVD.
% Pre: A d dimensional tensor TN a tensorized version of a matrix A , perscribed accuracy ep 
% Post: The cores of the tensor train approximation B of TN where the cores 
% correspond to the singular value decompostion of Tn. 

function [U,S,VTN]=TT_SVD(tt_a,ep)

[dTN,d]=size(tt_a);
% 1. initiialization: 
delta=(ep/sqrt(d-1))*norm(tt_a)
% 2. Temporary Tensor:


%3.
for i=1:dTN-1
    
end

end