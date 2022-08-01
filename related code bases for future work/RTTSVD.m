%generate a d dimensional ones
X=randn(1000,1000);
%Choose example target sketch rank
s=100;
ny = size(X,2);
%generate random matrix ny x s 
P = randn(ny,s);
size(P);
%convert to tensor train form with one node two legs
tt_tensor(P,.0001);

%Create a 4x4x4 tensor with matlab encoding 1-64


A = [1:4; 5:8; 9:12; 13:16];
A(:,:,2) = [17:20; 21:24; 25:28 ; 29:32];
A = cat(3,A,[33:36; 37:40; 41:44; 45:48]);
A = cat(3,A,[49:52; 53:56; 57:60 ; 61:64]);
display('size A:')
size(A)

%Create a 4x4x4 tensor with matlab encoding 2(1-64)
B = 2*[1:4; 5:8; 9:12; 13:16];
B(:,:,2) = 2*[17:20; 21:24; 25:28 ; 29:32];
B = cat(3,B,2*[33:36; 37:40; 41:44; 45:48]); 
B = cat(3,B,2*[49:52; 53:56; 57:60 ; 61:64]);
size(B)

%Create a 4x4x4 tensor with matlab encoding 3(1-64)
C = 3*[1:4; 5:8; 9:12; 13:16];
C(:,:,2) = 3*[17:20; 21:24; 25:28 ; 29:32];
C = cat(3,C,3*[33:36; 37:40; 41:44; 45:48]);   
C = cat(3,C,3*[49:52; 53:56; 57:60 ; 61:64]);
size(C)

%Create a 4x4x4 tensor with matlab encoding 4(1-64)
D = 4*[1:4; 5:8; 9:12; 13:16];
D(:,:,2) = 4*[17:20; 21:24; 25:28 ; 29:32];
D = cat(3,D,4*[33:36; 37:40; 41:44; 45:48]);  
D = cat(3,D,4*[49:52; 53:56; 57:60 ; 61:64]);
size(D)
%Create a 4x4x16 tensor to be represented as a 4 tensor 
E=cat(3,A,cat(3,B,cat(3,C,D)));
size(E)
% call tensor contstructor to make 4x4x4x4 tensor 
T=tensor(E,[4 4 4 4])
size(T)

% T_m=reshape(T,[4 64]) (Tensor version)
T_m=tenmat(T,[1],[2 3 4])
size(T_m)


%compute SVD then assemble SVD network. 
class(T_m)
[U S V]=svd(cast(T_m,"matrix"),'econ');
U_tt=tt_tensor(U);
S_tt=tt_tensor(S);
V_tt=tt_tensor(V);

% now that the matrix is constructed as an order 4 tensor step 2 of the svd
% algorithm for 4 tensors to TT is to matrisize the tensor to grouping
% n2n3n4 on the right and n1 on the left to set it up for svd 




% % convert it to a tensor train network using TT_toolbox
% B_tt=tt_tensor(B)
% % convert it to a tensor using tensor_toolbox
% B_t=tensor(B)
% size(B_t)
% % At this point b_t is our starting order 4 matrix
% reshape(B_t,[64 4 ])