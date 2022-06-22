% Author:Christian Camano 
% Implementation of Als Algorithm from page 9 of:
% "The Alternating Linear Scheme for Tensor Optimization in the TT format.
% Sebastian Holtz, Thorsten Rohwedder, et al. 


%examples will be tested on an order 4 tensor train and implemented with
%TT_toolbox. 

%Initialize tensor with dimensions 3x4x5
disp('------------------------------------------------')
A=randn(8,2);
T=reshape(A,[2 2 2 2]);

%construct an order 4 tensor with error threshold .1%
tt_y=tt_tensor(T,.001)
disp('Size of TT')
size(tt_y)

%right to left orthogonalization
  

%Line 2.

tt_x=tt_y;
core(tt_x,1)
tt_x.core
size(tt_x.core)

cr=tt_x.core;
ps=tt_x.ps;

%Line 3:
for n=4:-1:2
   disp('=====================')
   %Line 4
   n
   H=unfold_H(core(tt_x,n));
   [Q,R]=qr(H');
   disp('size Q')
   size(Q')
   tt_x.core( tt_x.ps(n): tt_x.ps(n+1)-1)=reshape(H,1,[]);
   
   %line 5
   V=unfold_V(core(tt_y,n-1))*R.'; %V* R^T
   tt_x.core( tt_x.ps(n-1): tt_x.ps(n)-1)=reshape(V,1,[]);
   disp('loop iteration complete')
end



%test for accuracy
disp('norm test')
norm(tt_y)-norm(tt_w)
disp('element test')
tt_y(2,2,2,1)
tt_w(2,2,2,1)
disp('------------------------------------------------')

function H=unfold_H(core)

  if(ismatrix(core))
      H=core;
    return;
  end
  
  [x y z]=size(core);
  H=reshape(core,[x y*z]);
  disp('size H')
  size(H)
end

function V=unfold_V(core)
  if(ismatrix(core))
    return;
  end
  [x y z]=size(core);
  V=reshape(core,[x*y z]);
end
