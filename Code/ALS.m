% Author:Christian Camano 
% Implementation of Als Algorithm from page 9 of:
% "The Alternating Linear Scheme for Tensor Optimization in the TT format.
% Sebastian Holtz, Thorsten Rohwedder, et al. 


%examples will be tested on an order 4 tensor train and implemented with
%TT_toolbox. 

%Initialize tensor with dimensions 3x4x5
A=randn(8,2)
T=reshape(A,[2 2 2 2])
size(T)

%construct an order 4 tensor with error threshold .1%
tt=tt_tensor(T,.001)
size(tt)

%right to left orthogonalization
cores=core2cell(tt);   


tt_x=tt;
coresy=core2cell(tt);
coresx=core2cell(tt_x);
for n=4:2
   H=unfold_H(coresx{n});
   [H,R]=qr(H.')
   V=unfold_V(coresy{n-1})*R.'
end

function T_H=fold_H(core)

end

function t_V=fold_v(core)

end
function H=unfold_H(core)
  if(ismatrix(core))
    return;
  end
  [x y z]=size(core)
  H=reshape(core,[x y*z])
end
function V=unfold_V(core)
  if(ismatrix(core))
    return;
  end
  [x y z]=size(core)
  V=reshape(core,[x*y z])
end
