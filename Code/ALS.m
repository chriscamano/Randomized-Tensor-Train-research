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

%Line 2.
tt_x=tt_y;
cr=tt_x.core;
ps=tt_x.ps;

core(tt_x)
%Line 3:

disp('╔═════════════════════════╗')
disp( "  Right to left Orthogonalization Algorithm")
disp('╚═════════════════════════╝')

for n=4:-1:2
    %---------------------------------------
    fprintf(' ━━━━━━━━━━━━━▼━━━━━━━━━━━━━\n')
    fprintf(' Current Core index:%d\n',n)
    fprintf(' current core dimensions: %s\n',mat2str(size(core(tt_x,n))))
    fprintf(' ━━━━━━━━━━━━━▲━━━━━━━━━━━━━\n\n')
    %---------------------------------------
    
    %Line 4
    %---------------------------------------
    H = unfold_H(core(tt_x,n));
    [Q,R] = qr(H',0);
    % update value of tensor array with orthogonal basis from QR
    % formula for array manip taken from tt_toolbox docs
    cr(ps(n) : ps(n+1)-1) = reshape(Q',1,[]); 
    %---------------------------------------
    
    %line 5
    %---------------------------------------
    V = unfold_V(core(tt_y,n-1))*R'; %V* R^T
    cr(ps(n-1) : ps(n)-1) = reshape(V,1,[]);
    %---------------------------------------
    
    %---------------------------------------
    fprintf(' loop iteration %d complete \n',n)
    fprintf(' -----------------------------\n\n')
    %---------------------------------------
end

disp('╔═════════════════════════╗')
disp( "           Error analysis step")
disp('╚═════════════════════════╝')


%core(tt_x,4)'*core(tt_x,4)
fprintf(' ━━━━━━━━━━━━━▼━━━━━━━━━━━━━\n')
disp   (' Orthogonality check G_3*G_3')
fprintf(' ━━━━━━━━━━━━━▲━━━━━━━━━━━━━\n\n')

A=unfold_H(core(tt_x,3))'*unfold_H(core(tt_x,3))
heatmap(A,'Colormap',bone)

fprintf(' ━━━━━━━━━━━━━▼━━━━━━━━━━━━━\n')
disp   ('       TT norm test')
fprintf(' ━━━━━━━━━━━━━▲━━━━━━━━━━━━━\n\n')

norm(tt_y-tt_x)

fprintf(' ━━━━━━━━━━━━━▼━━━━━━━━━━━━━\n')
disp('      Matrix norm test')
fprintf(' ━━━━━━━━━━━━━▲━━━━━━━━━━━━━\n\n')

norm(full(tt_y)-full(tt_x))

fprintf(' ━━━━━━━━━━━━━▼━━━━━━━━━━━━━\n')
disp(' Element comparison: (1,1,1,1)')
fprintf(' ━━━━━━━━━━━━━▲━━━━━━━━━━━━━\n\n')

tt_x(1,1,1,1)
tt_y(1,1,1,1)
disp('------------------------------------------------')

function H=unfold_H(core)
  if(ismatrix(core))
      H=core;
    return;
  end
  
  [x y z]=size(core);
  H=reshape(core,[x y*z]);
end

function V=unfold_V(core)
  if(ismatrix(core))
    return;
  end
  
  [x y z]=size(core);
  V=reshape(core,[x*y z]);
end
